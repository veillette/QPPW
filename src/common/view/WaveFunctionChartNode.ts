/**
 * WaveFunctionChartNode displays the wave function or probability density
 * for the selected energy state. This is the bottom chart in the One Well screen.
 */

import { Node, Line, Path, Text } from "scenerystack/scenery";
import { Shape } from "scenerystack/kite";
import { NumberProperty } from "scenerystack/axon";
import { Range } from "scenerystack/dot";
import { Orientation } from "scenerystack/phet-core";
import { ChartTransform, ChartRectangle, AxisLine, GridLineSet, TickMarkSet, TickLabelSet } from "scenerystack/bamboo";
import { OneWellModel } from "../../one-well/model/OneWellModel.js";
import { BoundStateResult } from "../model/PotentialFunction.js";
import QuantumConstants from "../model/QuantumConstants.js";
import QPPWColors from "../../QPPWColors.js";

// Chart axis range constant (shared with EnergyChartNode)
const X_AXIS_RANGE_NM = 4; // X-axis extends from -X_AXIS_RANGE_NM to +X_AXIS_RANGE_NM

export class WaveFunctionChartNode extends Node {
  private readonly model: OneWellModel;
  private readonly chartWidth: number;
  private readonly chartHeight: number;
  private readonly chartMargins = { left: 60, right: 20, top: 40, bottom: 50 };

  // Chart bounds in view coordinates
  private readonly plotWidth: number;
  private readonly plotHeight: number;

  // ChartTransform for model-to-view coordinate conversion
  private readonly chartTransform: ChartTransform;

  // View range properties (synchronized with EnergyChartNode X-axis)
  private readonly xMinProperty: NumberProperty;
  private readonly xMaxProperty: NumberProperty;
  private readonly yMinProperty: NumberProperty;
  private readonly yMaxProperty: NumberProperty;

  // Visual elements
  private readonly backgroundRect: ChartRectangle;
  private readonly plotContentNode: Node; // Clipped container for plot content
  private readonly realPartPath: Path;
  private readonly imaginaryPartPath: Path;
  private readonly magnitudePath: Path;
  private readonly probabilityDensityPath: Path;
  private readonly zeroLine: Line;
  private readonly axesNode: Node;
  private yAxisLabel!: Text;

  public constructor(model: OneWellModel, options?: { width?: number; height?: number }) {
    super();

    this.model = model;
    this.chartWidth = options?.width ?? 600;
    this.chartHeight = options?.height ?? 250;

    this.plotWidth = this.chartWidth - this.chartMargins.left - this.chartMargins.right;
    this.plotHeight = this.chartHeight - this.chartMargins.top - this.chartMargins.bottom;

    // Initialize view range (x-axis is fixed, y-axis will be updated based on data)
    this.xMinProperty = new NumberProperty(-X_AXIS_RANGE_NM);
    this.xMaxProperty = new NumberProperty(X_AXIS_RANGE_NM);
    this.yMinProperty = new NumberProperty(-1);
    this.yMaxProperty = new NumberProperty(1);

    // Create ChartTransform for model-to-view coordinate conversion
    this.chartTransform = new ChartTransform({
      viewWidth: this.plotWidth,
      viewHeight: this.plotHeight,
      modelXRange: new Range(this.xMinProperty.value, this.xMaxProperty.value),
      modelYRange: new Range(this.yMinProperty.value, this.yMaxProperty.value),
    });

    // Create background using ChartRectangle
    this.backgroundRect = new ChartRectangle(this.chartTransform, {
      fill: QPPWColors.backgroundColorProperty,
      stroke: QPPWColors.gridLineProperty,
      lineWidth: 1,
    });
    this.backgroundRect.x = this.chartMargins.left;
    this.backgroundRect.y = this.chartMargins.top;
    this.addChild(this.backgroundRect);

    // Create axes
    this.axesNode = this.createAxes();
    this.addChild(this.axesNode);

    // Create a clipped content node for all plot elements
    this.plotContentNode = new Node({
      clipArea: Shape.rectangle(
        this.chartMargins.left,
        this.chartMargins.top,
        this.plotWidth,
        this.plotHeight
      ),
    });
    this.addChild(this.plotContentNode);

    // Create zero line
    this.zeroLine = new Line(0, 0, 0, 0, {
      stroke: QPPWColors.gridLineProperty,
      lineWidth: 1,
      lineDash: [5, 5],
    });
    this.plotContentNode.addChild(this.zeroLine);

    // Create wave function paths
    this.realPartPath = new Path(null, {
      stroke: QPPWColors.wavefunctionRealProperty,
      lineWidth: 2,
      visible: false,
    });
    this.plotContentNode.addChild(this.realPartPath);

    this.imaginaryPartPath = new Path(null, {
      stroke: QPPWColors.wavefunctionImaginaryProperty,
      lineWidth: 2,
      visible: false,
    });
    this.plotContentNode.addChild(this.imaginaryPartPath);

    this.magnitudePath = new Path(null, {
      stroke: "purple",
      lineWidth: 2,
      visible: false,
    });
    this.plotContentNode.addChild(this.magnitudePath);

    this.probabilityDensityPath = new Path(null, {
      stroke: QPPWColors.wavefunctionProbabilityProperty,
      lineWidth: 2,
      fill: "rgba(255, 215, 0, 0.2)", // Semi-transparent fill
    });
    this.plotContentNode.addChild(this.probabilityDensityPath);

    // Link to model properties
    this.linkToModel();

    // Initial update
    this.update();
  }

  /**
   * Creates the axes (X and Y) with labels using bamboo components.
   */
  private createAxes(): Node {
    const axesNode = new Node();

    // Y-axis using bamboo AxisLine
    const yAxis = new AxisLine(this.chartTransform, Orientation.VERTICAL, {
      stroke: QPPWColors.axisProperty,
      lineWidth: 2,
    });
    yAxis.x = this.chartMargins.left;
    yAxis.y = this.chartMargins.top;
    axesNode.addChild(yAxis);

    // X-axis using bamboo AxisLine
    const xAxis = new AxisLine(this.chartTransform, Orientation.HORIZONTAL, {
      stroke: QPPWColors.axisProperty,
      lineWidth: 2,
    });
    xAxis.x = this.chartMargins.left;
    xAxis.y = this.chartMargins.top + this.plotHeight;
    axesNode.addChild(xAxis);

    // X-axis grid lines
    const xGridLines = new GridLineSet(this.chartTransform, Orientation.VERTICAL, 2, {
      stroke: QPPWColors.gridLineProperty,
      lineWidth: 1,
    });
    xGridLines.x = this.chartMargins.left;
    xGridLines.y = this.chartMargins.top;
    axesNode.addChild(xGridLines);

    // X-axis tick marks
    const xTickMarks = new TickMarkSet(this.chartTransform, Orientation.HORIZONTAL, 2, {
      edge: "min",
      extent: 8,
      stroke: QPPWColors.axisProperty,
      lineWidth: 1,
    });
    xTickMarks.x = this.chartMargins.left;
    xTickMarks.y = this.chartMargins.top + this.plotHeight;
    axesNode.addChild(xTickMarks);

    // X-axis tick labels
    const xTickLabels = new TickLabelSet(this.chartTransform, Orientation.HORIZONTAL, 2, {
      edge: "min",
      extent: 10,
      createLabel: (value: number) => new Text(value.toFixed(0), {
        font: "12px sans-serif",
        fill: QPPWColors.labelFillProperty,
      }),
    });
    xTickLabels.x = this.chartMargins.left;
    xTickLabels.y = this.chartMargins.top + this.plotHeight;
    axesNode.addChild(xTickLabels);

    // Y-axis label (will be updated based on display mode)
    this.yAxisLabel = new Text("Probability Density", {
      font: "14px sans-serif",
      fill: QPPWColors.labelFillProperty,
      rotation: -Math.PI / 2,
      centerX: this.chartMargins.left - 40,
      centerY: this.chartHeight / 2,
    });
    axesNode.addChild(this.yAxisLabel);

    // X-axis label
    const xLabel = new Text("Position (nm)", {
      font: "14px sans-serif",
      fill: QPPWColors.labelFillProperty,
      centerX: this.chartWidth / 2,
      centerY: this.chartHeight - 15,
    });
    axesNode.addChild(xLabel);

    return axesNode;
  }

  /**
   * Links chart updates to model property changes.
   */
  private linkToModel(): void {
    // Update when any parameter changes
    this.model.potentialTypeProperty.link(() => this.update());
    this.model.wellWidthProperty.link(() => this.update());
    this.model.wellDepthProperty.link(() => this.update());
    this.model.wellOffsetProperty.link(() => this.update());
    this.model.particleMassProperty.link(() => this.update());
    this.model.selectedEnergyLevelIndexProperty.link(() => this.update());
    this.model.displayModeProperty.link(() => {
      this.updateYAxisLabel();
      this.update();
    });
    this.model.timeProperty.link(() => this.updateTimeEvolution());

    // Update visibility of wave function components
    this.model.showRealPartProperty.link((show: boolean) => {
      this.realPartPath.visible = show && this.model.displayModeProperty.value === "waveFunction";
    });
    this.model.showImaginaryPartProperty.link((show: boolean) => {
      this.imaginaryPartPath.visible = show && this.model.displayModeProperty.value === "waveFunction";
    });
    this.model.showMagnitudeProperty.link((show: boolean) => {
      this.magnitudePath.visible = show && this.model.displayModeProperty.value === "waveFunction";
    });
  }

  /**
   * Updates the Y-axis label based on display mode.
   */
  private updateYAxisLabel(): void {
    const displayMode = this.model.displayModeProperty.value;
    if (displayMode === "probabilityDensity") {
      this.yAxisLabel.string = "Probability Density";
    } else {
      this.yAxisLabel.string = "Wave Function";
    }
  }

  /**
   * Main update method - recalculates and redraws everything.
   */
  private update(): void {
    const boundStates = this.model.getBoundStates();
    if (!boundStates) {
      return;
    }

    const selectedIndex = this.model.selectedEnergyLevelIndexProperty.value;
    if (selectedIndex < 0 || selectedIndex >= boundStates.wavefunctions.length) {
      return;
    }

    this.updateViewRange(boundStates);
    this.updateZeroLine();
    this.updateWaveFunction(boundStates, selectedIndex);
  }

  /**
   * Updates the view range based on the data and updates the ChartTransform.
   * Note: X-axis range is fixed, only Y-axis is updated dynamically.
   */
  private updateViewRange(boundStates: BoundStateResult): void {
    // Calculate Y range based on wave function values
    const selectedIndex = this.model.selectedEnergyLevelIndexProperty.value;
    if (selectedIndex >= 0 && selectedIndex < boundStates.wavefunctions.length) {
      const wavefunction = boundStates.wavefunctions[selectedIndex];
      const maxAbs = Math.max(...wavefunction.map((v) => Math.abs(v)));
      const displayMode = this.model.displayModeProperty.value;

      if (displayMode === "probabilityDensity") {
        this.yMinProperty.value = 0;
        this.yMaxProperty.value = maxAbs * maxAbs * 1.2; // 20% margin
      } else {
        this.yMinProperty.value = -maxAbs * 1.2;
        this.yMaxProperty.value = maxAbs * 1.2;
      }
    }

    // Update ChartTransform with new Y range (X range is fixed)
    this.chartTransform.setModelYRange(new Range(this.yMinProperty.value, this.yMaxProperty.value));
  }

  /**
   * Updates the zero-line reference.
   */
  private updateZeroLine(): void {
    const y = this.dataToViewY(0);
    this.zeroLine.x1 = this.chartMargins.left;
    this.zeroLine.y1 = y;
    this.zeroLine.x2 = this.chartWidth - this.chartMargins.right;
    this.zeroLine.y2 = y;
  }

  /**
   * Updates the wave function or probability density plot.
   */
  private updateWaveFunction(boundStates: BoundStateResult, selectedIndex: number): void {
    const xGrid = boundStates.xGrid;
    const wavefunction = boundStates.wavefunctions[selectedIndex];
    const displayMode = this.model.displayModeProperty.value;
    const time = this.model.timeProperty.value * 1e-15; // Convert fs to seconds
    const energy = boundStates.energies[selectedIndex];

    // Calculate time evolution phase
    const phase = -(energy * time) / QuantumConstants.HBAR;

    if (displayMode === "probabilityDensity") {
      // Plot probability density |ψ|²
      this.plotProbabilityDensity(xGrid, wavefunction);
      this.probabilityDensityPath.visible = true;
      this.realPartPath.visible = false;
      this.imaginaryPartPath.visible = false;
      this.magnitudePath.visible = false;
    } else {
      // Plot wave function components
      this.plotWaveFunctionComponents(xGrid, wavefunction, phase);
      this.probabilityDensityPath.visible = false;
      this.realPartPath.visible = this.model.showRealPartProperty.value;
      this.imaginaryPartPath.visible = this.model.showImaginaryPartProperty.value;
      this.magnitudePath.visible = this.model.showMagnitudeProperty.value;
    }
  }

  /**
   * Plots the probability density.
   */
  private plotProbabilityDensity(xGrid: number[], wavefunction: number[]): void {
    const shape = new Shape();

    // Start at zero on the left
    const x0 = this.dataToViewX(xGrid[0] * QuantumConstants.M_TO_NM);
    const y0 = this.dataToViewY(0);
    shape.moveTo(x0, y0);

    // Plot the curve
    for (let i = 0; i < xGrid.length; i++) {
      const x = this.dataToViewX(xGrid[i] * QuantumConstants.M_TO_NM);
      const probabilityDensity = wavefunction[i] * wavefunction[i];
      const y = this.dataToViewY(probabilityDensity);
      shape.lineTo(x, y);
    }

    // Close at zero on the right
    const xEnd = this.dataToViewX(xGrid[xGrid.length - 1] * QuantumConstants.M_TO_NM);
    shape.lineTo(xEnd, y0);
    shape.close();

    this.probabilityDensityPath.shape = shape;
  }

  /**
   * Plots wave function components (real, imaginary, magnitude).
   */
  private plotWaveFunctionComponents(xGrid: number[], wavefunction: number[], phase: number): void {
    const realShape = new Shape();
    const imagShape = new Shape();
    const magShape = new Shape();

    const cosPhi = Math.cos(phase);
    const sinPhi = Math.sin(phase);

    for (let i = 0; i < xGrid.length; i++) {
      const x = this.dataToViewX(xGrid[i] * QuantumConstants.M_TO_NM);
      // Invert the Y-axis: negate wavefunction values to render upright
      const psi = -wavefunction[i];

      // Apply time evolution: ψ(x,t) = ψ(x) * e^(-iEt/ℏ)
      const realPart = psi * cosPhi;
      const imagPart = -psi * sinPhi;
      const magnitude = Math.abs(psi);

      const yReal = this.dataToViewY(realPart);
      const yImag = this.dataToViewY(imagPart);
      const yMag = this.dataToViewY(magnitude);

      if (i === 0) {
        realShape.moveTo(x, yReal);
        imagShape.moveTo(x, yImag);
        magShape.moveTo(x, yMag);
      } else {
        realShape.lineTo(x, yReal);
        imagShape.lineTo(x, yImag);
        magShape.lineTo(x, yMag);
      }
    }

    this.realPartPath.shape = realShape;
    this.imaginaryPartPath.shape = imagShape;
    this.magnitudePath.shape = magShape;
  }

  /**
   * Updates only the time evolution (for animation).
   */
  private updateTimeEvolution(): void {
    if (this.model.isPlayingProperty.value) {
      const boundStates = this.model.getBoundStates();
      if (boundStates) {
        const selectedIndex = this.model.selectedEnergyLevelIndexProperty.value;
        if (selectedIndex >= 0 && selectedIndex < boundStates.wavefunctions.length) {
          this.updateWaveFunction(boundStates, selectedIndex);
        }
      }
    }
  }

  /**
   * Converts data X coordinate to view X coordinate using ChartTransform.
   */
  private dataToViewX(x: number): number {
    return this.chartMargins.left + this.chartTransform.modelToViewX(x);
  }

  /**
   * Converts data Y coordinate to view Y coordinate using ChartTransform.
   * Note: ChartTransform already handles the Y-axis inversion.
   */
  private dataToViewY(y: number): number {
    return this.chartMargins.top + (this.plotHeight - this.chartTransform.modelToViewY(y));
  }
}
