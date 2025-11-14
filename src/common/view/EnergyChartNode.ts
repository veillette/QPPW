/**
 * EnergyChartNode displays the potential energy landscape and energy levels.
 * This is the top chart in the One Well screen.
 */

import { Node, Rectangle, Line, Path, Text, VBox } from "scenerystack/scenery";
import { Shape } from "scenerystack/kite";
import { NumberProperty } from "scenerystack/axon";
import { Range } from "scenerystack/dot";
import { Orientation } from "scenerystack/phet-core";
import { Checkbox } from "scenerystack/sun";
import { ChartTransform, ChartRectangle, AxisLine, TickMarkSet, TickLabelSet } from "scenerystack/bamboo";
import { OneWellModel } from "../../one-well/model/OneWellModel.js";
import { PotentialType, BoundStateResult } from "../model/PotentialFunction.js";
import QuantumConstants from "../model/QuantumConstants.js";
import QPPWColors from "../../QPPWColors.js";
import { PhetFont } from "scenerystack/scenery-phet";

// Chart axis range constants
const X_AXIS_RANGE_NM = 4; // X-axis extends from -X_AXIS_RANGE_NM to +X_AXIS_RANGE_NM

// Energy axis ranges depend on potential type
// For potentials with V=0 at center (e.g., Harmonic Oscillator): -5 to 15 eV
// For potentials with V=0 at infinity (e.g., Finite Well, Coulomb): -15 to 5 eV
// For Asymmetric Triangle: -5 to 15 eV (special case with shifted potential)
function getEnergyAxisRange(potentialType: PotentialType): { min: number; max: number } {
  switch (potentialType) {
    case PotentialType.HARMONIC_OSCILLATOR:
      // V=0 at center, grows outward
      return { min: -5, max: 15 };
    case PotentialType.ASYMMETRIC_TRIANGLE:
      // Special case: V=10eV at infinity, V=0 at x=0
      return { min: -5, max: 15 };
    case PotentialType.INFINITE_WELL:
    case PotentialType.FINITE_WELL:
    case PotentialType.MORSE:
    case PotentialType.POSCHL_TELLER:
    case PotentialType.ROSEN_MORSE:
    case PotentialType.ECKART:
    case PotentialType.COULOMB_1D:
    case PotentialType.COULOMB_3D:
    case PotentialType.CUSTOM:
    default:
      // V=0 at infinity (wells with negative energy states)
      return { min: -15, max: 5 };
  }
}

export class EnergyChartNode extends Node {
  private readonly model: OneWellModel;
  private readonly chartWidth: number;
  private readonly chartHeight: number;
  private readonly chartMargins = { left: 60, right: 20, top: 40, bottom: 50 };

  // Chart bounds in view coordinates
  private readonly plotWidth: number;
  private readonly plotHeight: number;

  // ChartTransform for model-to-view coordinate conversion
  private readonly chartTransform: ChartTransform;

  // View range properties (can be dynamic for zooming)
  private readonly xMinProperty: NumberProperty;
  private readonly xMaxProperty: NumberProperty;
  private readonly yMinProperty: NumberProperty;
  private readonly yMaxProperty: NumberProperty;

  // Visual elements
  private readonly backgroundRect: ChartRectangle;
  private readonly plotContentNode: Node; // Clipped container for plot content
  private readonly potentialPath: Path;
  private readonly energyLevelNodes: Map<number, Line>;
  private readonly energyLabelNodes: Map<number, Text>;
  private readonly totalEnergyLine: Line;
  private readonly zeroLine: Line;
  private axesNode: Node;
  private readonly legendNode: Node;

  // Hover state
  private hoveredEnergyLevelIndex: number | null = null;

  public constructor(model: OneWellModel, options?: { width?: number; height?: number }) {
    super();

    this.model = model;
    this.chartWidth = options?.width ?? 600;
    this.chartHeight = options?.height ?? 300;

    this.plotWidth = this.chartWidth - this.chartMargins.left - this.chartMargins.right;
    this.plotHeight = this.chartHeight - this.chartMargins.top - this.chartMargins.bottom;

    // Initialize view range with values based on initial potential type
    const initialEnergyRange = getEnergyAxisRange(model.potentialTypeProperty.value);
    this.xMinProperty = new NumberProperty(-X_AXIS_RANGE_NM);
    this.xMaxProperty = new NumberProperty(X_AXIS_RANGE_NM);
    this.yMinProperty = new NumberProperty(initialEnergyRange.min);
    this.yMaxProperty = new NumberProperty(initialEnergyRange.max);

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
      stroke: "red",
      lineWidth: 1,
      lineDash: [5, 5],
    });
    this.plotContentNode.addChild(this.zeroLine);

    // Create potential energy path
    this.potentialPath = new Path(null, {
      stroke: QPPWColors.potentialWellProperty,
      lineWidth: 3,
    });
    this.plotContentNode.addChild(this.potentialPath);

    // Create energy level lines and labels containers
    this.energyLevelNodes = new Map();
    this.energyLabelNodes = new Map();

    // Create total energy line
    this.totalEnergyLine = new Line(0, 0, 0, 0, {
      stroke: QPPWColors.energyLevelProperty,
      lineWidth: 2,
      lineDash: [10, 5],
    });
    this.plotContentNode.addChild(this.totalEnergyLine);

    // Create legend (outside clipped area)
    this.legendNode = this.createLegend();
    this.addChild(this.legendNode);

    // Link to model properties
    this.linkToModel();

    // Initial update
    this.update();
  }

  /**
   * Creates the axes (X and Y) with labels - using bamboo components where possible.
   * Note: GridLineSet from bamboo causes infinite loops, so we use manual grid lines.
   */
  private createAxes(): Node {
    const axesNode = new Node();

    // Manual Y-axis grid lines (GridLineSet causes hang)
    const yMin = this.yMinProperty.value;
    const yMax = this.yMaxProperty.value;
    for (let energy = yMin; energy <= yMax; energy += 5) {
      if (energy !== yMin) {
        const y = this.chartMargins.top + this.chartTransform.modelToViewY(energy);
        const gridLine = new Line(
          this.chartMargins.left, y,
          this.chartMargins.left + this.plotWidth, y, {
          stroke: QPPWColors.gridLineProperty,
          lineWidth: 1,
        });
        axesNode.addChild(gridLine);
      }
    }

    // Manual X-axis grid lines (GridLineSet causes hang)
    for (let pos = -X_AXIS_RANGE_NM; pos <= X_AXIS_RANGE_NM; pos += 2) {
      if (pos !== -X_AXIS_RANGE_NM) {
        const x = this.chartMargins.left + this.chartTransform.modelToViewX(pos);
        const gridLine = new Line(
          x, this.chartMargins.top,
          x, this.chartMargins.top + this.plotHeight, {
          stroke: QPPWColors.gridLineProperty,
          lineWidth: 1,
        });
        axesNode.addChild(gridLine);
      }
    }

    // Y-axis using bamboo AxisLine
    const yAxisNode = new AxisLine(this.chartTransform, Orientation.VERTICAL, {
      stroke: QPPWColors.axisProperty,
      lineWidth: 2,
    });
    yAxisNode.x = this.chartMargins.left;
    yAxisNode.y = this.chartMargins.top;
    axesNode.addChild(yAxisNode);

    // X-axis using bamboo AxisLine
    const xAxisNode = new AxisLine(this.chartTransform, Orientation.HORIZONTAL, {
      stroke: QPPWColors.axisProperty,
      lineWidth: 2,
    });
    xAxisNode.x = this.chartMargins.left;
    xAxisNode.y = this.chartMargins.top + this.plotHeight;
    axesNode.addChild(xAxisNode);

    // Y-axis tick marks using bamboo TickMarkSet
    const yTickMarksNode = new TickMarkSet(this.chartTransform, Orientation.VERTICAL, 5, {
      edge: "min",
      extent: 8,
      stroke: QPPWColors.axisProperty,
      lineWidth: 1,
    });
    yTickMarksNode.x = this.chartMargins.left;
    yTickMarksNode.y = this.chartMargins.top;
    axesNode.addChild(yTickMarksNode);

    // X-axis tick marks using bamboo TickMarkSet
    const xTickMarksNode = new TickMarkSet(this.chartTransform, Orientation.HORIZONTAL, 2, {
      edge: "min",
      extent: 8,
      stroke: QPPWColors.axisProperty,
      lineWidth: 1,
    });
    xTickMarksNode.x = this.chartMargins.left;
    xTickMarksNode.y = this.chartMargins.top + this.plotHeight;
    axesNode.addChild(xTickMarksNode);

    // Y-axis tick labels using bamboo TickLabelSet
    const yTickLabelsNode = new TickLabelSet(this.chartTransform, Orientation.VERTICAL, 5, {
      edge: "min",
      createLabel: (value: number) => new Text(value.toFixed(0), {
        font: new PhetFont(12),
        fill: QPPWColors.labelFillProperty,
      }),
    });
    yTickLabelsNode.x = this.chartMargins.left;
    yTickLabelsNode.y = this.chartMargins.top;
    axesNode.addChild(yTickLabelsNode);

    // X-axis tick labels using bamboo TickLabelSet
    const xTickLabelsNode = new TickLabelSet(this.chartTransform, Orientation.HORIZONTAL, 2, {
      edge: "min",
      createLabel: (value: number) => new Text(value.toFixed(0), {
        font: new PhetFont(12),
        fill: QPPWColors.labelFillProperty,
      }),
    });
    xTickLabelsNode.x = this.chartMargins.left;
    xTickLabelsNode.y = this.chartMargins.top + this.plotHeight;
    axesNode.addChild(xTickLabelsNode);

    // Axis labels
    const yLabelText = new Text("Energy (eV)", {
      font: new PhetFont(14),
      fill: QPPWColors.labelFillProperty,
      rotation: -Math.PI / 2,
      centerX: this.chartMargins.left - 40,
      centerY: this.chartHeight / 2,
    });
    axesNode.addChild(yLabelText);

    const xLabelText = new Text("Position (nm)", {
      font: new PhetFont(14),
      fill: QPPWColors.labelFillProperty,
      centerX: this.chartWidth / 2,
      centerY: this.chartHeight - 15,
    });
    axesNode.addChild(xLabelText);

    return axesNode;
  }

  /**
   * Creates the legend for toggling visibility of elements.
   */
  private createLegend(): Node {
    const legendContentNode = new VBox({
      spacing: 5,
      align: "left",
      children: [
        new Checkbox(this.model.showTotalEnergyProperty, new Text("Total Energy", {
          font: "12px sans-serif",
          fill: QPPWColors.textFillProperty,
        }), {
          boxWidth: 15,
        }),
        new Checkbox(this.model.showPotentialEnergyProperty, new Text("Potential Energy", {
          font: "12px sans-serif",
          fill: QPPWColors.textFillProperty,
        }), {
          boxWidth: 15,
        }),
      ],
    });

    const legendPanelNode = new Rectangle(0, 0, legendContentNode.width + 20, legendContentNode.height + 15, 5, 5, {
      fill: "rgba(255, 255, 255, 0.85)",
      stroke: QPPWColors.gridLineProperty,
      lineWidth: 1,
    });

    const legendNode = new Node({
      children: [legendPanelNode, legendContentNode],
      left: this.chartMargins.left + 10,
      top: this.chartMargins.top + 10,
    });

    legendContentNode.left = 10;
    legendContentNode.top = 8;

    return legendNode;
  }

  /**
   * Links chart updates to model property changes.
   */
  private linkToModel(): void {
    // Update when any parameter changes
    this.model.potentialTypeProperty.link(() => {
      this.updateEnergyAxisRange();
      this.update();
    });
    this.model.wellWidthProperty.link(() => this.update());
    this.model.wellDepthProperty.link(() => this.update());
    this.model.wellOffsetProperty.link(() => this.update());
    this.model.particleMassProperty.link(() => this.update());
    this.model.selectedEnergyLevelIndexProperty.link(() => this.updateSelection());
    this.model.showTotalEnergyProperty.link((show) => {
      this.totalEnergyLine.visible = show;
    });
    this.model.showPotentialEnergyProperty.link((show) => {
      this.potentialPath.visible = show;
    });
  }

  /**
   * Updates the energy axis range based on the current potential type.
   */
  private updateEnergyAxisRange(): void {
    const energyRange = getEnergyAxisRange(this.model.potentialTypeProperty.value);
    this.yMinProperty.value = energyRange.min;
    this.yMaxProperty.value = energyRange.max;

    // Update the chart transform with new Y range
    this.chartTransform.setModelYRange(new Range(energyRange.min, energyRange.max));

    // Recreate axes with new range
    this.removeChild(this.axesNode);
    this.axesNode = this.createAxes();
    this.addChild(this.axesNode);
    this.axesNode.moveToBack();
    this.backgroundRect.moveToBack();
  }

  /**
   * Main update method - recalculates and redraws everything.
   */
  private update(): void {
    const boundStates = this.model.getBoundStates();
    if (!boundStates) {
      return;
    }

    this.updatePotentialCurve(boundStates);
    this.updateEnergyLevels(boundStates);
    this.updateZeroLine();
    this.updateTotalEnergyLine(boundStates);
  }

  /**
   * Updates the potential energy curve.
   */
  private updatePotentialCurve(boundStates: BoundStateResult): void {
    const shape = new Shape();
    const xGrid = boundStates.xGrid;

    // For simplicity, draw the potential based on the type
    const potentialType = this.model.potentialTypeProperty.value;
    const wellWidth = this.model.wellWidthProperty.value;
    const wellDepth = this.model.wellDepthProperty.value;

    // Calculate the center of the xGrid for alignment
    const xCenter = ((xGrid[0] + xGrid[xGrid.length - 1]) / 2) * QuantumConstants.M_TO_NM;

    if (potentialType === PotentialType.INFINITE_WELL) {
      // Draw square well centered at xCenter
      const x1 = this.dataToViewX(xCenter - wellWidth / 2);
      const x2 = this.dataToViewX(xCenter + wellWidth / 2);
      const y1 = this.dataToViewY(0);
      const yTop = this.chartMargins.top;
      const yBottom = this.chartHeight - this.chartMargins.bottom;

      // Left wall
      shape.moveTo(x1, yBottom);
      shape.lineTo(x1, yTop);

      // Top (outside well)
      shape.moveTo(this.chartMargins.left, yTop);
      shape.lineTo(x1, yTop);

      // Bottom (inside well)
      shape.moveTo(x1, y1);
      shape.lineTo(x2, y1);

      // Right wall
      shape.moveTo(x2, yTop);
      shape.lineTo(x2, yBottom);

      // Top right
      shape.moveTo(x2, yTop);
      shape.lineTo(this.chartWidth - this.chartMargins.right, yTop);
    } else if (potentialType === PotentialType.FINITE_WELL) {
      // Draw finite square well centered at xCenter
      // Invert the Y-axis: use +wellDepth instead of -wellDepth to create a valley
      const x1 = this.dataToViewX(xCenter - wellWidth / 2);
      const x2 = this.dataToViewX(xCenter + wellWidth / 2);
      const y0 = this.dataToViewY(0);
      const yDepth = this.dataToViewY(wellDepth);

      shape.moveTo(this.chartMargins.left, y0);
      shape.lineTo(x1, y0);
      shape.lineTo(x1, yDepth);
      shape.lineTo(x2, yDepth);
      shape.lineTo(x2, y0);
      shape.lineTo(this.chartWidth - this.chartMargins.right, y0);
    } else if (potentialType === PotentialType.HARMONIC_OSCILLATOR) {
      // Draw parabola centered at xCenter
      // Invert the Y-axis: negate V to create a valley
      const centerX = xCenter;
      const numPoints = 100;
      let firstPoint = true;

      for (let i = 0; i < numPoints; i++) {
        const x = (xGrid[0] + (xGrid[xGrid.length - 1] - xGrid[0]) * i / (numPoints - 1)) * QuantumConstants.M_TO_NM;
        const dx = x - centerX;
        const k = (2 * wellDepth) / (wellWidth * wellWidth / 4); // Spring constant
        const V = -(0.5 * k * dx * dx); // Negate to invert Y-axis
        const viewX = this.dataToViewX(x);
        const viewY = this.dataToViewY(V);

        if (firstPoint) {
          shape.moveTo(viewX, viewY);
          firstPoint = false;
        } else {
          shape.lineTo(viewX, viewY);
        }
      }
    } else if (potentialType === PotentialType.ASYMMETRIC_TRIANGLE) {
      // Draw asymmetric triangle potential: V(∞) = 10 eV, V(0) = 0, V(a) = 10 eV
      // V(x) = 10 eV for x < 0
      // V(x) = 10 - ba + bx for 0 < x < a (linear ramp from 0 to 10 eV)
      // V(x) = 10 eV for x > a
      const centerX = xCenter;
      const POTENTIAL_AT_INFINITY = 10; // eV

      const y10eV = this.dataToViewY(POTENTIAL_AT_INFINITY);
      const y0eV = this.dataToViewY(0);

      // Left region (x < 0): horizontal line at 10 eV
      shape.moveTo(this.chartMargins.left, y10eV);
      shape.lineTo(this.dataToViewX(centerX - wellWidth / 2), y10eV);

      // Triangle region (0 < x < a): linear ramp from 0 to 10 eV
      const x0 = centerX - wellWidth / 2;
      const xa = centerX + wellWidth / 2;
      shape.lineTo(this.dataToViewX(x0), y0eV);
      shape.lineTo(this.dataToViewX(xa), y10eV);

      // Right region (x > a): horizontal line at 10 eV
      shape.lineTo(this.chartWidth - this.chartMargins.right, y10eV);
    } else {
      // For other potential types, draw a placeholder
      const y0 = this.dataToViewY(0);
      shape.moveTo(this.chartMargins.left, y0);
      shape.lineTo(this.chartWidth - this.chartMargins.right, y0);
    }

    this.potentialPath.shape = shape;
  }

  /**
   * Updates the energy level lines.
   */
  private updateEnergyLevels(boundStates: BoundStateResult): void {
    // Remove old energy level nodes
    this.energyLevelNodes.forEach((line) => {
      this.plotContentNode.removeChild(line);
    });
    this.energyLevelNodes.clear();

    // Remove old energy label nodes
    this.energyLabelNodes.forEach((label) => {
      this.removeChild(label);
    });
    this.energyLabelNodes.clear();

    // Create new energy level lines
    const energies = boundStates.energies.map((e) => e * QuantumConstants.JOULES_TO_EV);
    const x1 = this.chartMargins.left;
    const x2 = this.chartWidth - this.chartMargins.right;

    energies.forEach((energy, index) => {
      const y = this.dataToViewY(energy);
      const isSelected = index === this.model.selectedEnergyLevelIndexProperty.value;
      const isHovered = index === this.hoveredEnergyLevelIndex;

      const line = new Line(x1, y, x2, y, {
        stroke: isSelected ? QPPWColors.energyLevelSelectedProperty : QPPWColors.energyLevelProperty,
        lineWidth: isSelected ? 3 : (isHovered ? 2 : 1),
        cursor: "pointer",
        opacity: isHovered ? 1 : 0.7,
      });

      // Add click handler
      line.addInputListener({
        down: () => {
          this.model.selectedEnergyLevelIndexProperty.value = index;
        },
        enter: () => {
          this.hoveredEnergyLevelIndex = index;
          this.updateEnergyLevels(boundStates);
        },
        exit: () => {
          this.hoveredEnergyLevelIndex = null;
          this.updateEnergyLevels(boundStates);
        },
      });

      this.energyLevelNodes.set(index, line);
      this.plotContentNode.addChild(line); // Add to clipped content

      // Add energy label (outside clipped area for visibility)
      const label = new Text(`E${index + 1} = ${energy.toFixed(3)} eV`, {
        font: "10px sans-serif",
        fill: QPPWColors.labelFillProperty,
        left: x2 + 5,
        centerY: y,
      });
      this.energyLabelNodes.set(index, label);
      this.addChild(label);
    });

    // Ensure total energy line and legend are on top
    this.totalEnergyLine.moveToFront();
    this.legendNode.moveToFront();
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
   * Updates the total energy line for the selected state.
   */
  private updateTotalEnergyLine(boundStates: BoundStateResult): void {
    const selectedIndex = this.model.selectedEnergyLevelIndexProperty.value;
    if (selectedIndex >= 0 && selectedIndex < boundStates.energies.length) {
      const energy = boundStates.energies[selectedIndex] * QuantumConstants.JOULES_TO_EV;
      const y = this.dataToViewY(energy);
      this.totalEnergyLine.x1 = this.chartMargins.left;
      this.totalEnergyLine.y1 = y;
      this.totalEnergyLine.x2 = this.chartWidth - this.chartMargins.right;
      this.totalEnergyLine.y2 = y;
      this.totalEnergyLine.visible = this.model.showTotalEnergyProperty.value;
    } else {
      this.totalEnergyLine.visible = false;
    }
  }

  /**
   * Updates the selection highlighting.
   */
  private updateSelection(): void {
    const boundStates = this.model.getBoundStates();
    if (boundStates) {
      this.updateEnergyLevels(boundStates);
      this.updateTotalEnergyLine(boundStates);
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
   * ChartTransform handles the Y-axis inversion (higher model Y = lower view Y).
   */
  private dataToViewY(y: number): number {
    return this.chartMargins.top + this.chartTransform.modelToViewY(y);
  }
}
