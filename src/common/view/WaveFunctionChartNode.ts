/**
 * WaveFunctionChartNode displays the wave function or probability density
 * for the selected energy state. This is the bottom chart in the One Well screen.
 */

import { Node, Line, Path, Text, Rectangle } from "scenerystack/scenery";
import { Shape } from "scenerystack/kite";
import { NumberProperty } from "scenerystack/axon";
import { Range } from "scenerystack/dot";
import { Orientation } from "scenerystack/phet-core";
import {
  ChartTransform,
  ChartRectangle,
  AxisLine,
  TickMarkSet,
  TickLabelSet,
} from "scenerystack/bamboo";
import { OneWellModel } from "../../one-well/model/OneWellModel.js";
import { TwoWellsModel } from "../../two-wells/model/TwoWellsModel.js";
import { BoundStateResult } from "../model/PotentialFunction.js";
import QuantumConstants from "../model/QuantumConstants.js";
import QPPWColors from "../../QPPWColors.js";
import { PhetFont } from "scenerystack/scenery-phet";
import { SuperpositionType } from "../model/SuperpositionType.js";
import stringManager from "../../i18n/StringManager.js";
import QPPWPreferences from "../../QPPWPreferences.js";

// Chart axis range constant (shared with EnergyChartNode)
const X_AXIS_RANGE_NM = 4; // X-axis extends from -X_AXIS_RANGE_NM to +X_AXIS_RANGE_NM

export class WaveFunctionChartNode extends Node {
  private readonly model:
    | OneWellModel
    | TwoWellsModel
    | import("../../many-wells/model/ManyWellsModel.js").ManyWellsModel;
  private readonly chartWidth: number;
  private readonly chartHeight: number;
  private readonly chartMargins = { left: 60, right: 20, top: 10, bottom: 40 };

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
  private readonly phaseColorNode: Node; // Container for phase-colored visualization
  private phaseColorStrips: Rectangle[]; // Pool of rectangles for phase visualization
  private readonly zeroLine: Line;
  private readonly axesNode: Node;
  private yAxisLabel!: Text;
  private readonly stateLabelNode: Text; // Label showing which wavefunction is displayed

  // Guard flag to prevent reentry during updates
  private isUpdating: boolean = false;
  // Flag to indicate if an update was requested while another update was in progress
  private updatePending: boolean = false;

  public constructor(
    model:
      | OneWellModel
      | TwoWellsModel
      | import("../../many-wells/model/ManyWellsModel.js").ManyWellsModel,
    options?: { width?: number; height?: number },
  ) {
    super();

    this.model = model;
    this.chartWidth = options?.width ?? 600;
    this.chartHeight = options?.height ?? 140;

    this.plotWidth =
      this.chartWidth - this.chartMargins.left - this.chartMargins.right;
    this.plotHeight =
      this.chartHeight - this.chartMargins.top - this.chartMargins.bottom;

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
      stroke: null, // Remove border to avoid line appearing below energy chart
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
        this.plotHeight,
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
      stroke: QPPWColors.wavefunctionMagnitudeProperty,
      lineWidth: 2,
      visible: false,
    });
    this.plotContentNode.addChild(this.magnitudePath);

    this.probabilityDensityPath = new Path(null, {
      stroke: QPPWColors.wavefunctionProbabilityProperty,
      lineWidth: 2,
      fill: QPPWColors.wavefunctionProbabilityFillProperty, // Semi-transparent fill
    });
    this.plotContentNode.addChild(this.probabilityDensityPath);

    // Create phase color visualization node
    this.phaseColorNode = new Node({
      visible: false,
    });
    this.plotContentNode.addChild(this.phaseColorNode);

    // Initialize rectangle pool for phase color visualization (will be populated on first use)
    this.phaseColorStrips = [];

    // Create state label in upper right corner (outside clipped area)
    this.stateLabelNode = new Text("", {
      font: new PhetFont({ size: 16, style: "italic" }),
      fill: QPPWColors.labelFillProperty,
      right: this.chartWidth - this.chartMargins.right - 55,
      top: this.chartMargins.top + 5,
    });
    this.addChild(this.stateLabelNode);

    // Ensure axes are on top of the clipped plot content
    this.axesNode.moveToFront();

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

    // Manual X-axis grid lines (GridLineSet causes hang)
    for (let pos = -X_AXIS_RANGE_NM; pos <= X_AXIS_RANGE_NM; pos += 2) {
      if (pos !== -X_AXIS_RANGE_NM) {
        const x =
          this.chartMargins.left + this.chartTransform.modelToViewX(pos);
        const gridLine = new Line(
          x,
          this.chartMargins.top,
          x,
          this.chartMargins.top + this.plotHeight,
          {
            stroke: QPPWColors.gridLineProperty,
            lineWidth: 1,
          },
        );
        axesNode.addChild(gridLine);
      }
    }

    // Y-axis at left edge using bamboo AxisLine (at model x=-4nm)
    const yAxisLeftNode = new AxisLine(
      this.chartTransform,
      Orientation.VERTICAL,
      {
        stroke: QPPWColors.axisProperty,
        lineWidth: 2,
        value: this.xMinProperty.value,
      },
    );
    yAxisLeftNode.x = this.chartMargins.left;
    yAxisLeftNode.y = this.chartMargins.top;
    axesNode.addChild(yAxisLeftNode);

    // Y-axis at origin using bamboo AxisLine (at model x=0)
    const yAxisNode = new AxisLine(this.chartTransform, Orientation.VERTICAL, {
      stroke: QPPWColors.axisProperty,
      lineWidth: 2,
      value: 0,
      opacity: 0.3,
    });
    yAxisNode.x = this.chartMargins.left;
    yAxisNode.y = this.chartMargins.top;
    axesNode.addChild(yAxisNode);

    // X-axis using bamboo AxisLine (at model y=0)
    const xAxisNode = new AxisLine(
      this.chartTransform,
      Orientation.HORIZONTAL,
      {
        stroke: QPPWColors.axisProperty,
        lineWidth: 2,
        value: 0,
      },
    );
    xAxisNode.x = this.chartMargins.left;
    xAxisNode.y = this.chartMargins.top;
    axesNode.addChild(xAxisNode);

    // X-axis tick marks using bamboo TickMarkSet
    const xTickMarksNode = new TickMarkSet(
      this.chartTransform,
      Orientation.HORIZONTAL,
      2,
      {
        edge: "max",
        extent: 8,
        stroke: QPPWColors.axisProperty,
        lineWidth: 1,
      },
    );
    xTickMarksNode.x = this.chartMargins.left;
    xTickMarksNode.y = this.chartMargins.top + this.plotHeight;
    axesNode.addChild(xTickMarksNode);

    // X-axis tick labels using bamboo TickLabelSet
    const xTickLabelsNode = new TickLabelSet(
      this.chartTransform,
      Orientation.HORIZONTAL,
      2,
      {
        edge: "max",
        createLabel: (value: number) =>
          new Text(value.toFixed(0), {
            font: new PhetFont(12),
            fill: QPPWColors.labelFillProperty,
          }),
      },
    );
    xTickLabelsNode.x = this.chartMargins.left;
    xTickLabelsNode.y = this.chartMargins.top + this.plotHeight;
    axesNode.addChild(xTickLabelsNode);

    // Y-axis label (will be updated based on display mode)
    this.yAxisLabel = new Text(
      stringManager.probabilityDensityAxisStringProperty,
      {
        font: new PhetFont(14),
        fill: QPPWColors.labelFillProperty,
        rotation: -Math.PI / 2,
        centerX: this.chartMargins.left - 40,
        centerY: this.chartHeight / 2,
      },
    );
    axesNode.addChild(this.yAxisLabel);

    // X-axis label
    const xLabelText = new Text(stringManager.positionNmAxisStringProperty, {
      font: new PhetFont(14),
      fill: QPPWColors.labelFillProperty,
      centerX: this.chartWidth / 2,
      centerY: this.chartHeight - 15,
    });
    axesNode.addChild(xLabelText);

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
    if ("wellOffsetProperty" in this.model) {
      this.model.wellOffsetProperty.link(() => this.update());
    }
    this.model.particleMassProperty.link(() => this.update());

    // Link to wellSeparationProperty if available (TwoWellsModel only)
    if ("wellSeparationProperty" in this.model) {
      (this.model as TwoWellsModel).wellSeparationProperty.link(() =>
        this.update(),
      );
    }

    // Update when grid points preference changes (affects wavefunction resolution)
    QPPWPreferences.gridPointsProperty.link(() => {
      // Force model to invalidate its cache by notifying it of the change
      // The BaseModel listener will handle cache invalidation
      this.update();
    });

    this.model.selectedEnergyLevelIndexProperty.link(() => {
      this.updateStateLabel();
      this.update();
    });
    this.model.displayModeProperty.link(() => {
      this.updateYAxisLabel();
      this.updateStateLabel();
      this.update();
    });
    this.model.timeProperty.link(() => this.updateTimeEvolution());

    // Update when superposition configuration changes
    this.model.superpositionTypeProperty.link(() => {
      this.updateStateLabel();
      this.update();
    });
    this.model.superpositionConfigProperty.link(() => {
      this.update();
    });

    // Update visibility of wave function components
    this.model.showRealPartProperty.link((show: boolean) => {
      this.realPartPath.visible =
        show && this.model.displayModeProperty.value === "waveFunction";
    });
    this.model.showImaginaryPartProperty.link((show: boolean) => {
      this.imaginaryPartPath.visible =
        show && this.model.displayModeProperty.value === "waveFunction";
    });
    this.model.showMagnitudeProperty.link((show: boolean) => {
      this.magnitudePath.visible =
        show && this.model.displayModeProperty.value === "waveFunction";
    });
  }

  /**
   * Updates the Y-axis label based on display mode.
   */
  private updateYAxisLabel(): void {
    const displayMode = this.model.displayModeProperty.value;
    if (displayMode === "probabilityDensity") {
      this.yAxisLabel.string = "Probability Density";
    } else if (displayMode === "phaseColor") {
      this.yAxisLabel.string = "Wave Function Magnitude";
    } else {
      this.yAxisLabel.string = "Wave Function";
    }
  }

  /**
   * Updates the state label showing which wavefunction is displayed.
   */
  private updateStateLabel(): void {
    const displayMode = this.model.displayModeProperty.value;
    const superpositionType = this.model.superpositionTypeProperty.value;

    // Check if we're in a superposition state (not PSI_K which is single eigenstate)
    const isSuperposition = superpositionType !== SuperpositionType.SINGLE;

    if (isSuperposition) {
      // Display superposition label
      let label = "";
      switch (superpositionType) {
        case SuperpositionType.PSI_I_PSI_J:
          label = "ψ₀+ψ₁";
          break;
        case SuperpositionType.LOCALIZED_NARROW:
          label = "Localized";
          break;
        case SuperpositionType.LOCALIZED_WIDE:
          label = "Localized";
          break;
        case SuperpositionType.COHERENT:
          label = "Coherent";
          break;
        case SuperpositionType.CUSTOM:
          label = "Custom";
          break;
        default:
          label = "Superposition";
      }

      if (displayMode === "probabilityDensity") {
        this.stateLabelNode.string = `|${label}|²`;
      } else if (displayMode === "phaseColor") {
        this.stateLabelNode.string = `|${label}|`;
      } else {
        this.stateLabelNode.string = label;
      }
    } else {
      // Display single eigenstate label
      const selectedIndex = this.model.selectedEnergyLevelIndexProperty.value;

      if (selectedIndex < 0) {
        this.stateLabelNode.string = "";
        return;
      }

      // Convert index to subscript (n starts from 0 for display)
      const n = selectedIndex;
      const subscript = this.toSubscript(n);

      if (displayMode === "probabilityDensity") {
        // Probability density: |ψ_n|²
        this.stateLabelNode.string = `|ψ${subscript}|²`;
      } else if (displayMode === "phaseColor") {
        // Magnitude: |ψ_n|
        this.stateLabelNode.string = `|ψ${subscript}|`;
      } else {
        // Wavefunction: ψ_n
        this.stateLabelNode.string = `ψ${subscript}`;
      }
    }
  }

  /**
   * Converts a number to Unicode subscript characters.
   */
  private toSubscript(num: number): string {
    const subscriptDigits = ["₀", "₁", "₂", "₃", "₄", "₅", "₆", "₇", "₈", "₉"];
    return num
      .toString()
      .split("")
      .map((digit) => subscriptDigits[parseInt(digit)])
      .join("");
  }

  /**
   * Main update method - recalculates and redraws everything.
   * Called automatically when model properties change, but can also be called explicitly (e.g., during reset).
   */
  public update(): void {
    // Prevent reentry - if an update is in progress, mark that another update is pending
    if (this.isUpdating) {
      this.updatePending = true;
      return;
    }

    this.isUpdating = true;
    try {
      // Keep updating until no more updates are pending
      do {
        this.updatePending = false;

        const boundStates = this.model.getBoundStates();
        if (!boundStates) {
          return;
        }

        // Check if we're displaying a superposition or a single eigenstate
        const superpositionType = this.model.superpositionTypeProperty.value;
        const isSuperposition = superpositionType !== SuperpositionType.SINGLE;

        if (isSuperposition && "superpositionConfigProperty" in this.model) {
          // Display superposition wavefunction with proper time evolution
          this.updateViewRangeForSuperpositionFromModel();
          this.updateZeroLine();
          this.updateSuperpositionWavefunction();
          this.updateStateLabel();
        } else {
          // Display single eigenstate
          const selectedIndex =
            this.model.selectedEnergyLevelIndexProperty.value;
          if (
            selectedIndex < 0 ||
            selectedIndex >= boundStates.wavefunctions.length
          ) {
            return;
          }

          this.updateViewRange(boundStates, selectedIndex);
          this.updateZeroLine();
          this.updateWaveFunction(boundStates, selectedIndex);
          this.updateStateLabel();
        }
      } while (this.updatePending);
    } finally {
      this.isUpdating = false;
      this.updatePending = false;
    }
  }

  /**
   * Updates the view range based on the data and updates the ChartTransform.
   * Note: X-axis range is fixed, only Y-axis is updated dynamically.
   */
  private updateViewRange(
    boundStates: BoundStateResult,
    selectedIndex: number,
  ): void {
    // Calculate Y range based on wave function values
    if (
      selectedIndex >= 0 &&
      selectedIndex < boundStates.wavefunctions.length
    ) {
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
    this.chartTransform.setModelYRange(
      new Range(this.yMinProperty.value, this.yMaxProperty.value),
    );
  }

  /**
   * Updates the view range for superposition wavefunctions by computing from model.
   */
  private updateViewRangeForSuperpositionFromModel(): void {
    const boundStates = this.model.getBoundStates();
    if (!boundStates || !("superpositionConfigProperty" in this.model)) {
      return;
    }

    // Guard against empty wavefunctions
    if (boundStates.wavefunctions.length === 0) {
      console.warn("No bound states available for superposition view range");
      return;
    }

    const config = (this.model as OneWellModel).superpositionConfigProperty
      .value;
    const numPoints = boundStates.xGrid.length;

    // Compute the current superposition wavefunction magnitude
    let maxAbs = 0;
    for (let i = 0; i < numPoints; i++) {
      let real = 0;
      let imag = 0;

      for (let n = 0; n < config.amplitudes.length; n++) {
        const amplitude = config.amplitudes[n];
        if (amplitude === 0 || n >= boundStates.wavefunctions.length) {
          continue;
        }

        const eigenfunction = boundStates.wavefunctions[n];
        const phase = config.phases[n];

        real += amplitude * Math.cos(phase) * eigenfunction[i];
        imag += amplitude * Math.sin(phase) * eigenfunction[i];
      }

      const magnitude = Math.sqrt(real * real + imag * imag);
      maxAbs = Math.max(maxAbs, magnitude);
    }

    // Guard against zero maxAbs (no valid states for superposition)
    if (maxAbs === 0) {
      console.warn("No valid states for superposition, using default range");
      maxAbs = 1; // Use default range
    }

    const displayMode = this.model.displayModeProperty.value;

    if (displayMode === "probabilityDensity") {
      this.yMinProperty.value = 0;
      this.yMaxProperty.value = maxAbs * maxAbs * 1.2; // 20% margin
    } else {
      this.yMinProperty.value = -maxAbs * 1.2;
      this.yMaxProperty.value = maxAbs * 1.2;
    }

    // Update ChartTransform with new Y range (X range is fixed)
    this.chartTransform.setModelYRange(
      new Range(this.yMinProperty.value, this.yMaxProperty.value),
    );
  }

  /**
   * Updates the wavefunction display for a superposition state with proper time evolution.
   * Each eigenstate component evolves with its own energy: ψ(x,t) = Σ c_n * ψ_n(x) * e^(-iE_n*t/ℏ)
   */
  private updateSuperpositionWavefunction(): void {
    const boundStates = this.model.getBoundStates();
    if (!boundStates || !("superpositionConfigProperty" in this.model)) {
      return;
    }

    const config = (this.model as OneWellModel).superpositionConfigProperty
      .value;
    const time = this.model.timeProperty.value * 1e-15; // Convert fs to seconds
    const displayMode = this.model.displayModeProperty.value;
    const numPoints = boundStates.xGrid.length;

    // Compute time-evolved superposition: ψ(x,t) = Σ c_n * e^(iφ_n) * ψ_n(x) * e^(-iE_n*t/ℏ)
    const realPart = new Array(numPoints).fill(0);
    const imagPart = new Array(numPoints).fill(0);

    for (let n = 0; n < config.amplitudes.length; n++) {
      const amplitude = config.amplitudes[n];
      const initialPhase = config.phases[n];

      if (amplitude === 0 || n >= boundStates.wavefunctions.length) {
        continue;
      }

      const eigenfunction = boundStates.wavefunctions[n];
      const energy = boundStates.energies[n];

      // Time evolution phase for this eigenstate: -E_n*t/ℏ
      const timePhase = -(energy * time) / QuantumConstants.HBAR;

      // Total phase: initial phase + time evolution phase
      const totalPhase = initialPhase + timePhase;

      // Complex coefficient: c_n * e^(i*totalPhase) = c_n * (cos(totalPhase) + i*sin(totalPhase))
      const realCoeff = amplitude * Math.cos(totalPhase);
      const imagCoeff = amplitude * Math.sin(totalPhase);

      // Add contribution to superposition
      for (let i = 0; i < numPoints; i++) {
        realPart[i] += realCoeff * eigenfunction[i];
        imagPart[i] += imagCoeff * eigenfunction[i];
      }
    }

    // Now display based on display mode
    if (displayMode === "probabilityDensity") {
      // |ψ|² = Real² + Imag²
      const probabilityDensity = new Array(numPoints);
      for (let i = 0; i < numPoints; i++) {
        probabilityDensity[i] =
          realPart[i] * realPart[i] + imagPart[i] * imagPart[i];
      }
      this.plotProbabilityDensityFromArray(
        boundStates.xGrid,
        probabilityDensity,
      );
      this.probabilityDensityPath.visible = true;
      this.realPartPath.visible = false;
      this.imaginaryPartPath.visible = false;
      this.magnitudePath.visible = false;
      this.phaseColorNode.visible = false;
    } else if (displayMode === "phaseColor") {
      // Magnitude and phase for coloring
      this.plotPhaseColoredSuperposition(boundStates.xGrid, realPart, imagPart);
      this.probabilityDensityPath.visible = false;
      this.realPartPath.visible = false;
      this.imaginaryPartPath.visible = false;
      this.magnitudePath.visible = false;
      this.phaseColorNode.visible = true;
    } else {
      // Display real, imaginary, and magnitude components
      const magnitude = new Array(numPoints);
      for (let i = 0; i < numPoints; i++) {
        magnitude[i] = Math.sqrt(
          realPart[i] * realPart[i] + imagPart[i] * imagPart[i],
        );
      }
      this.plotSuperpositionComponents(
        boundStates.xGrid,
        realPart,
        imagPart,
        magnitude,
      );
      this.probabilityDensityPath.visible = false;
      this.realPartPath.visible = this.model.showRealPartProperty.value;
      this.imaginaryPartPath.visible =
        this.model.showImaginaryPartProperty.value;
      this.magnitudePath.visible = this.model.showMagnitudeProperty.value;
      this.phaseColorNode.visible = false;
    }
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
  private updateWaveFunction(
    boundStates: BoundStateResult,
    selectedIndex: number,
  ): void {
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
      this.phaseColorNode.visible = false;
    } else if (displayMode === "phaseColor") {
      // Plot magnitude with phase-colored fill
      this.plotPhaseColoredWavefunction(xGrid, wavefunction, phase);
      this.probabilityDensityPath.visible = false;
      this.realPartPath.visible = false;
      this.imaginaryPartPath.visible = false;
      this.magnitudePath.visible = false;
      this.phaseColorNode.visible = true;
    } else {
      // Plot wave function components
      this.plotWaveFunctionComponents(xGrid, wavefunction, phase);
      this.probabilityDensityPath.visible = false;
      this.realPartPath.visible = this.model.showRealPartProperty.value;
      this.imaginaryPartPath.visible =
        this.model.showImaginaryPartProperty.value;
      this.magnitudePath.visible = this.model.showMagnitudeProperty.value;
      this.phaseColorNode.visible = false;
    }
  }

  /**
   * Plots the probability density with smooth curves.
   */
  private plotProbabilityDensity(
    xGrid: number[],
    wavefunction: number[],
  ): void {
    const shape = new Shape();

    // Start at zero on the left
    const x0 = this.dataToViewX(xGrid[0] * QuantumConstants.M_TO_NM);
    const y0 = this.dataToViewY(0);
    shape.moveTo(x0, y0);

    // Build points array
    const points: { x: number; y: number }[] = [];
    for (let i = 0; i < xGrid.length; i++) {
      const x = this.dataToViewX(xGrid[i] * QuantumConstants.M_TO_NM);
      const probabilityDensity = wavefunction[i] * wavefunction[i];
      const y = this.dataToViewY(probabilityDensity);
      points.push({ x, y });
    }

    // Draw smooth curve using quadratic bezier curves
    if (points.length > 0) {
      shape.lineTo(points[0].x, points[0].y);

      for (let i = 0; i < points.length - 1; i++) {
        const p0 = points[i];
        const p1 = points[i + 1];

        // Control point is midpoint for simple smoothing
        const cpX = (p0.x + p1.x) / 2;
        const cpY = (p0.y + p1.y) / 2;

        shape.quadraticCurveTo(cpX, cpY, p1.x, p1.y);
      }
    }

    // Close at zero on the right
    const xEnd = this.dataToViewX(
      xGrid[xGrid.length - 1] * QuantumConstants.M_TO_NM,
    );
    shape.lineTo(xEnd, y0);
    shape.close();

    this.probabilityDensityPath.shape = shape;
  }

  /**
   * Plots wave function components (real, imaginary, magnitude) with smooth curves.
   */
  private plotWaveFunctionComponents(
    xGrid: number[],
    wavefunction: number[],
    phase: number,
  ): void {
    const realShape = new Shape();
    const imagShape = new Shape();
    const magShape = new Shape();

    const cosPhi = Math.cos(phase);
    const sinPhi = Math.sin(phase);

    // Build points arrays
    const realPoints: { x: number; y: number }[] = [];
    const imagPoints: { x: number; y: number }[] = [];
    const magPoints: { x: number; y: number }[] = [];

    for (let i = 0; i < xGrid.length; i++) {
      const x = this.dataToViewX(xGrid[i] * QuantumConstants.M_TO_NM);
      const psi = wavefunction[i];

      // Apply time evolution: ψ(x,t) = ψ(x) * e^(-iEt/ℏ)
      const realPart = psi * cosPhi;
      const imagPart = -psi * sinPhi;
      const magnitude = Math.abs(psi);

      const yReal = this.dataToViewY(realPart);
      const yImag = this.dataToViewY(imagPart);
      const yMag = this.dataToViewY(magnitude);

      realPoints.push({ x, y: yReal });
      imagPoints.push({ x, y: yImag });
      magPoints.push({ x, y: yMag });
    }

    // Draw smooth curves using quadratic bezier curves
    const drawSmoothCurve = (
      shape: Shape,
      points: { x: number; y: number }[],
    ) => {
      if (points.length > 0) {
        shape.moveTo(points[0].x, points[0].y);

        for (let i = 0; i < points.length - 1; i++) {
          const p0 = points[i];
          const p1 = points[i + 1];

          // Control point is midpoint for simple smoothing
          const cpX = (p0.x + p1.x) / 2;
          const cpY = (p0.y + p1.y) / 2;

          shape.quadraticCurveTo(cpX, cpY, p1.x, p1.y);
        }
      }
    };

    drawSmoothCurve(realShape, realPoints);
    drawSmoothCurve(imagShape, imagPoints);
    drawSmoothCurve(magShape, magPoints);

    this.realPartPath.shape = realShape;
    this.imaginaryPartPath.shape = imagShape;
    this.magnitudePath.shape = magShape;
  }

  /**
   * Plots the wavefunction magnitude with rainbow coloring based on phase.
   * This creates a visualization where vertical strips are colored according to the local phase.
   * Reuses rectangle nodes from a pool to avoid creating/destroying nodes each frame.
   */
  private plotPhaseColoredWavefunction(
    xGrid: number[],
    wavefunction: number[],
    globalPhase: number,
  ): void {
    const y0 = this.dataToViewY(0);
    const numStrips = xGrid.length - 1;

    // Ensure we have enough rectangles in the pool
    while (this.phaseColorStrips.length < numStrips) {
      const rect = new Rectangle(0, 0, 1, 1, {
        fill: "white",
        stroke: null,
      });
      this.phaseColorStrips.push(rect);
      this.phaseColorNode.addChild(rect);
    }

    // Update each strip
    for (let i = 0; i < numStrips; i++) {
      const x1 = this.dataToViewX(xGrid[i] * QuantumConstants.M_TO_NM);
      const x2 = this.dataToViewX(xGrid[i + 1] * QuantumConstants.M_TO_NM);
      const stripWidth = x2 - x1;

      const psi = wavefunction[i];
      const magnitude = Math.abs(psi);

      // Apply global time evolution phase
      const cosPhi = Math.cos(globalPhase);
      const sinPhi = Math.sin(globalPhase);

      // Apply time evolution: ψ(x,t) = ψ(x) * e^(-iEt/ℏ)
      const realPart = psi * cosPhi;
      const imagPart = -psi * sinPhi;

      // Calculate local phase: arg(ψ) = atan2(Im(ψ), Re(ψ))
      const localPhase = Math.atan2(imagPart, realPart);

      // Normalize phase to [0, 1] for hue (0 to 360 degrees)
      // phase ranges from -π to π, so we normalize to 0 to 1
      const normalizedPhase = (localPhase + Math.PI) / (2 * Math.PI);
      const hue = Math.round(normalizedPhase * 360); // 0 to 360 degrees (must be integer for CSS)

      // Create color using HSL: hue varies with phase, saturation and lightness are fixed
      const saturation = 80; // 80% saturation for vibrant colors
      const lightness = 60; // 60% lightness for good visibility
      const color = `hsl(${hue}, ${saturation}%, ${lightness}%)`;

      // Height of the strip is proportional to magnitude
      const yTop = this.dataToViewY(magnitude);
      const stripHeight = Math.abs(y0 - yTop);

      // Update the rectangle from the pool
      const strip = this.phaseColorStrips[i];
      strip.setRect(x1, Math.min(y0, yTop), stripWidth, stripHeight);
      strip.fill = color;
      strip.visible = stripHeight > 0.1; // Only show if visible
    }

    // Hide any extra strips we're not using
    for (let i = numStrips; i < this.phaseColorStrips.length; i++) {
      this.phaseColorStrips[i].visible = false;
    }
  }

  /**
   * Plots probability density from a pre-computed array.
   */
  private plotProbabilityDensityFromArray(
    xGrid: number[],
    probabilityDensity: number[],
  ): void {
    const shape = new Shape();

    // Start at zero on the left
    const x0 = this.dataToViewX(xGrid[0] * QuantumConstants.M_TO_NM);
    const y0 = this.dataToViewY(0);
    shape.moveTo(x0, y0);

    // Build points array
    const points: { x: number; y: number }[] = [];
    for (let i = 0; i < xGrid.length; i++) {
      const x = this.dataToViewX(xGrid[i] * QuantumConstants.M_TO_NM);
      const y = this.dataToViewY(probabilityDensity[i]);
      points.push({ x, y });
    }

    // Draw smooth curve using quadratic bezier curves
    if (points.length > 0) {
      shape.lineTo(points[0].x, points[0].y);

      for (let i = 0; i < points.length - 1; i++) {
        const p0 = points[i];
        const p1 = points[i + 1];

        // Control point is midpoint for simple smoothing
        const cpX = (p0.x + p1.x) / 2;
        const cpY = (p0.y + p1.y) / 2;

        shape.quadraticCurveTo(cpX, cpY, p1.x, p1.y);
      }
    }

    // Close at zero on the right
    const xEnd = this.dataToViewX(
      xGrid[xGrid.length - 1] * QuantumConstants.M_TO_NM,
    );
    shape.lineTo(xEnd, y0);
    shape.close();

    this.probabilityDensityPath.shape = shape;
  }

  /**
   * Plots superposition components (real, imaginary, magnitude) from arrays.
   */
  private plotSuperpositionComponents(
    xGrid: number[],
    realPart: number[],
    imagPart: number[],
    magnitude: number[],
  ): void {
    const realShape = new Shape();
    const imagShape = new Shape();
    const magShape = new Shape();

    // Build points arrays
    const realPoints: { x: number; y: number }[] = [];
    const imagPoints: { x: number; y: number }[] = [];
    const magPoints: { x: number; y: number }[] = [];

    for (let i = 0; i < xGrid.length; i++) {
      const x = this.dataToViewX(xGrid[i] * QuantumConstants.M_TO_NM);
      const yReal = this.dataToViewY(realPart[i]);
      const yImag = this.dataToViewY(imagPart[i]);
      const yMag = this.dataToViewY(magnitude[i]);

      realPoints.push({ x, y: yReal });
      imagPoints.push({ x, y: yImag });
      magPoints.push({ x, y: yMag });
    }

    // Draw smooth curves using quadratic bezier curves
    const drawSmoothCurve = (
      shape: Shape,
      points: { x: number; y: number }[],
    ) => {
      if (points.length > 0) {
        shape.moveTo(points[0].x, points[0].y);

        for (let i = 0; i < points.length - 1; i++) {
          const p0 = points[i];
          const p1 = points[i + 1];

          // Control point is midpoint for simple smoothing
          const cpX = (p0.x + p1.x) / 2;
          const cpY = (p0.y + p1.y) / 2;

          shape.quadraticCurveTo(cpX, cpY, p1.x, p1.y);
        }
      }
    };

    drawSmoothCurve(realShape, realPoints);
    drawSmoothCurve(imagShape, imagPoints);
    drawSmoothCurve(magShape, magPoints);

    this.realPartPath.shape = realShape;
    this.imaginaryPartPath.shape = imagShape;
    this.magnitudePath.shape = magShape;
  }

  /**
   * Plots phase-colored superposition from real and imaginary parts.
   */
  private plotPhaseColoredSuperposition(
    xGrid: number[],
    realPart: number[],
    imagPart: number[],
  ): void {
    const y0 = this.dataToViewY(0);
    const numStrips = xGrid.length - 1;

    // Ensure we have enough rectangles in the pool
    while (this.phaseColorStrips.length < numStrips) {
      const rect = new Rectangle(0, 0, 1, 1, {
        fill: "white",
        stroke: null,
      });
      this.phaseColorStrips.push(rect);
      this.phaseColorNode.addChild(rect);
    }

    // Update each strip
    for (let i = 0; i < numStrips; i++) {
      const x1 = this.dataToViewX(xGrid[i] * QuantumConstants.M_TO_NM);
      const x2 = this.dataToViewX(xGrid[i + 1] * QuantumConstants.M_TO_NM);
      const stripWidth = x2 - x1;

      const real = realPart[i];
      const imag = imagPart[i];
      const magnitude = Math.sqrt(real * real + imag * imag);

      // Calculate local phase: arg(ψ) = atan2(Im(ψ), Re(ψ))
      const localPhase = Math.atan2(imag, real);

      // Normalize phase to [0, 1] for hue (0 to 360 degrees)
      // phase ranges from -π to π, so we normalize to 0 to 1
      const normalizedPhase = (localPhase + Math.PI) / (2 * Math.PI);
      const hue = Math.round(normalizedPhase * 360); // 0 to 360 degrees (must be integer for CSS)

      // Create color using HSL: hue varies with phase, saturation and lightness are fixed
      const saturation = 80; // 80% saturation for vibrant colors
      const lightness = 60; // 60% lightness for good visibility
      const color = `hsl(${hue}, ${saturation}%, ${lightness}%)`;

      // Height of the strip is proportional to magnitude
      const yTop = this.dataToViewY(magnitude);
      const stripHeight = Math.abs(y0 - yTop);

      // Update the rectangle from the pool
      const strip = this.phaseColorStrips[i];
      strip.setRect(x1, Math.min(y0, yTop), stripWidth, stripHeight);
      strip.fill = color;
      strip.visible = stripHeight > 0.1; // Only show if visible
    }

    // Hide any extra strips we're not using
    for (let i = numStrips; i < this.phaseColorStrips.length; i++) {
      this.phaseColorStrips[i].visible = false;
    }
  }

  /**
   * Updates only the time evolution (for animation).
   */
  private updateTimeEvolution(): void {
    // Prevent reentry - if an update is in progress, mark that another update is pending
    if (this.isUpdating) {
      this.updatePending = true;
      return;
    }

    this.isUpdating = true;
    try {
      if (this.model.isPlayingProperty.value) {
        const boundStates = this.model.getBoundStates();
        if (boundStates) {
          // Check if we're displaying a superposition or a single eigenstate
          const superpositionType = this.model.superpositionTypeProperty.value;
          const isSuperposition =
            superpositionType !== SuperpositionType.SINGLE;

          if (isSuperposition && "superpositionConfigProperty" in this.model) {
            // Display superposition wavefunction with proper time evolution
            this.updateSuperpositionWavefunction();
          } else {
            // Display single eigenstate
            const selectedIndex =
              this.model.selectedEnergyLevelIndexProperty.value;
            if (
              selectedIndex >= 0 &&
              selectedIndex < boundStates.wavefunctions.length
            ) {
              this.updateWaveFunction(boundStates, selectedIndex);
            }
          }
        }
      }
    } finally {
      this.isUpdating = false;
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
