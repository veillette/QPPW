/**
 * WaveFunctionChartNode displays the wave function or probability density
 * for the selected energy state. This is the bottom chart in the One Well screen.
 */

import { Node, Line, Path, Text, VBox } from "scenerystack/scenery";
import { DerivedProperty } from "scenerystack/axon";
import { Orientation } from "scenerystack/phet-core";
import {
  AxisLine,
  TickMarkSet,
  TickLabelSet,
} from "scenerystack/bamboo";
import { Checkbox } from "scenerystack/sun";
import type { ScreenModel } from "../model/ScreenModels.js";
import {
  hasWellOffset,
  hasWellSeparation,
  hasSuperpositionConfig,
} from "../model/ModelTypeGuards.js";
import { BoundStateResult, PotentialType } from "../model/PotentialFunction.js";
import QuantumConstants from "../model/QuantumConstants.js";
import QPPWColors from "../../QPPWColors.js";
import { PhetFont } from "scenerystack/scenery-phet";
import { SuperpositionType } from "../model/SuperpositionType.js";
import QPPWPreferences from "../../QPPWPreferences.js";
import type { ScreenViewState } from "./ScreenViewStates.js";
import { AreaMeasurementTool } from "./chart-tools/AreaMeasurementTool.js";
import { CurvatureTool } from "./chart-tools/CurvatureTool.js";
import { DerivativeTool } from "./chart-tools/DerivativeTool.js";
import { ZerosVisualization } from "./chart-tools/ZerosVisualization.js";
import { PhaseColorVisualization } from "./chart-tools/PhaseColorVisualization.js";
import { ClassicalProbabilityOverlay } from "./chart-tools/ClassicalProbabilityOverlay.js";
import { ChartToolRegistry } from "./chart-tools/ChartToolRegistry.js";
import { BaseChartNode, ChartOptions } from "./BaseChartNode.js";
import {
  WaveFunctionDisplayStrategy,
  type WavefunctionData,
  type RenderContext,
  type RMSIndicatorContext,
} from "./wavefunction-strategies/WaveFunctionDisplayStrategy.js";
import { ProbabilityDensityDisplayStrategy } from "./wavefunction-strategies/ProbabilityDensityDisplayStrategy.js";
import { WaveFunctionComponentsDisplayStrategy } from "./wavefunction-strategies/WaveFunctionComponentsDisplayStrategy.js";
import { PhaseColorDisplayStrategy } from "./wavefunction-strategies/PhaseColorDisplayStrategy.js";

// Chart axis range constant (shared with EnergyChartNode)
const X_AXIS_RANGE_NM = 4; // X-axis extends from -X_AXIS_RANGE_NM to +X_AXIS_RANGE_NM

export class WaveFunctionChartNode extends BaseChartNode {
  // Visual elements specific to wavefunction chart
  private readonly realPartPath: Path;
  private readonly imaginaryPartPath: Path;
  private readonly magnitudePath: Path;
  private readonly probabilityDensityPath: Path;
  private readonly avgPositionIndicator: Line; // Vertical line indicator for average position
  private readonly rmsPositionIndicator: Path; // Double arrow indicator for RMS position
  private yAxisLabel!: Text;
  private readonly stateLabelNode: Text; // Label showing which wavefunction is displayed
  private readonly avgPositionLabel: Text;
  private readonly rmsPositionLabel: Text;

  // Tool components
  private readonly areaMeasurementTool: AreaMeasurementTool;
  private readonly curvatureTool: CurvatureTool;
  private readonly derivativeTool: DerivativeTool;
  private readonly zerosVisualization: ZerosVisualization;
  private readonly phaseColorVisualization: PhaseColorVisualization;
  private readonly classicalProbabilityOverlay: ClassicalProbabilityOverlay;

  // Tool registry for centralized management
  private readonly toolRegistry: ChartToolRegistry;

  // Display strategy for current display mode
  private displayStrategy: WaveFunctionDisplayStrategy;

  // Guard flag to prevent reentry during updates
  private isUpdating: boolean = false;
  // Flag to indicate if an update was requested while another update was in progress
  private updatePending: boolean = false;

  // Optional fixed display mode (overrides model's display mode)
  private readonly fixedDisplayMode?:
    | "probabilityDensity"
    | "waveFunction"
    | "phaseColor";

  // Public getters for tool show properties (for control panels)
  public get showAreaToolProperty() {
    return this.areaMeasurementTool.showProperty;
  }

  public get showCurvatureToolProperty() {
    return this.curvatureTool.showProperty;
  }

  public get showDerivativeToolProperty() {
    return this.derivativeTool.showProperty;
  }

  public constructor(
    model: ScreenModel,
    viewState: ScreenViewState,
    options?: {
      width?: number;
      height?: number;
      fixedDisplayMode?: "probabilityDensity" | "waveFunction" | "phaseColor";
      showToolCheckboxes?: boolean; // Whether to show curvature/derivative checkboxes (intro screen only)
    },
  ) {
    // Call BaseChartNode constructor with chart options
    const chartOptions: ChartOptions = {
      width: options?.width ?? 600,
      height: options?.height ?? 140,
      margins: { left: 60, right: 20, top: 10, bottom: 40 },
      xRange: { min: -X_AXIS_RANGE_NM, max: X_AXIS_RANGE_NM },
      yRange: { min: -1, max: 1 },
      showZeroLine: true,
    };

    super(model, viewState, chartOptions);

    // PDOM - make wavefunction chart accessible
    this.tagName = "div";
    this.labelTagName = "h3";
    this.labelContent = "Wavefunction Visualization";
    this.descriptionTagName = "p";

    // Set up accessible description after this.model is initialized
    this.descriptionContent = new DerivedProperty(
      [
        model.selectedEnergyLevelIndexProperty,
        model.potentialTypeProperty,
        model.superpositionTypeProperty,
      ],
      (
        selectedIndex: number,
        potentialType: PotentialType,
        superpositionType: SuperpositionType,
      ) => {
        return this.createWavefunctionDescription(
          selectedIndex,
          potentialType,
          superpositionType,
        );
      },
    );

    this.fixedDisplayMode = options?.fixedDisplayMode;

    // Create axes
    this.axesNode = this.createAxes();
    this.addChild(this.axesNode);

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

    // Create average position indicator (vertical line)
    this.avgPositionIndicator = new Line(0, 0, 0, 0, {
      stroke: QPPWColors.energyLevelSelectedProperty,
      lineWidth: 2,
      lineDash: [8, 4],
    });
    this.plotContentNode.addChild(this.avgPositionIndicator);

    // Create RMS position indicator (double arrow)
    this.rmsPositionIndicator = new Path(null, {
      stroke: QPPWColors.energyLevelSelectedProperty,
      lineWidth: 2,
      fill: QPPWColors.energyLevelSelectedProperty,
    });
    this.plotContentNode.addChild(this.rmsPositionIndicator);

    // Initialize tool components
    const toolOptions = {
      chartMargins: this.chartMargins,
      plotWidth: this.plotWidth,
      plotHeight: this.plotHeight,
      chartWidth: this.chartWidth,
      xMinProperty: this.xMinProperty,
      xMaxProperty: this.xMaxProperty,
      yMinProperty: this.yMinProperty,
      yMaxProperty: this.yMaxProperty,
      dataToViewX: this.dataToViewX.bind(this),
      dataToViewY: this.dataToViewY.bind(this),
      viewToDataX: this.viewToDataX.bind(this),
      parentNode: this,
    };

    // Create classical probability overlay
    this.classicalProbabilityOverlay = new ClassicalProbabilityOverlay(
      model,
      toolOptions,
    );
    this.plotContentNode.addChild(this.classicalProbabilityOverlay);

    // Create phase color visualization
    this.phaseColorVisualization = new PhaseColorVisualization({
      dataToViewX: this.dataToViewX.bind(this),
      dataToViewY: this.dataToViewY.bind(this),
    });
    this.plotContentNode.addChild(this.phaseColorVisualization);

    // Create zeros visualization
    this.zerosVisualization = new ZerosVisualization({
      dataToViewX: this.dataToViewX.bind(this),
      dataToViewY: this.dataToViewY.bind(this),
    });
    this.plotContentNode.addChild(this.zerosVisualization);

    // Create area measurement tool
    this.areaMeasurementTool = new AreaMeasurementTool(
      model,
      () => this.getEffectiveDisplayMode(),
      toolOptions,
    );
    this.plotContentNode.addChild(this.areaMeasurementTool);

    // Create curvature tool
    this.curvatureTool = new CurvatureTool(
      model,
      () => this.getEffectiveDisplayMode(),
      toolOptions,
    );
    this.plotContentNode.addChild(this.curvatureTool);

    // Create derivative tool
    this.derivativeTool = new DerivativeTool(
      model,
      () => this.getEffectiveDisplayMode(),
      toolOptions,
    );
    this.plotContentNode.addChild(this.derivativeTool);

    // Initialize tool registry and register all tools
    this.toolRegistry = new ChartToolRegistry();
    this.toolRegistry.registerTool("area", this.areaMeasurementTool);
    this.toolRegistry.registerTool("curvature", this.curvatureTool);
    this.toolRegistry.registerTool("derivative", this.derivativeTool);

    // Initialize display strategy based on current display mode
    this.displayStrategy = this.createDisplayStrategy(
      this.getEffectiveDisplayMode(),
    );

    // Create state label in upper right corner (outside clipped area)
    this.stateLabelNode = new Text("", {
      font: new PhetFont({ size: 16, style: "italic" }),
      fill: QPPWColors.labelFillProperty,
      right: this.chartWidth - this.chartMargins.right - 55,
      top: this.chartMargins.top + 5,
    });
    this.addChild(this.stateLabelNode);

    // Create labels for average and RMS position
    this.avgPositionLabel = new Text("", {
      font: new PhetFont(12),
      fill: QPPWColors.labelFillProperty,
      left: this.chartMargins.left + 10,
      top: this.chartMargins.top + 5,
    });
    this.addChild(this.avgPositionLabel);

    this.rmsPositionLabel = new Text("", {
      font: new PhetFont(12),
      fill: QPPWColors.labelFillProperty,
      left: this.chartMargins.left + 10,
      top: this.chartMargins.top + 25,
    });
    this.addChild(this.rmsPositionLabel);

    // Create checkboxes for derivative and curvature tools (only on intro screen)
    if (options?.showToolCheckboxes) {
      const curvatureCheckbox = new Checkbox(
        this.curvatureTool.showProperty,
        new Text("Show Curvature", {
          font: new PhetFont(12),
          fill: QPPWColors.textFillProperty,
        }),
        { boxWidth: 14 },
      );

      const derivativeCheckbox = new Checkbox(
        this.derivativeTool.showProperty,
        new Text("Show Derivative", {
          font: new PhetFont(12),
          fill: QPPWColors.textFillProperty,
        }),
        { boxWidth: 14 },
      );

      // Group checkboxes in a VBox
      const toolCheckboxes = new VBox({
        children: [curvatureCheckbox, derivativeCheckbox],
        spacing: 5,
        align: "left",
        right: this.chartWidth - this.chartMargins.right - 5,
        top: this.chartMargins.top + 5,
      });
      this.addChild(toolCheckboxes);

      // Only show checkboxes in wavefunction mode
      this.viewState.displayModeProperty.link((displayMode) => {
        const effectiveMode =
          this.fixedDisplayMode !== undefined
            ? this.fixedDisplayMode
            : displayMode;
        toolCheckboxes.visible = effectiveMode === "waveFunction";
      });
    }

    // Ensure axes are on top of the clipped plot content
    this.axesNode.moveToFront();

    // Link to model properties
    this.linkToModel();

    // Note: Initial update is now done asynchronously inside linkToModel()
    // to prevent blocking the page load
  }

  /**
   * Creates the appropriate display strategy based on the display mode.
   */
  private createDisplayStrategy(displayMode: string): WaveFunctionDisplayStrategy {
    switch (displayMode) {
      case "probabilityDensity":
        return new ProbabilityDensityDisplayStrategy();
      case "phaseColor":
        return new PhaseColorDisplayStrategy(this.phaseColorVisualization);
      case "waveFunction":
      default:
        return new WaveFunctionComponentsDisplayStrategy();
    }
  }

  /**
   * Creates an accessible description of the wavefunction chart based on current state.
   * This provides screen reader users with meaningful information about the visualization.
   */
  private createWavefunctionDescription(
    selectedIndex: number,
    _potentialType: PotentialType,
    superpositionType: SuperpositionType,
  ): string {
    const boundStates = this.model.getBoundStates();
    if (!boundStates || boundStates.energies.length === 0) {
      return "No wavefunction data available.";
    }

    const displayMode = this.getEffectiveDisplayMode();
    const isSuperposition = superpositionType !== SuperpositionType.SINGLE;

    let description = "";

    // Describe what is being shown
    if (isSuperposition) {
      description += `Superposition state wavefunction. `;
    } else {
      description += `Wavefunction for energy level ${selectedIndex + 1}. `;
    }

    // Display mode description
    if (displayMode === "probabilityDensity") {
      description += `Showing probability density |ψ(x)|². `;
      description += `This indicates where the particle is most likely to be found. `;
    } else if (displayMode === "phaseColor") {
      description += `Showing phase angle of complex wavefunction with color coding. `;
    } else {
      description += `Showing real and imaginary components of ψ(x). `;
    }

    // Get statistical properties if available
    const nmData = isSuperposition
      ? this.model.getTimeEvolvedSuperpositionInNmUnits(
          this.model.timeProperty.value * 1e-15,
        )
      : this.model.getWavefunctionInNmUnits(selectedIndex + 1);

    if (nmData) {
      const xGrid = boundStates.xGrid.map((x) => x * 1e9); // Convert to nm
      const { avg, rms } = calculateRMSStatistics(
        xGrid,
        nmData.probabilityDensity,
      );

      description += `\n\nPosition statistics: `;
      description += `Average position: ${avg.toFixed(2)} nanometers. `;
      description += `Position uncertainty (RMS, Δx): ${rms.toFixed(2)} nm. `;
    }

    // Node count for single eigenstates
    if (!isSuperposition && selectedIndex >= 0) {
      const nodes = selectedIndex; // Quantum number n-1
      description += `\n\nWavefunction has ${nodes} node${nodes !== 1 ? "s" : ""} `;
      description += `(zero crossing${nodes !== 1 ? "s" : ""}). `;
    }

    return description;
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
            lineDash: [5, 5],
          },
        );
        axesNode.addChild(gridLine);
      }
    }

    // X-axis line at the bottom
    const xAxis = new AxisLine(this.chartTransform, Orientation.HORIZONTAL);
    xAxis.x = this.chartMargins.left;
    xAxis.y = this.chartMargins.top + this.plotHeight;
    axesNode.addChild(xAxis);

    // X-axis tick marks
    const xTickMarkSet = new TickMarkSet(
      this.chartTransform,
      Orientation.HORIZONTAL,
      2, // spacing
      {
        edge: "max",
        extent: 8,
        stroke: QPPWColors.labelFillProperty,
        lineWidth: 1,
      },
    );
    xTickMarkSet.x = this.chartMargins.left;
    xTickMarkSet.y = this.chartMargins.top + this.plotHeight;
    axesNode.addChild(xTickMarkSet);

    // X-axis tick labels
    const xTickLabelSet = new TickLabelSet(
      this.chartTransform,
      Orientation.HORIZONTAL,
      2, // spacing
      {
        edge: "max",
        createLabel: (value: number) =>
          new Text(value.toString(), {
            font: new PhetFont(14),
            fill: QPPWColors.labelFillProperty,
          }),
      },
    );
    xTickLabelSet.x = this.chartMargins.left;
    xTickLabelSet.y = this.chartMargins.top + this.plotHeight;
    axesNode.addChild(xTickLabelSet);

    // Y-axis line on the left
    const yAxis = new AxisLine(this.chartTransform, Orientation.VERTICAL);
    yAxis.x = this.chartMargins.left;
    yAxis.y = this.chartMargins.top;
    axesNode.addChild(yAxis);

    // Y-axis tick marks using bamboo TickMarkSet
    // Use fixed spacing of 0.5 (works well for typical nm^-1 and nm^-1/2 values)
    const yTickMarkSet = new TickMarkSet(
      this.chartTransform,
      Orientation.VERTICAL,
      0.5, // spacing in nm^-1 or nm^-1/2 units
      {
        edge: "min",
        extent: 8,
        stroke: QPPWColors.labelFillProperty,
        lineWidth: 1,
      },
    );
    yTickMarkSet.x = this.chartMargins.left;
    yTickMarkSet.y = this.chartMargins.top;
    axesNode.addChild(yTickMarkSet);

    // Y-axis tick labels using bamboo TickLabelSet
    const yTickLabelSet = new TickLabelSet(
      this.chartTransform,
      Orientation.VERTICAL,
      0.5, // spacing in nm^-1 or nm^-1/2 units
      {
        edge: "min",
        createLabel: (value: number) =>
          new Text(this.formatYTickLabel(value), {
            font: new PhetFont(12),
            fill: QPPWColors.labelFillProperty,
          }),
      },
    );
    yTickLabelSet.x = this.chartMargins.left;
    yTickLabelSet.y = this.chartMargins.top;
    axesNode.addChild(yTickLabelSet);

    // Y-axis label - centered on the visual center of the plot area
    // Account for labels at top of chart (state label, avg/rms labels take ~40px)
    const labelAreaHeight = 40;
    const visualTop = this.chartMargins.top + labelAreaHeight;
    const visualBottom = this.chartMargins.top + this.plotHeight;
    this.yAxisLabel = new Text("", {
      font: new PhetFont(14),
      fill: QPPWColors.labelFillProperty,
      centerX: 15,
      centerY: (visualTop + visualBottom) / 2,
      rotation: -Math.PI / 2,
    });
    axesNode.addChild(this.yAxisLabel);

    // X-axis label
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
   * Links chart updates to model property changes.
   */
  private linkToModel(): void {
    // Update when any parameter changes
    this.model.potentialTypeProperty.lazyLink(() => this.update());
    this.model.wellWidthProperty.lazyLink(() => this.update());
    this.model.wellDepthProperty.lazyLink(() => this.update());
    if (hasWellOffset(this.model)) {
      this.model.wellOffsetProperty.lazyLink(() => this.update());
    }
    this.model.particleMassProperty.lazyLink(() => this.update());

    // Link to wellSeparationProperty if available (TwoWellsModel and ManyWellsModel)
    if (hasWellSeparation(this.model)) {
      this.model.wellSeparationProperty.lazyLink(() => this.update());
    }

    // Update when grid points preference changes (affects wavefunction resolution)
    QPPWPreferences.gridPointsProperty.lazyLink(() => {
      // Force model to invalidate its cache by notifying it of the change
      // The BaseModel listener will handle cache invalidation
      this.update();
    });

    this.model.selectedEnergyLevelIndexProperty.lazyLink(() => {
      this.updateStateLabel();
      this.update();
    });
    this.viewState.displayModeProperty.lazyLink(() => {
      // Recreate display strategy when mode changes
      const newMode = this.getEffectiveDisplayMode();
      const oldStrategy = this.displayStrategy;
      this.displayStrategy = this.createDisplayStrategy(newMode);

      // Hide old strategy's visualization if needed (for phase color)
      if (oldStrategy instanceof PhaseColorDisplayStrategy) {
        oldStrategy.hide();
      }

      this.updateYAxisLabel();
      this.updateStateLabel();
      this.update();
    });
    this.model.timeProperty.lazyLink(() => this.updateTimeEvolution());

    // Update when superposition configuration changes
    this.model.superpositionTypeProperty.lazyLink(() => {
      this.updateStateLabel();
      this.update();
    });
    this.model.superpositionConfigProperty.lazyLink(() => {
      this.update();
    });

    // Update visibility of wave function components
    this.viewState.showRealPartProperty.lazyLink((show: boolean) => {
      this.realPartPath.visible =
        show && this.getEffectiveDisplayMode() === "waveFunction";
    });
    this.viewState.showImaginaryPartProperty.lazyLink((show: boolean) => {
      this.imaginaryPartPath.visible =
        show && this.getEffectiveDisplayMode() === "waveFunction";
    });
    this.viewState.showMagnitudeProperty.lazyLink((show: boolean) => {
      this.magnitudePath.visible =
        show && this.getEffectiveDisplayMode() === "waveFunction";
    });

    // Update visibility of classical probability
    this.viewState.showClassicalProbabilityProperty.lazyLink(() => {
      this.update();
    });

    // Update visibility of zeros
    this.viewState.showZerosProperty.lazyLink(() => {
      this.update();
    });

    // Update visibility of RMS indicators if the property exists (IntroViewState only)
    if ("showRMSIndicatorProperty" in this.viewState) {
      this.viewState.showRMSIndicatorProperty.lazyLink(() => {
        this.update();
      });
    }

    // Link tool updates to model changes using the tool registry
    const updateTools = () => {
      const displayMode = this.getEffectiveDisplayMode();
      const boundStates = this.model.getBoundStates();
      const selectedIndex = this.model.selectedEnergyLevelIndexProperty.value;

      this.toolRegistry.updateAllTools({
        displayMode,
        boundStates,
        selectedIndex,
      });
    };

    // Update tools when display mode changes
    this.viewState.displayModeProperty.lazyLink(() => {
      updateTools();
    });

    // Update tools when time changes (for superposition states)
    this.model.timeProperty.lazyLink(() => {
      if (this.model.isPlayingProperty.value) {
        updateTools();
      }
    });

    // Update tools when potential or well parameters change
    this.model.potentialTypeProperty.lazyLink(updateTools);
    this.model.wellWidthProperty.lazyLink(updateTools);
    this.model.wellDepthProperty.lazyLink(updateTools);

    if ("barrierHeightProperty" in this.model) {
      this.model.barrierHeightProperty.lazyLink(updateTools);
    }

    if ("potentialOffsetProperty" in this.model) {
      this.model.potentialOffsetProperty.lazyLink(updateTools);
    }

    // Update tools when selected energy level changes
    this.model.selectedEnergyLevelIndexProperty.lazyLink(updateTools);

    if ("wellSeparationProperty" in this.model) {
      this.model.wellSeparationProperty.lazyLink(updateTools);
    }

    if ("numberOfWellsProperty" in this.model) {
      this.model.numberOfWellsProperty.lazyLink(updateTools);
    }

    if ("electricFieldProperty" in this.model) {
      this.model.electricFieldProperty.lazyLink(updateTools);
    }

    // Initialize labels (important for fixed display mode charts)
    this.updateYAxisLabel();
    this.updateStateLabel();

    // Perform initial update and tool updates asynchronously (after construction completes)
    // This prevents blocking the page load with expensive calculations
    setTimeout(() => {
      this.update();
      updateTools();
    }, 0);
  }

  /**
   * Gets the effective display mode, using the fixed display mode if set,
   * otherwise falling back to the model's display mode.
   */
  private getEffectiveDisplayMode(): string {
    return this.fixedDisplayMode || this.viewState.displayModeProperty.value;
  }

  /**
   * Checks if RMS indicators should be shown based on viewState property.
   * Returns true if the property doesn't exist (for backwards compatibility with other screens).
   */
  private shouldShowRMSIndicators(): boolean {
    if ("showRMSIndicatorProperty" in this.viewState) {
      return this.viewState.showRMSIndicatorProperty.value;
    }
    return true; // Show by default if property doesn't exist
  }

  /**
   * Updates the Y-axis label based on display mode.
   */
  private updateYAxisLabel(): void {
    this.yAxisLabel.string = this.displayStrategy.getYAxisLabel();
  }

  /**
   * Updates the state label showing which wavefunction is displayed.
   */
  private updateStateLabel(): void {
    const superpositionType = this.model.superpositionTypeProperty.value;
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

      // Format label using display strategy
      this.stateLabelNode.string = this.displayStrategy.getStateLabel(label);
    } else {
      // Display single eigenstate label
      const selectedIndex = this.model.selectedEnergyLevelIndexProperty.value;
      const boundStates = this.model.getBoundStates();

      if (!boundStates || selectedIndex < 0) {
        this.stateLabelNode.string = "";
        return;
      }

      if (selectedIndex >= boundStates.wavefunctions.length) {
        this.stateLabelNode.string = "";
        return;
      }

      // Format: ψ₁, ψ₂, etc. (display strategy handles mode-specific formatting)
      const stateNumber = selectedIndex + 1;
      const stateLabel = `ψ${this.toSubscript(stateNumber)}`;

      // Format label using display strategy
      this.stateLabelNode.string = this.displayStrategy.getStateLabel(stateLabel);
    }
  }

  /**
   * Converts a number to subscript Unicode characters.
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
   * Formats a y-axis tick label value for display.
   * Uses appropriate precision based on the magnitude of the value.
   */
  private formatYTickLabel(value: number): string {
    // For values close to zero, show as 0.00
    if (Math.abs(value) < 1e-10) {
      return "0.00";
    }

    // Use consistent 2 decimal places for all tick values
    return value.toFixed(2);
  }

  public update(): void {
    // Prevent reentry - if an update is in progress, mark that another update is pending
    if (this.isUpdating) {
      this.updatePending = true;
      return;
    }

    this.isUpdating = true;
    try {
      // Keep updating until no more updates are pending
      let loopCount = 0;
      do {
        this.updatePending = false;
        loopCount++;
        if (loopCount > 10) {
          console.warn(
            "[WaveFunctionChartNode] Infinite loop detected! Breaking after",
            loopCount,
            "iterations",
          );
          break;
        }

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
          // Hide classical probability visualization for superpositions
          this.classicalProbabilityOverlay.hide();
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

          // Update classical probability visualization
          // Show forbidden regions on both wavefunction and probability density charts
          // But only show the classical probability curve on probability density chart
          const displayMode = this.getEffectiveDisplayMode();
          const showOverlay =
            this.viewState.showClassicalProbabilityProperty.value;
          const showCurve = displayMode === "probabilityDensity";
          this.classicalProbabilityOverlay.update(
            boundStates,
            selectedIndex,
            showOverlay,
            showCurve,
          );
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
   * Uses nm units (nm^-1/2 for wavefunction, nm^-1 for probability density).
   */
  private updateViewRange(
    boundStates: BoundStateResult,
    selectedIndex: number,
  ): void {
    // Calculate Y range based on wave function values
    if (
      selectedIndex < 0 ||
      selectedIndex >= boundStates.wavefunctions.length
    ) {
      return;
    }

    // Get wavefunction and probability density in nm units
    const nmData = this.model.getWavefunctionInNmUnits(selectedIndex + 1);
    if (!nmData) {
      return;
    }

    // Create wavefunction data structure
    const wavefunctionData: WavefunctionData = {
      xGrid: boundStates.xGrid,
      realPart: nmData.wavefunction.map((psi) => psi), // For single state, wavefunction is real
      imagPart: nmData.wavefunction.map(() => 0), // Single eigenstate has no imaginary part at t=0
      magnitude: nmData.wavefunction.map(Math.abs),
      probabilityDensity: nmData.probabilityDensity,
      maxMagnitude: Math.max(...nmData.wavefunction.map(Math.abs)),
    };

    // Use display strategy to calculate Y range
    const range = this.displayStrategy.calculateYRange(wavefunctionData);

    // Update properties
    this.yMinProperty.value = range.min;
    this.yMaxProperty.value = range.max;

    // Update chart transform
    this.chartTransform.setModelYRange(range);
  }

  /**
   * Updates the view range for superposition wavefunctions.
   * Uses nm units (nm^-1/2 for wavefunction, nm^-1 for probability density).
   */
  private updateViewRangeForSuperpositionFromModel(): void {
    const boundStates = this.model.getBoundStates();
    if (!boundStates || !hasSuperpositionConfig(this.model)) {
      return;
    }

    const time = this.model.timeProperty.value * 1e-15; // Convert fs to seconds

    // Get time-evolved superposition in nm units
    const nmData = this.model.getTimeEvolvedSuperpositionInNmUnits(time);
    if (!nmData) {
      return;
    }

    // Create wavefunction data structure
    const wavefunctionData: WavefunctionData = {
      xGrid: boundStates.xGrid,
      realPart: nmData.realPart,
      imagPart: nmData.imagPart,
      magnitude: nmData.realPart.map((re, i) => {
        const im = nmData.imagPart[i];
        return Math.sqrt(re * re + im * im);
      }),
      probabilityDensity: nmData.probabilityDensity,
      maxMagnitude: nmData.maxMagnitude,
    };

    // Use display strategy to calculate Y range
    const range = this.displayStrategy.calculateYRange(wavefunctionData);

    // Update properties
    this.yMinProperty.value = range.min;
    this.yMaxProperty.value = range.max;

    // Update chart transform
    this.chartTransform.setModelYRange(range);
  }

  /**
   * Updates the superposition wavefunction visualization.
   * Uses nm units (nm^-1/2 for wavefunction, nm^-1 for probability density).
   */
  private updateSuperpositionWavefunction(): void {
    const boundStates = this.model.getBoundStates();
    if (!boundStates || !hasSuperpositionConfig(this.model)) {
      return;
    }

    const time = this.model.timeProperty.value * 1e-15; // Convert fs to seconds

    // Get time-evolved superposition in nm units
    const nmData = this.model.getTimeEvolvedSuperpositionInNmUnits(time);
    if (!nmData) {
      return;
    }

    const xGrid = boundStates.xGrid;

    // Create wavefunction data structure
    const wavefunctionData: WavefunctionData = {
      xGrid,
      realPart: nmData.realPart,
      imagPart: nmData.imagPart,
      magnitude: nmData.realPart.map((re, i) => {
        const im = nmData.imagPart[i];
        return Math.sqrt(re * re + im * im);
      }),
      probabilityDensity: nmData.probabilityDensity,
      maxMagnitude: nmData.maxMagnitude,
    };

    // Create render context
    const renderContext: RenderContext = {
      realPartPath: this.realPartPath,
      imaginaryPartPath: this.imaginaryPartPath,
      magnitudePath: this.magnitudePath,
      probabilityDensityPath: this.probabilityDensityPath,
      dataToViewX: this.dataToViewX.bind(this),
      dataToViewY: this.dataToViewY.bind(this),
      chartMargins: this.chartMargins,
      yMinProperty: this.yMinProperty,
      yMaxProperty: this.yMaxProperty,
      showRealPart: this.viewState.showRealPartProperty.value,
      showImaginaryPart: this.viewState.showImaginaryPartProperty.value,
      showMagnitude: this.viewState.showMagnitudeProperty.value,
    };

    // Delegate rendering to display strategy
    this.displayStrategy.render(wavefunctionData, renderContext);

    // Create RMS indicator context
    const rmsContext: RMSIndicatorContext = {
      avgPositionIndicator: this.avgPositionIndicator,
      rmsPositionIndicator: this.rmsPositionIndicator,
      avgPositionLabel: this.avgPositionLabel,
      rmsPositionLabel: this.rmsPositionLabel,
      shouldShow: this.shouldShowRMSIndicators(),
      dataToViewX: this.dataToViewX.bind(this),
      dataToViewY: this.dataToViewY.bind(this),
      yMaxValue: this.yMaxProperty.value,
    };

    // Delegate RMS indicator updates to display strategy
    this.displayStrategy.updateRMSIndicators(wavefunctionData, rmsContext);

    // Update zeros visualization if enabled and supported by display mode
    // Get SI units for zeros visualization (which still uses SI)
    const siData = this.model.getTimeEvolvedSuperposition(time);
    const realPartSI = siData ? siData.realPart : nmData.realPart;

    if (
      this.viewState.showZerosProperty.value &&
      this.displayStrategy.shouldShowZeros()
    ) {
      this.zerosVisualization.showProperty.value = true;
      this.zerosVisualization.update(xGrid, realPartSI);
    } else {
      this.zerosVisualization.showProperty.value = false;
    }
  }

  /**
   * Updates the zero line position.
   */
  private updateZeroLine(): void {
    const y = this.dataToViewY(0);
    this.zeroLine.setLine(
      this.chartMargins.left,
      y,
      this.chartMargins.left + this.plotWidth,
      y,
    );
  }

  /**
   * Updates the wave function visualization for a single eigenstate.
   */
  private updateWaveFunction(
    boundStates: BoundStateResult,
    selectedIndex: number,
  ): void {
    if (
      selectedIndex < 0 ||
      selectedIndex >= boundStates.wavefunctions.length
    ) {
      return;
    }

    // Get wavefunction and probability density in nm units
    const nmData = this.model.getWavefunctionInNmUnits(selectedIndex + 1);
    if (!nmData) {
      return;
    }

    const xGrid = boundStates.xGrid;
    const wavefunctionNm = nmData.wavefunction;
    const wavefunctionSI = boundStates.wavefunctions[selectedIndex];

    // Calculate time evolution phase for single eigenstates
    const energy = boundStates.energies[selectedIndex];
    const time = this.model.timeProperty.value * 1e-15; // Convert fs to seconds
    const globalPhase = -(energy * time) / QuantumConstants.HBAR;

    // Apply time evolution to get real and imaginary parts (for wavefunction mode)
    const realPart = wavefunctionNm.map((psi) => psi * Math.cos(globalPhase));
    const imagPart = wavefunctionNm.map((psi) => -psi * Math.sin(globalPhase));

    // Create wavefunction data structure
    const wavefunctionData: WavefunctionData = {
      xGrid,
      realPart,
      imagPart,
      magnitude: wavefunctionNm.map(Math.abs),
      probabilityDensity: nmData.probabilityDensity,
      maxMagnitude: Math.max(...wavefunctionNm.map(Math.abs)),
    };

    // Create render context
    const renderContext: RenderContext = {
      realPartPath: this.realPartPath,
      imaginaryPartPath: this.imaginaryPartPath,
      magnitudePath: this.magnitudePath,
      probabilityDensityPath: this.probabilityDensityPath,
      dataToViewX: this.dataToViewX.bind(this),
      dataToViewY: this.dataToViewY.bind(this),
      chartMargins: this.chartMargins,
      yMinProperty: this.yMinProperty,
      yMaxProperty: this.yMaxProperty,
      showRealPart: this.viewState.showRealPartProperty.value,
      showImaginaryPart: this.viewState.showImaginaryPartProperty.value,
      showMagnitude: this.viewState.showMagnitudeProperty.value,
    };

    // Delegate rendering to display strategy
    // Phase color mode needs special handling for single eigenstate with global phase
    if (this.displayStrategy instanceof PhaseColorDisplayStrategy) {
      this.displayStrategy.renderSingleState(
        wavefunctionData,
        renderContext,
        globalPhase,
      );
    } else {
      this.displayStrategy.render(wavefunctionData, renderContext);
    }

    // Create RMS indicator context
    const rmsContext: RMSIndicatorContext = {
      avgPositionIndicator: this.avgPositionIndicator,
      rmsPositionIndicator: this.rmsPositionIndicator,
      avgPositionLabel: this.avgPositionLabel,
      rmsPositionLabel: this.rmsPositionLabel,
      shouldShow: this.shouldShowRMSIndicators(),
      dataToViewX: this.dataToViewX.bind(this),
      dataToViewY: this.dataToViewY.bind(this),
      yMaxValue: this.yMaxProperty.value,
    };

    // Delegate RMS indicator updates to display strategy
    this.displayStrategy.updateRMSIndicators(wavefunctionData, rmsContext);

    // Update zeros visualization if enabled and supported by display mode
    if (
      this.viewState.showZerosProperty.value &&
      this.displayStrategy.shouldShowZeros()
    ) {
      this.zerosVisualization.showProperty.value = true;
      this.zerosVisualization.update(xGrid, wavefunctionSI);
    } else {
      this.zerosVisualization.showProperty.value = false;
    }
  }

  /**
   * Updates the wave function visualization during time evolution.
   * This is optimized to only update what changes with time.
   */
  private updateTimeEvolution(): void {
    // Only update if we're in a mode that shows time-dependent changes
    const displayMode = this.getEffectiveDisplayMode();
    const superpositionType = this.model.superpositionTypeProperty.value;
    const isSuperposition = superpositionType !== SuperpositionType.SINGLE;

    if (isSuperposition) {
      // Superposition wavefunctions evolve in time
      this.updateSuperpositionWavefunction();
    } else if (displayMode === "waveFunction" || displayMode === "phaseColor") {
      // Single eigenstates show time evolution in waveFunction and phaseColor modes
      const boundStates = this.model.getBoundStates();
      if (!boundStates) {
        return;
      }

      const selectedIndex = this.model.selectedEnergyLevelIndexProperty.value;
      if (
        selectedIndex < 0 ||
        selectedIndex >= boundStates.wavefunctions.length
      ) {
        return;
      }

      this.updateWaveFunction(boundStates, selectedIndex);
    }
    // For probability density mode with single eigenstates, nothing changes with time
  }
}
