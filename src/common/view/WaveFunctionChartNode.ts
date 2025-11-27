/**
 * WaveFunctionChartNode displays the wave function or probability density
 * for the selected energy state. This is the bottom chart in the One Well screen.
 */

import { Node, Line, Path, Text, Rectangle, Circle, DragListener } from "scenerystack/scenery";
import { Shape } from "scenerystack/kite";
import { NumberProperty, BooleanProperty } from "scenerystack/axon";
import { Range } from "scenerystack/dot";
import { Orientation } from "scenerystack/phet-core";
import {
  ChartTransform,
  ChartRectangle,
  AxisLine,
  TickMarkSet,
  TickLabelSet,
} from "scenerystack/bamboo";
import type { ScreenModel } from "../model/ScreenModels.js";
import {
  hasWellOffset,
  hasWellSeparation,
  hasSuperpositionConfig,
  hasClassicallyForbiddenProbability,
  hasClassicalTurningPoints,
} from "../model/ModelTypeGuards.js";
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
  private readonly model: ScreenModel;
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
  private readonly classicalProbabilityPath: Path;
  private readonly phaseColorNode: Node; // Container for phase-colored visualization
  private phaseColorStrips: Rectangle[]; // Pool of rectangles for phase visualization
  private readonly zeroLine: Line;
  private readonly axesNode: Node;
  private yAxisLabel!: Text;
  private readonly stateLabelNode: Text; // Label showing which wavefunction is displayed

  // Classical probability visualization
  private readonly leftTurningPointLine: Line;
  private readonly rightTurningPointLine: Line;
  private readonly leftForbiddenRegion: Rectangle;
  private readonly rightForbiddenRegion: Rectangle;
  private readonly forbiddenProbabilityLabel: Text; // Label showing the forbidden probability percentage

  // Zeros visualization
  private readonly zerosNode: Node; // Container for zero markers

  // Area measurement tool
  public readonly showAreaToolProperty: BooleanProperty;
  private readonly areaToolContainer: Node; // Container for all area tool elements
  private readonly leftMarker: Line;
  private readonly rightMarker: Line;
  private readonly leftMarkerHandle: Circle;
  private readonly rightMarkerHandle: Circle;
  private readonly areaBackgroundRegion: Rectangle; // Faint rectangular background
  private readonly areaRegion: Path; // Area under the curve
  private readonly areaLabel: Text;
  private leftMarkerXProperty: NumberProperty; // X position in nm
  private rightMarkerXProperty: NumberProperty; // X position in nm

  // Curvature visualization tool
  public readonly showCurvatureToolProperty: BooleanProperty;
  private readonly curvatureToolContainer: Node; // Container for all curvature tool elements
  private readonly curvatureMarker: Line;
  private readonly curvatureMarkerHandle: Circle;
  private readonly curvatureParabola: Path; // Parabola showing the curvature
  private readonly curvatureLabel: Text;
  private curvatureXProperty: NumberProperty; // X position in nm

  // Guard flag to prevent reentry during updates
  private isUpdating: boolean = false;
  // Flag to indicate if an update was requested while another update was in progress
  private updatePending: boolean = false;

  // Optional fixed display mode (overrides model's display mode)
  private readonly fixedDisplayMode?:
    | "probabilityDensity"
    | "waveFunction"
    | "phaseColor";

  public constructor(
    model: ScreenModel,
    options?: {
      width?: number;
      height?: number;
      fixedDisplayMode?: "probabilityDensity" | "waveFunction" | "phaseColor";
    },
  ) {
    super();

    this.model = model;
    this.chartWidth = options?.width ?? 600;
    this.chartHeight = options?.height ?? 140;
    this.fixedDisplayMode = options?.fixedDisplayMode;

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

    // Create classical turning point lines
    this.leftTurningPointLine = new Line(0, 0, 0, 0, {
      stroke: QPPWColors.energyLevelSelectedProperty,
      lineWidth: 2,
      lineDash: [8, 4],
      visible: false,
    });
    this.plotContentNode.addChild(this.leftTurningPointLine);

    this.rightTurningPointLine = new Line(0, 0, 0, 0, {
      stroke: QPPWColors.energyLevelSelectedProperty,
      lineWidth: 2,
      lineDash: [8, 4],
      visible: false,
    });
    this.plotContentNode.addChild(this.rightTurningPointLine);

    // Create classically forbidden regions (shaded areas)
    this.leftForbiddenRegion = new Rectangle(0, 0, 1, 1, {
      fill: "rgba(255, 200, 200, 0.3)",
      stroke: null,
      visible: false,
    });
    this.plotContentNode.addChild(this.leftForbiddenRegion);

    this.rightForbiddenRegion = new Rectangle(0, 0, 1, 1, {
      fill: "rgba(255, 200, 200, 0.3)",
      stroke: null,
      visible: false,
    });
    this.plotContentNode.addChild(this.rightForbiddenRegion);

    // Create forbidden probability label (appears on hover)
    this.forbiddenProbabilityLabel = new Text("", {
      font: new PhetFont(14),
      fill: QPPWColors.labelFillProperty,
      visible: false,
      centerX: this.chartWidth / 2,
      top: this.chartMargins.top + 30,
    });
    this.addChild(this.forbiddenProbabilityLabel);

    // Add hover listeners to forbidden regions to show probability percentage
    const showForbiddenProbability = () => {
      const selectedIndex = this.model.selectedEnergyLevelIndexProperty.value;
      if (
        hasClassicallyForbiddenProbability(this.model) &&
        selectedIndex >= 0
      ) {
        const percentage =
          this.model.getClassicallyForbiddenProbability(selectedIndex);
        this.forbiddenProbabilityLabel.string = `Classically Forbidden: ${percentage.toFixed(1)}%`;
        this.forbiddenProbabilityLabel.visible = true;
      }
    };

    const hideForbiddenProbability = () => {
      this.forbiddenProbabilityLabel.visible = false;
    };

    this.leftForbiddenRegion.addInputListener({
      enter: showForbiddenProbability,
      exit: hideForbiddenProbability,
    });

    this.rightForbiddenRegion.addInputListener({
      enter: showForbiddenProbability,
      exit: hideForbiddenProbability,
    });

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

    this.classicalProbabilityPath = new Path(null, {
      stroke: QPPWColors.classicalProbabilityProperty,
      lineWidth: 2,
      lineDash: [5, 3],
      visible: false,
    });
    this.plotContentNode.addChild(this.classicalProbabilityPath);

    // Create phase color visualization node
    this.phaseColorNode = new Node({
      visible: false,
    });
    this.plotContentNode.addChild(this.phaseColorNode);

    // Initialize rectangle pool for phase color visualization (will be populated on first use)
    this.phaseColorStrips = [];

    // Create zeros visualization node
    this.zerosNode = new Node({
      visible: false,
    });
    this.plotContentNode.addChild(this.zerosNode);

    // Initialize area measurement tool
    this.showAreaToolProperty = new BooleanProperty(false);

    // Initialize marker positions (default to -1nm and +1nm)
    this.leftMarkerXProperty = new NumberProperty(-1);
    this.rightMarkerXProperty = new NumberProperty(1);

    // Create area tool container
    this.areaToolContainer = new Node({
      visible: false,
    });
    this.plotContentNode.addChild(this.areaToolContainer);

    // Create faint background rectangle showing the measurement region
    this.areaBackgroundRegion = new Rectangle(0, 0, 1, 1, {
      fill: "rgba(100, 150, 255, 0.1)", // Very faint blue
      stroke: null,
    });
    this.areaToolContainer.addChild(this.areaBackgroundRegion);

    // Create shaded region between markers (follows the curve)
    this.areaRegion = new Path(null, {
      fill: "rgba(100, 150, 255, 0.3)", // Semi-transparent blue
      stroke: null,
    });
    this.areaToolContainer.addChild(this.areaRegion);

    // Create left marker line
    this.leftMarker = new Line(0, 0, 0, 0, {
      stroke: QPPWColors.energyLevelSelectedProperty,
      lineWidth: 2,
      lineDash: [6, 4],
    });
    this.areaToolContainer.addChild(this.leftMarker);

    // Create right marker line
    this.rightMarker = new Line(0, 0, 0, 0, {
      stroke: QPPWColors.energyLevelSelectedProperty,
      lineWidth: 2,
      lineDash: [6, 4],
    });
    this.areaToolContainer.addChild(this.rightMarker);

    // Create left marker handle (draggable circle at top)
    this.leftMarkerHandle = new Circle(8, {
      fill: QPPWColors.energyLevelSelectedProperty,
      stroke: QPPWColors.backgroundColorProperty,
      lineWidth: 2,
      cursor: "ew-resize",
    });
    this.areaToolContainer.addChild(this.leftMarkerHandle);

    // Create right marker handle (draggable circle at top)
    this.rightMarkerHandle = new Circle(8, {
      fill: QPPWColors.energyLevelSelectedProperty,
      stroke: QPPWColors.backgroundColorProperty,
      lineWidth: 2,
      cursor: "ew-resize",
    });
    this.areaToolContainer.addChild(this.rightMarkerHandle);

    // Create area percentage label
    this.areaLabel = new Text("", {
      font: new PhetFont({ size: 16, weight: "bold" }),
      fill: QPPWColors.labelFillProperty,
      visible: false,
    });
    this.addChild(this.areaLabel);

    // Setup drag listeners for markers
    this.setupAreaToolDragListeners();

    // Link area tool visibility to property
    this.showAreaToolProperty.link((show: boolean) => {
      this.areaToolContainer.visible = show;
      this.areaLabel.visible = show;
      if (show) {
        this.updateAreaTool();
      }
    });

    // Initialize curvature visualization tool
    this.showCurvatureToolProperty = new BooleanProperty(false);

    // Initialize marker position (default to 0nm, center)
    this.curvatureXProperty = new NumberProperty(0);

    // Create curvature tool container
    this.curvatureToolContainer = new Node({
      visible: false,
    });
    this.plotContentNode.addChild(this.curvatureToolContainer);

    // Create marker line
    this.curvatureMarker = new Line(0, 0, 0, 0, {
      stroke: "rgba(255, 100, 100, 0.8)", // Semi-transparent red
      lineWidth: 2,
      lineDash: [4, 3],
    });
    this.curvatureToolContainer.addChild(this.curvatureMarker);

    // Create marker handle (draggable circle)
    this.curvatureMarkerHandle = new Circle(8, {
      fill: "rgba(255, 100, 100, 0.9)",
      stroke: QPPWColors.backgroundColorProperty,
      lineWidth: 2,
      cursor: "ew-resize",
    });
    this.curvatureToolContainer.addChild(this.curvatureMarkerHandle);

    // Create parabola path
    this.curvatureParabola = new Path(null, {
      stroke: "rgba(255, 100, 100, 0.9)",
      lineWidth: 3,
      fill: null,
    });
    this.curvatureToolContainer.addChild(this.curvatureParabola);

    // Create curvature label
    this.curvatureLabel = new Text("", {
      font: new PhetFont({ size: 14, weight: "bold" }),
      fill: "rgba(255, 100, 100, 1)",
      visible: false,
    });
    this.addChild(this.curvatureLabel);

    // Setup drag listener for marker
    this.setupCurvatureToolDragListener();

    // Link curvature tool visibility to property
    this.showCurvatureToolProperty.link((show: boolean) => {
      this.curvatureToolContainer.visible = show;
      this.curvatureLabel.visible = show;
      if (show) {
        this.updateCurvatureTool();
      }
    });

    // Update curvature tool when marker position changes
    this.curvatureXProperty.link(() => {
      if (this.showCurvatureToolProperty.value) {
        this.updateCurvatureTool();
      }
    });

    // Update area tool when marker positions change
    this.leftMarkerXProperty.link(() => {
      if (this.showAreaToolProperty.value) {
        this.updateAreaTool();
      }
    });

    this.rightMarkerXProperty.link(() => {
      if (this.showAreaToolProperty.value) {
        this.updateAreaTool();
      }
    });

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
    if (hasWellOffset(this.model)) {
      this.model.wellOffsetProperty.link(() => this.update());
    }
    this.model.particleMassProperty.link(() => this.update());

    // Link to wellSeparationProperty if available (TwoWellsModel and ManyWellsModel)
    if (hasWellSeparation(this.model)) {
      this.model.wellSeparationProperty.link(() => this.update());
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
        show && this.getEffectiveDisplayMode() === "waveFunction";
    });
    this.model.showImaginaryPartProperty.link((show: boolean) => {
      this.imaginaryPartPath.visible =
        show && this.getEffectiveDisplayMode() === "waveFunction";
    });
    this.model.showMagnitudeProperty.link((show: boolean) => {
      this.magnitudePath.visible =
        show && this.getEffectiveDisplayMode() === "waveFunction";
    });

    // Update visibility of classical probability (only for OneWellModel)
    if ("showClassicalProbabilityProperty" in this.model) {
      this.model.showClassicalProbabilityProperty.link(() => {
        this.update();
      });
    }

    // Update visibility of zeros
    this.model.showZerosProperty.link(() => {
      this.update();
    });

    // Update area tool when display mode changes or when data updates
    this.model.displayModeProperty.link(() => {
      if (this.showAreaToolProperty.value) {
        this.updateAreaTool();
      }
    });

    // Update area tool when time changes (for superposition states)
    this.model.timeProperty.link(() => {
      if (this.showAreaToolProperty.value && this.model.isPlayingProperty.value) {
        this.updateAreaTool();
      }
    });

    // Update area tool when potential or well parameters change
    this.model.potentialTypeProperty.link(() => {
      if (this.showAreaToolProperty.value) {
        this.updateAreaTool();
      }
    });

    this.model.wellWidthProperty.link(() => {
      if (this.showAreaToolProperty.value) {
        this.updateAreaTool();
      }
    });

    this.model.wellDepthProperty.link(() => {
      if (this.showAreaToolProperty.value) {
        this.updateAreaTool();
      }
    });

    if ("barrierHeightProperty" in this.model) {
      this.model.barrierHeightProperty.link(() => {
        if (this.showAreaToolProperty.value) {
          this.updateAreaTool();
        }
      });
    }

    if ("potentialOffsetProperty" in this.model) {
      this.model.potentialOffsetProperty.link(() => {
        if (this.showAreaToolProperty.value) {
          this.updateAreaTool();
        }
      });
    }

    // Update area tool when selected energy level changes
    this.model.selectedEnergyLevelIndexProperty.link(() => {
      if (this.showAreaToolProperty.value) {
        this.updateAreaTool();
      }
      if (this.showCurvatureToolProperty.value) {
        this.updateCurvatureTool();
      }
    });

    // Update curvature tool when potential or well parameters change
    this.model.potentialTypeProperty.link(() => {
      if (this.showCurvatureToolProperty.value) {
        this.updateCurvatureTool();
      }
    });

    this.model.wellWidthProperty.link(() => {
      if (this.showCurvatureToolProperty.value) {
        this.updateCurvatureTool();
      }
    });

    this.model.wellDepthProperty.link(() => {
      if (this.showCurvatureToolProperty.value) {
        this.updateCurvatureTool();
      }
    });

    if ("barrierHeightProperty" in this.model) {
      this.model.barrierHeightProperty.link(() => {
        if (this.showCurvatureToolProperty.value) {
          this.updateCurvatureTool();
        }
      });
    }

    if ("potentialOffsetProperty" in this.model) {
      this.model.potentialOffsetProperty.link(() => {
        if (this.showCurvatureToolProperty.value) {
          this.updateCurvatureTool();
        }
      });
    }

    if ("wellSeparationProperty" in this.model) {
      this.model.wellSeparationProperty.link(() => {
        if (this.showCurvatureToolProperty.value) {
          this.updateCurvatureTool();
        }
      });
    }

    if ("numberOfWellsProperty" in this.model) {
      this.model.numberOfWellsProperty.link(() => {
        if (this.showCurvatureToolProperty.value) {
          this.updateCurvatureTool();
        }
      });
    }

    if ("electricFieldProperty" in this.model) {
      this.model.electricFieldProperty.link(() => {
        if (this.showCurvatureToolProperty.value) {
          this.updateCurvatureTool();
        }
      });
    }

    // Initialize labels (important for fixed display mode charts)
    this.updateYAxisLabel();
    this.updateStateLabel();
  }

  /**
   * Gets the effective display mode, using the fixed display mode if set,
   * otherwise falling back to the model's display mode.
   */
  private getEffectiveDisplayMode(): string {
    return this.fixedDisplayMode || this.model.displayModeProperty.value;
  }

  /**
   * Updates the Y-axis label based on display mode.
   */
  private updateYAxisLabel(): void {
    const displayMode = this.getEffectiveDisplayMode();
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
    const displayMode = this.getEffectiveDisplayMode();
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
          // Hide classical probability visualization for superpositions
          this.hideClassicalProbabilityVisualization();
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
          this.updateClassicalProbabilityVisualization(
            boundStates,
            selectedIndex,
          );
        }
      } while (this.updatePending);
    } finally {
      this.isUpdating = false;
      this.updatePending = false;
    }
  }

  /**
   * Hides the classical probability visualization.
   */
  private hideClassicalProbabilityVisualization(): void {
    this.leftTurningPointLine.visible = false;
    this.rightTurningPointLine.visible = false;
    this.leftForbiddenRegion.visible = false;
    this.rightForbiddenRegion.visible = false;
  }

  /**
   * Updates the classical probability visualization (turning points and forbidden regions).
   */
  private updateClassicalProbabilityVisualization(
    boundStates: BoundStateResult,
    selectedIndex: number,
  ): void {
    // Early return if conditions aren't met
    if (
      !this.model.showClassicalProbabilityProperty.value ||
      !hasClassicalTurningPoints(this.model) ||
      selectedIndex < 0 ||
      selectedIndex >= boundStates.energies.length
    ) {
      this.hideClassicalProbabilityVisualization();
      return;
    }

    // TypeScript now knows this.model has getClassicalTurningPoints
    const turningPoints = this.model.getClassicalTurningPoints(selectedIndex);

    if (!turningPoints) {
      this.hideClassicalProbabilityVisualization();
      return;
    }

    // Draw vertical lines at turning points
    const xLeft = this.dataToViewX(turningPoints.left);
    const xRight = this.dataToViewX(turningPoints.right);
    const yTop = this.chartMargins.top;
    const yBottom = this.chartMargins.top + this.plotHeight;

    this.leftTurningPointLine.setLine(xLeft, yTop, xLeft, yBottom);
    this.leftTurningPointLine.visible = true;

    this.rightTurningPointLine.setLine(xRight, yTop, xRight, yBottom);
    this.rightTurningPointLine.visible = true;

    // Draw shaded forbidden regions
    // Left forbidden region (from left edge to left turning point)
    const leftRegionX = this.chartMargins.left;
    const leftRegionWidth = xLeft - leftRegionX;
    this.leftForbiddenRegion.setRect(
      leftRegionX,
      yTop,
      leftRegionWidth,
      this.plotHeight,
    );
    this.leftForbiddenRegion.visible = true;

    // Right forbidden region (from right turning point to right edge)
    const rightRegionX = xRight;
    const rightRegionWidth = this.chartMargins.left + this.plotWidth - xRight;
    this.rightForbiddenRegion.setRect(
      rightRegionX,
      yTop,
      rightRegionWidth,
      this.plotHeight,
    );
    this.rightForbiddenRegion.visible = true;
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
      const displayMode = this.getEffectiveDisplayMode();

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
    if (!boundStates || !hasSuperpositionConfig(this.model)) {
      return;
    }

    // Guard against empty wavefunctions
    if (boundStates.wavefunctions.length === 0) {
      console.warn("No bound states available for superposition view range");
      return;
    }

    const config = this.model.superpositionConfigProperty.value;
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

    const displayMode = this.getEffectiveDisplayMode();

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
    if (!boundStates || !hasSuperpositionConfig(this.model)) {
      return;
    }

    const config = this.model.superpositionConfigProperty.value;
    const time = this.model.timeProperty.value * 1e-15; // Convert fs to seconds
    const displayMode = this.getEffectiveDisplayMode();
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
      this.classicalProbabilityPath.visible = false; // Hide classical for superposition
      this.realPartPath.visible = false;
      this.imaginaryPartPath.visible = false;
      this.magnitudePath.visible = false;
      this.phaseColorNode.visible = false;
    } else if (displayMode === "phaseColor") {
      // Magnitude and phase for coloring
      this.plotPhaseColoredSuperposition(boundStates.xGrid, realPart, imagPart);
      this.probabilityDensityPath.visible = false;
      this.classicalProbabilityPath.visible = false;
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
      this.classicalProbabilityPath.visible = false;
      this.realPartPath.visible = this.model.showRealPartProperty.value;
      this.imaginaryPartPath.visible =
        this.model.showImaginaryPartProperty.value;
      this.magnitudePath.visible = this.model.showMagnitudeProperty.value;
      this.phaseColorNode.visible = false;
    }

    // Update zeros visualization using the real part of the superposition
    this.updateZerosVisualization(boundStates.xGrid, realPart);
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
    const displayMode = this.getEffectiveDisplayMode();
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

      // Plot classical probability if enabled
      if (
        "showClassicalProbabilityProperty" in this.model &&
        this.model.showClassicalProbabilityProperty.value
      ) {
        const classicalProbability =
          this.model.getClassicalProbabilityDensity(selectedIndex);
        if (classicalProbability) {
          this.plotClassicalProbabilityDensity(xGrid, classicalProbability);
          this.classicalProbabilityPath.visible = true;
        } else {
          this.classicalProbabilityPath.visible = false;
        }
      } else {
        this.classicalProbabilityPath.visible = false;
      }
    } else if (displayMode === "phaseColor") {
      // Plot magnitude with phase-colored fill
      this.plotPhaseColoredWavefunction(xGrid, wavefunction, phase);
      this.probabilityDensityPath.visible = false;
      this.classicalProbabilityPath.visible = false;
      this.realPartPath.visible = false;
      this.imaginaryPartPath.visible = false;
      this.magnitudePath.visible = false;
      this.phaseColorNode.visible = true;
    } else {
      // Plot wave function components
      this.plotWaveFunctionComponents(xGrid, wavefunction, phase);
      this.probabilityDensityPath.visible = false;
      this.classicalProbabilityPath.visible = false;
      this.realPartPath.visible = this.model.showRealPartProperty.value;
      this.imaginaryPartPath.visible =
        this.model.showImaginaryPartProperty.value;
      this.magnitudePath.visible = this.model.showMagnitudeProperty.value;
      this.phaseColorNode.visible = false;
    }

    // Update zeros visualization (works for all display modes)
    this.updateZerosVisualization(xGrid, wavefunction);
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
   * Plots the classical probability density with smooth curves.
   */
  private plotClassicalProbabilityDensity(
    xGrid: number[],
    classicalProbability: number[],
  ): void {
    const shape = new Shape();

    // Build points array
    const points: { x: number; y: number }[] = [];
    for (let i = 0; i < xGrid.length; i++) {
      const x = this.dataToViewX(xGrid[i] * QuantumConstants.M_TO_NM);
      const y = this.dataToViewY(classicalProbability[i]);
      points.push({ x, y });
    }

    // Draw  curve
    if (points.length > 0) {
      shape.moveTo(points[0].x, points[0].y);

      for (let i = 0; i < points.length - 1; i++) {
        shape.lineTo(points[i].x, points[i].y);
      }
    }

    this.classicalProbabilityPath.shape = shape;
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
   * Finds zeros (nodes) in the wavefunction data.
   * Zeros occur where the wavefunction changes sign.
   * @param xGrid - Array of x positions in meters
   * @param wavefunction - Array of wavefunction values
   * @returns Array of x positions (in nm) where zeros occur
   */
  private findZeros(xGrid: number[], wavefunction: number[]): number[] {
    const zeros: number[] = [];

    for (let i = 0; i < wavefunction.length - 1; i++) {
      const y1 = wavefunction[i];
      const y2 = wavefunction[i + 1];

      // Check for sign change (zero crossing)
      if (y1 * y2 < 0) {
        // Linear interpolation to find more accurate zero position
        const x1 = xGrid[i];
        const x2 = xGrid[i + 1];
        const zeroX = x1 - y1 * (x2 - x1) / (y2 - y1);
        zeros.push(zeroX * QuantumConstants.M_TO_NM); // Convert to nm
      }
    }

    return zeros;
  }

  /**
   * Updates the visualization of wavefunction zeros (nodes).
   * @param xGrid - Array of x positions in meters
   * @param wavefunction - Array of wavefunction values (real part for complex wavefunctions)
   */
  private updateZerosVisualization(xGrid: number[], wavefunction: number[]): void {
    // Clear existing zeros
    this.zerosNode.removeAllChildren();

    // Only show if enabled
    if (!this.model.showZerosProperty.value) {
      this.zerosNode.visible = false;
      return;
    }

    // Find zeros
    const zeros = this.findZeros(xGrid, wavefunction);

    // Create circles at each zero position
    zeros.forEach(zeroX => {
      const x = this.dataToViewX(zeroX);
      const y = this.dataToViewY(0); // Zeros are at y=0

      const circle = new Circle(4, {
        fill: QPPWColors.energyLevelSelectedProperty,
        stroke: QPPWColors.backgroundColorProperty,
        lineWidth: 1.5,
        centerX: x,
        centerY: y,
      });

      this.zerosNode.addChild(circle);
    });

    this.zerosNode.visible = zeros.length > 0;
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

  /**
   * Converts view X coordinate to data X coordinate using ChartTransform.
   */
  private viewToDataX(x: number): number {
    return this.chartTransform.viewToModelX(x - this.chartMargins.left);
  }

  /**
   * Sets up drag listeners for the area measurement tool markers.
   */
  private setupAreaToolDragListeners(): void {
    // Left marker drag listener
    this.leftMarkerHandle.addInputListener(
      new DragListener({
        drag: (event) => {
          const parentPoint = this.globalToLocalPoint(event.pointer.point);
          let newX = this.viewToDataX(parentPoint.x);

          // Constrain to chart bounds and ensure left marker stays left of right marker
          newX = Math.max(this.xMinProperty.value, newX);
          newX = Math.min(this.rightMarkerXProperty.value - 0.1, newX);

          this.leftMarkerXProperty.value = newX;
        },
      }),
    );

    // Right marker drag listener
    this.rightMarkerHandle.addInputListener(
      new DragListener({
        drag: (event) => {
          const parentPoint = this.globalToLocalPoint(event.pointer.point);
          let newX = this.viewToDataX(parentPoint.x);

          // Constrain to chart bounds and ensure right marker stays right of left marker
          newX = Math.min(this.xMaxProperty.value, newX);
          newX = Math.max(this.leftMarkerXProperty.value + 0.1, newX);

          this.rightMarkerXProperty.value = newX;
        },
      }),
    );
  }

  /**
   * Updates the area measurement tool visualization and calculates the probability.
   */
  private updateAreaTool(): void {
    if (!this.showAreaToolProperty.value) {
      return;
    }

    const leftX = this.leftMarkerXProperty.value;
    const rightX = this.rightMarkerXProperty.value;

    // Convert to view coordinates
    const leftViewX = this.dataToViewX(leftX);
    const rightViewX = this.dataToViewX(rightX);
    const yTop = this.chartMargins.top;
    const yBottom = this.chartMargins.top + this.plotHeight;

    // Update marker lines
    this.leftMarker.setLine(leftViewX, yTop, leftViewX, yBottom);
    this.rightMarker.setLine(rightViewX, yTop, rightViewX, yBottom);

    // Update marker handles (circles at top)
    this.leftMarkerHandle.centerX = leftViewX;
    this.leftMarkerHandle.centerY = yTop + 15;

    this.rightMarkerHandle.centerX = rightViewX;
    this.rightMarkerHandle.centerY = yTop + 15;

    // Update faint background rectangle (shows full measurement region)
    this.areaBackgroundRegion.setRect(
      leftViewX,
      yTop,
      rightViewX - leftViewX,
      this.plotHeight,
    );

    // Create shape that follows the probability density curve
    const shape = this.createAreaShape(leftX, rightX);
    this.areaRegion.shape = shape;

    // Calculate probability in the selected region
    const probability = this.calculateProbabilityInRegion(leftX, rightX);

    // Update label
    if (probability !== null) {
      this.areaLabel.string = `${probability.toFixed(1)}%`;
      // Position label at the center between markers, near the top
      this.areaLabel.centerX = (leftViewX + rightViewX) / 2;
      this.areaLabel.top = this.chartMargins.top + 35;
    } else {
      this.areaLabel.string = "N/A";
    }
  }

  /**
   * Creates a shape that highlights the area under the probability density curve
   * between the specified x positions.
   * @param xStart - Start x position in nm
   * @param xEnd - End x position in nm
   * @returns Shape representing the area under the curve
   */
  private createAreaShape(xStart: number, xEnd: number): Shape {
    const shape = new Shape();

    const boundStates = this.model.getBoundStates();
    if (!boundStates) {
      return shape; // Return empty shape
    }

    // Only show area for probability density mode
    const displayMode = this.getEffectiveDisplayMode();
    if (displayMode !== "probabilityDensity") {
      return shape; // Return empty shape
    }

    // Get probability density data
    const superpositionType = this.model.superpositionTypeProperty.value;
    const isSuperposition = superpositionType !== SuperpositionType.SINGLE;

    let probabilityDensity: number[];

    if (isSuperposition && hasSuperpositionConfig(this.model)) {
      // Calculate superposition probability density
      const config = this.model.superpositionConfigProperty.value;
      const time = this.model.timeProperty.value * 1e-15; // Convert fs to seconds
      const numPoints = boundStates.xGrid.length;

      // Compute time-evolved superposition
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

        // Complex coefficient
        const realCoeff = amplitude * Math.cos(totalPhase);
        const imagCoeff = amplitude * Math.sin(totalPhase);

        // Add contribution to superposition
        for (let i = 0; i < numPoints; i++) {
          realPart[i] += realCoeff * eigenfunction[i];
          imagPart[i] += imagCoeff * eigenfunction[i];
        }
      }

      // Calculate |ψ|²
      probabilityDensity = new Array(numPoints);
      for (let i = 0; i < numPoints; i++) {
        probabilityDensity[i] =
          realPart[i] * realPart[i] + imagPart[i] * imagPart[i];
      }
    } else {
      // Single eigenstate
      const selectedIndex = this.model.selectedEnergyLevelIndexProperty.value;
      if (
        selectedIndex < 0 ||
        selectedIndex >= boundStates.wavefunctions.length
      ) {
        return shape; // Return empty shape
      }

      const wavefunction = boundStates.wavefunctions[selectedIndex];
      probabilityDensity = wavefunction.map((psi) => psi * psi);
    }

    // Build points array for the region
    const xGrid = boundStates.xGrid;
    const points: { x: number; y: number }[] = [];

    // Get the y-coordinate for y=0 (baseline)
    const y0 = this.dataToViewY(0);

    // Add points within the region
    for (let i = 0; i < xGrid.length; i++) {
      const xData = xGrid[i] * QuantumConstants.M_TO_NM; // Convert to nm

      if (xData >= xStart && xData <= xEnd) {
        const x = this.dataToViewX(xData);
        const y = this.dataToViewY(probabilityDensity[i]);
        points.push({ x, y });
      }
    }

    if (points.length === 0) {
      return shape; // Return empty shape if no points in region
    }

    // Create the shape: start at baseline on left, trace curve, return to baseline on right
    const leftViewX = this.dataToViewX(xStart);
    const rightViewX = this.dataToViewX(xEnd);

    // Start at bottom-left (baseline at left marker)
    shape.moveTo(leftViewX, y0);

    // Draw line up to the first point's y-coordinate
    shape.lineTo(leftViewX, points[0].y);

    // Trace along the curve
    for (let i = 0; i < points.length; i++) {
      shape.lineTo(points[i].x, points[i].y);
    }

    // Draw line down from last point to baseline
    shape.lineTo(rightViewX, points[points.length - 1].y);
    shape.lineTo(rightViewX, y0);

    // Close the shape back to starting point
    shape.close();

    return shape;
  }

  /**
   * Calculates the probability (area under the curve) in the specified region.
   * Uses trapezoidal integration of the probability density.
   * @param xStart - Start x position in nm
   * @param xEnd - End x position in nm
   * @returns Probability as a percentage (0-100), or null if not available
   */
  private calculateProbabilityInRegion(
    xStart: number,
    xEnd: number,
  ): number | null {
    const boundStates = this.model.getBoundStates();
    if (!boundStates) {
      return null;
    }

    // Only calculate for probability density mode
    const displayMode = this.getEffectiveDisplayMode();
    if (displayMode !== "probabilityDensity") {
      return null;
    }

    // Check if we're displaying a superposition or a single eigenstate
    const superpositionType = this.model.superpositionTypeProperty.value;
    const isSuperposition = superpositionType !== SuperpositionType.SINGLE;

    let probabilityDensity: number[];

    if (isSuperposition && hasSuperpositionConfig(this.model)) {
      // Calculate superposition probability density
      const config = this.model.superpositionConfigProperty.value;
      const time = this.model.timeProperty.value * 1e-15; // Convert fs to seconds
      const numPoints = boundStates.xGrid.length;

      // Compute time-evolved superposition
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

        // Complex coefficient
        const realCoeff = amplitude * Math.cos(totalPhase);
        const imagCoeff = amplitude * Math.sin(totalPhase);

        // Add contribution to superposition
        for (let i = 0; i < numPoints; i++) {
          realPart[i] += realCoeff * eigenfunction[i];
          imagPart[i] += imagCoeff * eigenfunction[i];
        }
      }

      // Calculate |ψ|²
      probabilityDensity = new Array(numPoints);
      for (let i = 0; i < numPoints; i++) {
        probabilityDensity[i] =
          realPart[i] * realPart[i] + imagPart[i] * imagPart[i];
      }
    } else {
      // Single eigenstate
      const selectedIndex = this.model.selectedEnergyLevelIndexProperty.value;
      if (
        selectedIndex < 0 ||
        selectedIndex >= boundStates.wavefunctions.length
      ) {
        return null;
      }

      const wavefunction = boundStates.wavefunctions[selectedIndex];
      probabilityDensity = wavefunction.map((psi) => psi * psi);
    }

    // Integrate probability density in the region using trapezoidal rule
    const xGrid = boundStates.xGrid;
    let probability = 0;

    for (let i = 0; i < xGrid.length - 1; i++) {
      const x1 = xGrid[i] * QuantumConstants.M_TO_NM; // Convert to nm
      const x2 = xGrid[i + 1] * QuantumConstants.M_TO_NM;

      // Check if this segment overlaps with our region
      if (x2 >= xStart && x1 <= xEnd) {
        // Calculate overlap
        const segmentStart = Math.max(x1, xStart);
        const segmentEnd = Math.min(x2, xEnd);

        if (segmentEnd > segmentStart) {
          // Linearly interpolate probability density values at segment boundaries
          const t1 = (segmentStart - x1) / (x2 - x1);
          const t2 = (segmentEnd - x1) / (x2 - x1);

          const p1 =
            probabilityDensity[i] * (1 - t1) +
            probabilityDensity[i + 1] * t1;
          const p2 =
            probabilityDensity[i] * (1 - t2) +
            probabilityDensity[i + 1] * t2;

          // Trapezoidal rule
          const dx = segmentEnd - segmentStart;
          probability += 0.5 * (p1 + p2) * dx * QuantumConstants.NM_TO_M;
        }
      }
    }

    // Convert to percentage
    return probability * 100;
  }

  /**
   * Sets up the drag listener for the curvature tool marker.
   */
  private setupCurvatureToolDragListener(): void {
    const dragListener = new DragListener({
      drag: (event) => {
        const parentPoint = this.globalToParentPoint(event.pointer.point);
        const localX = parentPoint.x - this.chartMargins.left;
        const dataX = this.viewToDataX(localX);

        // Clamp to chart bounds
        const clampedX = Math.max(
          -X_AXIS_RANGE_NM,
          Math.min(X_AXIS_RANGE_NM, dataX),
        );
        this.curvatureXProperty.value = clampedX;
      },
    });

    this.curvatureMarkerHandle.addInputListener(dragListener);
  }

  /**
   * Updates the curvature tool visualization.
   */
  private updateCurvatureTool(): void {
    if (!this.showCurvatureToolProperty.value) {
      return;
    }

    const x = this.curvatureXProperty.value;
    const viewX = this.dataToViewX(x);
    const yTop = this.chartMargins.top;
    const yBottom = this.chartMargins.top + this.plotHeight;

    // Update marker line
    this.curvatureMarker.setLine(viewX, yTop, viewX, yBottom);

    // Update marker handle
    this.curvatureMarkerHandle.centerX = viewX;
    this.curvatureMarkerHandle.centerY = yTop + 20;

    // Calculate and display derivatives
    const derivatives = this.calculateDerivatives(x);

    if (derivatives !== null) {
      // Create parabola shape with first and second derivatives
      const parabolaShape = this.createCurvatureParabola(
        x,
        derivatives.firstDerivative,
        derivatives.secondDerivative,
        derivatives.wavefunctionValue,
      );
      this.curvatureParabola.shape = parabolaShape;

      // Update label with proper units (nm^-5/2)
      this.curvatureLabel.string = `d²ψ/dx² = ${derivatives.secondDerivative.toExponential(2)} nm⁻⁵ᐟ²`;
      this.curvatureLabel.centerX = viewX;
      this.curvatureLabel.bottom = yTop - 5;
    } else {
      this.curvatureParabola.shape = null;
      this.curvatureLabel.string = "N/A";
    }
  }

  /**
   * Calculates the first and second derivatives of the wavefunction at a given x position.
   * Uses the analytical solution from the model when available for better accuracy,
   * falls back to finite difference method otherwise.
   * @param xData - X position in nm
   * @returns Object with first derivative, second derivative, and wavefunction value, or null if not available
   */
  private calculateDerivatives(
    xData: number,
  ): {
    firstDerivative: number;
    secondDerivative: number;
    wavefunctionValue: number;
  } | null {
    const boundStates = this.model.getBoundStates();
    if (!boundStates) {
      return null;
    }

    // Only calculate for wavefunction mode
    const displayMode = this.getEffectiveDisplayMode();
    if (displayMode !== "waveFunction") {
      return null;
    }

    const selectedIndex = this.model.selectedEnergyLevelIndexProperty.value;
    if (
      selectedIndex < 0 ||
      selectedIndex >= boundStates.wavefunctions.length
    ) {
      return null;
    }

    const xGrid = boundStates.xGrid;
    const wavefunction = boundStates.wavefunctions[selectedIndex];

    // Convert xData from nm to m
    const xInM = xData * QuantumConstants.NM_TO_M;

    // Find the grid points surrounding xData for interpolation
    let i1 = -1;
    for (let i = 0; i < xGrid.length - 1; i++) {
      if (xGrid[i] <= xInM && xInM <= xGrid[i + 1]) {
        i1 = i;
        break;
      }
    }

    if (i1 === -1) {
      return null;
    }

    // Need at least one point on each side for derivatives
    if (i1 === 0 || i1 >= xGrid.length - 2) {
      return null;
    }

    // Interpolate wavefunction value at xData
    const t = (xInM - xGrid[i1]) / (xGrid[i1 + 1] - xGrid[i1]);
    const psi = wavefunction[i1] * (1 - t) + wavefunction[i1 + 1] * t;

    // Grid spacing (assumes uniform grid)
    const h = xGrid[1] - xGrid[0];

    // Calculate first derivative using central difference
    // f'(x) ≈ (f(x+h) - f(x-h)) / (2h)
    const psi_left = wavefunction[i1];
    const psi_right = wavefunction[i1 + 1];
    const firstDerivativeInM = (psi_right - psi_left) / (2 * h);

    // Convert first derivative from m^(-3/2) to nm^(-3/2)
    // dψ/dx in nm^(-3/2) = dψ/dx in m^(-3/2) * (m/nm)
    const firstDerivativeInNm =
      firstDerivativeInM * QuantumConstants.M_TO_NM;

    // Try to use analytical solution for second derivative (more accurate)
    const secondDerivativeArray = this.model.getWavefunctionSecondDerivative(
      selectedIndex + 1,
      [xInM],
    );

    let secondDerivativeInNm: number;

    if (secondDerivativeArray && secondDerivativeArray.length > 0) {
      // Analytical solution available for second derivative
      const secondDerivativeInM = secondDerivativeArray[0];

      // Convert from m^(-5/2) to nm^(-5/2)
      secondDerivativeInNm =
        secondDerivativeInM * Math.pow(QuantumConstants.M_TO_NM, 2);
    } else {
      // Fall back to finite difference for second derivative
      // f''(x) ≈ (f(x-h) - 2f(x) + f(x+h)) / h²
      const secondDerivativeInM = (psi_left - 2 * psi + psi_right) / (h * h);

      // Convert from m^(-5/2) to nm^(-5/2)
      secondDerivativeInNm =
        secondDerivativeInM * Math.pow(QuantumConstants.M_TO_NM, 2);
    }

    return {
      firstDerivative: firstDerivativeInNm,
      secondDerivative: secondDerivativeInNm,
      wavefunctionValue: psi,
    };
  }

  /**
   * Creates a parabola shape that matches the curvature at a point.
   * Uses Taylor expansion: y = y₀ + y'₀(x - x₀) + (1/2)y''₀(x - x₀)²
   * @param xData - X position in nm
   * @param firstDerivative - First derivative value at the point (nm^-3/2)
   * @param secondDerivative - Second derivative value at the point (nm^-5/2)
   * @param wavefunctionValue - Value of wavefunction at the point (nm^-1/2)
   * @returns Shape representing the parabola
   */
  private createCurvatureParabola(
    xData: number,
    firstDerivative: number,
    secondDerivative: number,
    wavefunctionValue: number,
  ): Shape {
    const shape = new Shape();

    // Create a small parabola centered at (xData, wavefunctionValue)
    // Using Taylor expansion: y = y₀ + y'₀ * (x - x₀) + (1/2) * y''₀ * (x - x₀)²
    const centerX = xData;
    const centerY = wavefunctionValue;

    // Width of parabola in nm (±0.3 nm from center)
    const parabolaHalfWidth = 0.3;

    const points: { x: number; y: number }[] = [];

    // Generate parabola points
    for (
      let dx = -parabolaHalfWidth;
      dx <= parabolaHalfWidth;
      dx += parabolaHalfWidth / 10
    ) {
      const x = centerX + dx;
      // Taylor expansion with first and second derivative terms
      const y =
        centerY + firstDerivative * dx + 0.5 * secondDerivative * dx * dx;

      const viewX = this.dataToViewX(x);
      const viewY = this.dataToViewY(y);

      points.push({ x: viewX, y: viewY });
    }

    if (points.length > 0) {
      shape.moveTo(points[0].x, points[0].y);
      for (let i = 1; i < points.length; i++) {
        shape.lineTo(points[i].x, points[i].y);
      }
    }

    return shape;
  }
}
