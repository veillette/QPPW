/**
 * EnergyChartNode displays the potential energy landscape and energy levels.
 * This is the top chart in the One Well screen.
 */

import {
  Rectangle,
  Line,
  Path,
  Text,
  VBox,
  SceneryEvent,
} from "scenerystack/scenery";
import { Shape } from "scenerystack/kite";
import { Range } from "scenerystack/dot";
import { Orientation } from "scenerystack/phet-core";
import { Checkbox } from "scenerystack/sun";
import { AxisLine, TickMarkSet, TickLabelSet } from "scenerystack/bamboo";
import { DerivedProperty } from "scenerystack/axon";
import { BaseChartNode, ChartOptions } from "./BaseChartNode.js";
import type { ScreenModel } from "../model/ScreenModels.js";
import {
  isTwoWellsModel,
  isManyWellsModel,
  hasBarrierHeight,
  hasPotentialOffset,
  hasWellSeparation,
  hasElectricField,
  hasClassicalTurningPoints,
} from "../model/ModelTypeGuards.js";
import { PotentialType, BoundStateResult } from "../model/PotentialFunction.js";
import QuantumConstants from "../model/QuantumConstants.js";
import QPPWColors from "../../QPPWColors.js";
import { PhetFont } from "scenerystack/scenery-phet";
import stringManager from "../../i18n/StringManager.js";
import type { ScreenViewState } from "./ScreenViewStates.js";
import { Node } from "scenerystack/scenery";
import { PotentialRenderer, RenderContext } from "./potential-renderers/PotentialRenderer.js";
import { InfiniteWellRenderer } from "./potential-renderers/InfiniteWellRenderer.js";
import { FiniteWellRenderer } from "./potential-renderers/FiniteWellRenderer.js";
import { HarmonicOscillatorRenderer } from "./potential-renderers/HarmonicOscillatorRenderer.js";
import { AsymmetricTriangleRenderer } from "./potential-renderers/AsymmetricTriangleRenderer.js";
import { TriangularRenderer } from "./potential-renderers/TriangularRenderer.js";
import { CoulombRenderer } from "./potential-renderers/CoulombRenderer.js";
import { DoubleSquareWellRenderer } from "./potential-renderers/DoubleSquareWellRenderer.js";
import { MorseRenderer } from "./potential-renderers/MorseRenderer.js";
import { PoschlTellerRenderer } from "./potential-renderers/PoschlTellerRenderer.js";
import { RosenMorseRenderer } from "./potential-renderers/RosenMorseRenderer.js";
import { EckartRenderer } from "./potential-renderers/EckartRenderer.js";
import { MultiSquareWellRenderer } from "./potential-renderers/MultiSquareWellRenderer.js";
import { MultiCoulomb1DRenderer } from "./potential-renderers/MultiCoulomb1DRenderer.js";

// Chart axis range constants
const X_AXIS_RANGE_NM = 4; // X-axis extends from -X_AXIS_RANGE_NM to +X_AXIS_RANGE_NM

// Energy axis ranges depend on potential type
// For potentials with V=0 at center (e.g., Harmonic Oscillator, Infinite Well): -5 to 15 eV
// For potentials with V=0 at infinity (e.g., Finite Well, Coulomb): -15 to 5 eV
// For Asymmetric Triangle: -5 to 15 eV (special case with shifted potential)
function getEnergyAxisRange(potentialType: PotentialType): {
  min: number;
  max: number;
} {
  switch (potentialType) {
    case PotentialType.HARMONIC_OSCILLATOR:
      // V=0 at center, grows outward
      return { min: -5, max: 15 };
    case PotentialType.ASYMMETRIC_TRIANGLE:
      // Special case: V=10eV at infinity, V=0 at x=0
      return { min: -5, max: 15 };
    case PotentialType.TRIANGULAR:
      // Triangular well: offset can be -5 to 15 eV, height adds to that
      return { min: -5, max: 15 };
    case PotentialType.INFINITE_WELL:
      // V=0 inside well (centered at x=0), V=∞ outside (displayed as 15 eV)
      return { min: -5, max: 15 };
    case PotentialType.DOUBLE_SQUARE_WELL:
      // Double square well with barrier between wells
      // Energy reference: V=0 in wells, V=wellDepth in barrier
      // Bound states have positive energies (between 0 and wellDepth)
      return { min: 0, max: 20 };
    case PotentialType.MULTI_SQUARE_WELL:
      // Multi-square well (generalization of double square well)
      // Energy reference: V=0 in wells, V=wellDepth in barrier
      return { min: -5, max: 15 };
    case PotentialType.MULTI_COULOMB_1D:
      // Multi-Coulomb 1D potential with V=0 at infinity
      return { min: -15, max: 5 };
    case PotentialType.MORSE:
      // Morse potential: V=0 at dissociation limit (infinity), V=-De at bottom
      return { min: -15, max: 5 };
    case PotentialType.FINITE_WELL:
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

export class EnergyChartNode extends BaseChartNode {
  // Visual elements specific to energy chart
  private readonly potentialPath: Path;
  private readonly energyLevelNodes: Map<number, Line>;
  private readonly energyLabelNodes: Map<number, Text>;
  private readonly totalEnergyLine: Line;
  private readonly legendNode: Node;

  // Classical turning point lines
  private readonly leftTurningPointLine: Line;
  private readonly rightTurningPointLine: Line;

  // Potential renderer registry
  private readonly renderers: Map<PotentialType, PotentialRenderer>;

  // Hover state
  private hoveredEnergyLevelIndex: number | null = null;

  public constructor(
    model: ScreenModel,
    viewState: ScreenViewState,
    options?: { width?: number; height?: number },
  ) {
    // Initialize view range with values based on initial potential type
    const initialEnergyRange = getEnergyAxisRange(
      model.potentialTypeProperty.value,
    );

    // Call super constructor with chart options
    const chartOptions: ChartOptions = {
      width: options?.width ?? 600,
      height: options?.height ?? 300,
      margins: { left: 60, right: 20, top: 40, bottom: 50 },
      xRange: { min: -X_AXIS_RANGE_NM, max: X_AXIS_RANGE_NM },
      yRange: { min: initialEnergyRange.min, max: initialEnergyRange.max },
      showZeroLine: true,
    };

    super(model, viewState, chartOptions);

    // PDOM - make energy chart accessible with dynamic description
    this.tagName = "div";
    this.labelTagName = "h3";
    this.labelContent = "Energy Level Diagram";
    this.descriptionTagName = "p";
    this.descriptionContent = new DerivedProperty(
      [
        model.potentialTypeProperty,
        model.wellWidthProperty,
        model.wellDepthProperty,
        model.selectedEnergyLevelIndexProperty,
      ],
      (
        potentialType: PotentialType,
        width: number,
        depth: number,
        selectedIndex: number,
      ) => {
        return this.createEnergyChartDescription(
          potentialType,
          width,
          depth,
          selectedIndex,
        );
      },
    );

    // Reconfigure zero line with energy chart specific styling
    this.zeroLine.stroke = QPPWColors.potentialBarrierProperty;
    this.zeroLine.lineWidth = 1;
    this.zeroLine.lineDash = [5, 5];

    // Initialize potential renderer registry
    this.renderers = new Map([
      [PotentialType.INFINITE_WELL, new InfiniteWellRenderer()],
      [PotentialType.FINITE_WELL, new FiniteWellRenderer()],
      [PotentialType.HARMONIC_OSCILLATOR, new HarmonicOscillatorRenderer()],
      [PotentialType.ASYMMETRIC_TRIANGLE, new AsymmetricTriangleRenderer()],
      [PotentialType.TRIANGULAR, new TriangularRenderer()],
      [PotentialType.COULOMB_1D, new CoulombRenderer()],
      [PotentialType.COULOMB_3D, new CoulombRenderer()],
      [PotentialType.DOUBLE_SQUARE_WELL, new DoubleSquareWellRenderer()],
      [PotentialType.MORSE, new MorseRenderer()],
      [PotentialType.POSCHL_TELLER, new PoschlTellerRenderer()],
      [PotentialType.ROSEN_MORSE, new RosenMorseRenderer()],
      [PotentialType.ECKART, new EckartRenderer()],
      [PotentialType.MULTI_SQUARE_WELL, new MultiSquareWellRenderer()],
      [PotentialType.MULTI_COULOMB_1D, new MultiCoulomb1DRenderer()],
    ]);

    // Create axes
    this.axesNode = this.createAxes();
    this.addChild(this.axesNode);

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

    // Create legend (outside clipped area)
    this.legendNode = this.createLegend();
    this.addChild(this.legendNode);

    // Link to model properties
    this.linkToModel();

    // Set up keyboard navigation for energy levels
    this.setupKeyboardNavigation();

    // Note: Initial update is now done asynchronously inside linkToModel()
    // to prevent blocking the page load
  }

  /**
   * Creates an accessible description of the energy chart based on current state.
   * This provides screen reader users with meaningful information about the visualization.
   */
  private createEnergyChartDescription(
    potentialType: PotentialType,
    width: number,
    depth: number,
    selectedIndex: number,
  ): string {
    const boundStates = this.model.getBoundStates();
    if (!boundStates || boundStates.energies.length === 0) {
      return "No bound states found for current potential configuration.";
    }

    const numLevels = boundStates.energies.length;
    const energies = boundStates.energies.map(
      (e) => e * QuantumConstants.JOULES_TO_EV,
    );
    const groundEnergy = energies[0];

    let description = `${potentialType} potential well. `;
    description += `Width: ${width.toFixed(2)} nanometers. `;

    // Add depth information if applicable
    if (
      potentialType !== PotentialType.INFINITE_WELL &&
      potentialType !== PotentialType.HARMONIC_OSCILLATOR
    ) {
      description += `Depth: ${depth.toFixed(2)} electron volts. `;
    }

    description += `\n\n`;
    description += `Found ${numLevels} bound state${numLevels !== 1 ? "s" : ""}. `;
    description += `Ground state energy: ${groundEnergy.toFixed(3)} eV. `;

    if (numLevels > 1) {
      const firstExcited = energies[1];
      const spacing = firstExcited - groundEnergy;
      description += `First excited state: ${firstExcited.toFixed(3)} eV. `;
      description += `Energy spacing: ${spacing.toFixed(3)} eV. `;
    }

    // Selected level information
    if (selectedIndex >= 0 && selectedIndex < numLevels) {
      const selectedEnergy = energies[selectedIndex];
      description += `\n\n`;
      description += `Currently viewing level ${selectedIndex + 1} `;
      description += `with energy ${selectedEnergy.toFixed(3)} eV.`;

      // Classical turning points if available
      if (hasClassicalTurningPoints(this.model)) {
        const turningPoints =
          this.model.getClassicalTurningPoints(selectedIndex);
        if (turningPoints) {
          description += ` Classical turning points at `;
          description += `${turningPoints.left.toFixed(2)} nm and `;
          description += `${turningPoints.right.toFixed(2)} nm.`;
        }
      }
    }

    return description;
  }

  /**
   * Sets up keyboard navigation for energy level selection.
   * Arrow keys navigate between levels, Home/End jump to first/last.
   */
  private setupKeyboardNavigation(): void {
    // Add keyboard listener to the entire chart node
    this.addInputListener({
      keydown: (event: SceneryEvent<KeyboardEvent>) => {
        if (!event.domEvent) {
          return;
        }

        const boundStates = this.model.getBoundStates();
        if (!boundStates || boundStates.energies.length === 0) {
          return;
        }

        const currentLevel = this.model.selectedEnergyLevelIndexProperty.value;
        const maxLevel = boundStates.energies.length - 1;
        let newLevel = currentLevel;

        switch (event.domEvent.key) {
          case "ArrowUp":
          case "ArrowRight":
            // Move to higher energy level (higher index)
            if (currentLevel < maxLevel) {
              newLevel = currentLevel + 1;
              event.domEvent.preventDefault();
            }
            break;

          case "ArrowDown":
          case "ArrowLeft":
            // Move to lower energy level (lower index)
            if (currentLevel > 0) {
              newLevel = currentLevel - 1;
              event.domEvent.preventDefault();
            }
            break;

          case "Home":
            // Jump to ground state (level 0)
            newLevel = 0;
            event.domEvent.preventDefault();
            break;

          case "End":
            // Jump to highest energy level
            newLevel = maxLevel;
            event.domEvent.preventDefault();
            break;
        }

        // Update selection if changed
        if (newLevel !== currentLevel) {
          this.model.selectedEnergyLevelIndexProperty.value = newLevel;
        }
      },
    });
  }

  /**
   * Creates the axes (X and Y) with labels - using bamboo components where possible.
   * Note: GridLineSet from bamboo causes infinite loops, so we use manual grid lines.
   */
  protected createAxes(): Node {
    const axesNode = new Node();

    // Manual Y-axis grid lines (GridLineSet causes hang)
    const yMin = this.yMinProperty.value;
    const yMax = this.yMaxProperty.value;
    for (let energy = yMin; energy <= yMax; energy += 5) {
      if (energy !== yMin) {
        const y =
          this.chartMargins.top + this.chartTransform.modelToViewY(energy);
        const gridLine = new Line(
          this.chartMargins.left,
          y,
          this.chartMargins.left + this.plotWidth,
          y,
          {
            stroke: QPPWColors.gridLineProperty,
            lineWidth: 1,
          },
        );
        axesNode.addChild(gridLine);
      }
    }

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

    // X-axis at bottom using bamboo AxisLine (at model y=yMin)
    const xAxisBottomNode = new AxisLine(
      this.chartTransform,
      Orientation.HORIZONTAL,
      {
        stroke: QPPWColors.axisProperty,
        lineWidth: 2,
        value: this.yMinProperty.value,
      },
    );
    xAxisBottomNode.x = this.chartMargins.left;
    xAxisBottomNode.y = this.chartMargins.top;
    axesNode.addChild(xAxisBottomNode);

    // Y-axis tick marks using bamboo TickMarkSet
    const yTickMarksNode = new TickMarkSet(
      this.chartTransform,
      Orientation.VERTICAL,
      5,
      {
        edge: "min",
        extent: 8,
        stroke: QPPWColors.axisProperty,
        lineWidth: 1,
      },
    );
    yTickMarksNode.x = this.chartMargins.left;
    yTickMarksNode.y = this.chartMargins.top;
    axesNode.addChild(yTickMarksNode);

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

    // Y-axis tick labels using bamboo TickLabelSet
    const yTickLabelsNode = new TickLabelSet(
      this.chartTransform,
      Orientation.VERTICAL,
      5,
      {
        edge: "min",
        createLabel: (value: number) =>
          new Text(value.toFixed(0), {
            font: new PhetFont(12),
            fill: QPPWColors.labelFillProperty,
          }),
      },
    );
    yTickLabelsNode.x = this.chartMargins.left;
    yTickLabelsNode.y = this.chartMargins.top;
    axesNode.addChild(yTickLabelsNode);

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

    // Axis labels
    const yLabelText = new Text(stringManager.energyEvStringProperty, {
      font: new PhetFont(14),
      fill: QPPWColors.labelFillProperty,
      rotation: -Math.PI / 2,
      centerX: this.chartMargins.left - 40,
      centerY: this.chartHeight / 2,
    });
    axesNode.addChild(yLabelText);

    const xLabelText = new Text(stringManager.positionNmStringProperty, {
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
    const legendContentVBox = new VBox({
      spacing: 5,
      align: "left",
      children: [
        new Checkbox(
          this.viewState.showTotalEnergyProperty,
          new Text(stringManager.totalEnergyStringProperty, {
            font: "12px sans-serif",
            fill: QPPWColors.textFillProperty,
          }),
          {
            boxWidth: 15,
          },
        ),
        new Checkbox(
          this.viewState.showPotentialEnergyProperty,
          new Text(stringManager.potentialEnergyStringProperty, {
            font: "12px sans-serif",
            fill: QPPWColors.textFillProperty,
          }),
          {
            boxWidth: 15,
          },
        ),
      ],
    });

    const legendPanelRectangle = new Rectangle(
      0,
      0,
      legendContentVBox.width + 20,
      legendContentVBox.height + 15,
      5,
      5,
      {
        fill: QPPWColors.panelFillProperty,
        stroke: QPPWColors.gridLineProperty,
        lineWidth: 1,
      },
    );

    const legendNode = new Node({
      children: [legendPanelRectangle, legendContentVBox],
      left: this.chartMargins.left + 10,
      top: this.chartMargins.top + 10,
    });

    legendContentVBox.left = 10;
    legendContentVBox.top = 8;

    return legendNode;
  }

  /**
   * Links chart updates to model property changes.
   */
  protected linkToModel(): void {
    // Update when any parameter changes
    this.model.potentialTypeProperty.lazyLink(() => {
      this.updateEnergyAxisRange();
      this.update();
    });
    this.model.wellWidthProperty.lazyLink(() => this.update());
    this.model.wellDepthProperty.lazyLink(() => this.update());
    if ("wellOffsetProperty" in this.model) {
      this.model.wellOffsetProperty.lazyLink(() => this.update());
    }
    this.model.particleMassProperty.lazyLink(() => this.update());
    this.model.selectedEnergyLevelIndexProperty.lazyLink(() =>
      this.updateSelection(),
    );
    this.viewState.showTotalEnergyProperty.lazyLink((show: boolean) => {
      this.totalEnergyLine.visible = show;
    });
    this.viewState.showPotentialEnergyProperty.lazyLink((show: boolean) => {
      this.potentialPath.visible = show;
    });

    // Link to wellSeparationProperty if available (TwoWellsModel and ManyWellsModel)
    if (hasWellSeparation(this.model)) {
      this.model.wellSeparationProperty.lazyLink(() => this.update());
    }

    // Link to numberOfWellsProperty and electricFieldProperty if available (ManyWellsModel only)
    if (isManyWellsModel(this.model)) {
      this.model.numberOfWellsProperty.lazyLink(() => this.update());
      if (hasElectricField(this.model)) {
        this.model.electricFieldProperty.lazyLink(() => this.update());
      }
    }

    // Link to barrierHeightProperty and potentialOffsetProperty if available (OneWellModel only)
    if (hasBarrierHeight(this.model)) {
      this.model.barrierHeightProperty.lazyLink(() => this.update());
    }
    if (hasPotentialOffset(this.model)) {
      this.model.potentialOffsetProperty.lazyLink(() => this.update());
    }

    // Update classical probability visualization when property changes
    this.viewState.showClassicalProbabilityProperty.lazyLink(() =>
      this.update(),
    );

    // Perform initial updates asynchronously (after construction completes)
    // This prevents blocking the page load with expensive calculations
    setTimeout(() => {
      this.updateEnergyAxisRange();
      this.update();
      this.updateSelection();
    }, 0);
  }

  /**
   * Updates the energy axis range based on the current potential type.
   */
  private updateEnergyAxisRange(): void {
    const energyRange = getEnergyAxisRange(
      this.model.potentialTypeProperty.value,
    );
    this.yMinProperty.value = energyRange.min;
    this.yMaxProperty.value = energyRange.max;

    // Update the chart transform with new Y range
    this.chartTransform.setModelYRange(
      new Range(energyRange.min, energyRange.max),
    );

    // Recreate axes with new range
    this.removeChild(this.axesNode);
    this.axesNode = this.createAxes();
    this.addChild(this.axesNode);
    this.axesNode.moveToBack();
    this.backgroundRect.moveToBack();
  }

  /**
   * Main update method - recalculates and redraws everything.
   * Called automatically when model properties change, but can also be called explicitly (e.g., during reset).
   */
  public update(): void {
    const boundStates = this.model.getBoundStates();
    if (!boundStates) {
      // Clear the chart if bound states cannot be calculated
      this.clearChart();
      return;
    }

    this.updatePotentialCurve(boundStates);
    this.updateEnergyLevels(boundStates);
    this.updateZeroLine();
    this.updateTotalEnergyLine(boundStates);
    this.updateClassicalTurningPoints(boundStates);
  }

  /**
   * Updates the classical turning point lines.
   */
  private updateClassicalTurningPoints(boundStates: BoundStateResult): void {
    const selectedIndex = this.model.selectedEnergyLevelIndexProperty.value;

    // Early return if conditions aren't met
    if (
      !this.viewState.showClassicalProbabilityProperty.value ||
      !hasClassicalTurningPoints(this.model) ||
      selectedIndex < 0 ||
      selectedIndex >= boundStates.energies.length
    ) {
      this.leftTurningPointLine.visible = false;
      this.rightTurningPointLine.visible = false;
      return;
    }

    // TypeScript now knows this.model has getClassicalTurningPoints
    const turningPoints = this.model.getClassicalTurningPoints(selectedIndex);

    if (!turningPoints) {
      this.leftTurningPointLine.visible = false;
      this.rightTurningPointLine.visible = false;
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
  }

  /**
   * Clears all chart elements (potential curve, energy levels, etc.)
   * Called when bound states cannot be calculated.
   */
  private clearChart(): void {
    // Clear potential curve
    this.potentialPath.shape = null;

    // Clear energy level lines
    this.energyLevelNodes.forEach((line) => {
      this.plotContentNode.removeChild(line);
    });
    this.energyLevelNodes.clear();

    // Clear energy labels
    this.energyLabelNodes.forEach((label) => {
      this.removeChild(label);
    });
    this.energyLabelNodes.clear();

    // Hide total energy line
    this.totalEnergyLine.visible = false;
  }

  /**
   * Updates the potential energy curve using the renderer registry.
   * This method has been refactored to use the Strategy pattern with PotentialRenderer classes.
   * The old 520-line implementation with 13 if/else blocks has been extracted into separate renderer files.
   */
  private updatePotentialCurve(boundStates: BoundStateResult): void {
    // Clear the old shape explicitly
    this.potentialPath.shape = null;

    const xGrid = boundStates.xGrid;
    const potentialType = this.model.potentialTypeProperty.value;
    const wellWidth = this.model.wellWidthProperty.value;
    const wellDepth = this.model.wellDepthProperty.value;

    // Calculate the center of the xGrid for alignment
    const xCenter =
      ((xGrid[0] + xGrid[xGrid.length - 1]) / 2) * QuantumConstants.M_TO_NM;

    // Get the renderer for this potential type
    const renderer = this.renderers.get(potentialType);

    if (renderer) {
      // Create the render context with all necessary information
      const context: RenderContext = {
        xGrid,
        xCenter,
        wellWidth,
        wellDepth,
        dataToViewX: this.dataToViewX.bind(this),
        dataToViewY: this.dataToViewY.bind(this),
        chartMargins: this.chartMargins,
        chartWidth: this.chartWidth,
        chartHeight: this.chartHeight,
        model: this.model,
      };

      // Render the potential using the appropriate renderer
      const shape = renderer.render(context);
      this.potentialPath.shape = shape;
    } else {
      // Fallback for unknown potential types: draw a flat line at y=0
      const shape = new Shape();
      const y0 = this.dataToViewY(0);
      shape.moveTo(this.chartMargins.left, y0);
      shape.lineTo(this.chartWidth - this.chartMargins.right, y0);
      this.potentialPath.shape = shape;
    }
  }

  // [OLD 520-LINE METHOD REMOVED]
  // The previous implementation contained 13 massive if/else blocks that have been
  // extracted into separate renderer classes in: src/common/view/potential-renderers/

  /**
   * Updates just the styling (hover/selection state) of existing energy levels without recalculating positions.
   */
  private updateEnergyLevelStyling(): void {
    this.energyLevelNodes.forEach((line, index) => {
      const isSelected =
        index === this.model.selectedEnergyLevelIndexProperty.value;
      const isHovered = index === this.hoveredEnergyLevelIndex;

      line.stroke = isSelected
        ? QPPWColors.energyLevelSelectedProperty
        : QPPWColors.energyLevelProperty;
      line.lineWidth = isSelected ? 4 : isHovered ? 3 : 2;
      line.opacity = isHovered ? 1 : 0.7;
    });

    this.energyLabelNodes.forEach((label, index) => {
      const isHovered = index === this.hoveredEnergyLevelIndex;
      label.visible = isHovered;
    });
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
    const energies = boundStates.energies.map(
      (e) => e * QuantumConstants.JOULES_TO_EV,
    );
    const x1 = this.chartMargins.left;
    const x2 = this.chartWidth - this.chartMargins.right;

    energies.forEach((energy, index) => {
      const y = this.dataToViewY(energy);
      const isSelected =
        index === this.model.selectedEnergyLevelIndexProperty.value;
      const isHovered = index === this.hoveredEnergyLevelIndex;

      const line = new Line(x1, y, x2, y, {
        stroke: isSelected
          ? QPPWColors.energyLevelSelectedProperty
          : QPPWColors.energyLevelProperty,
        lineWidth: isSelected ? 4 : isHovered ? 3 : 2,
        cursor: "pointer",
        opacity: isHovered ? 1 : 0.7,
      });

      // Create a wider invisible hit area to make the line easier to grab
      const hitAreaHeight = 5; // pixels above and below the line
      const hitArea = new Rectangle(
        x1,
        y - hitAreaHeight,
        x2 - x1,
        hitAreaHeight * 2,
        {
          fill: "transparent",
          cursor: "pointer",

          // PDOM - make energy level selection keyboard accessible
          tagName: "button",
          ariaRole: "radio",
          innerContent: `Level ${index + 1}`,
          accessibleName: `Energy Level ${index + 1}`,
          descriptionContent: `Energy: ${energy.toFixed(3)} electron volts. ${index} node${index !== 1 ? "s" : ""}.`,
          focusable: true,
        },
      );

      // Update aria-checked when selection changes
      const updateAriaChecked = () => {
        const isSelected =
          index === this.model.selectedEnergyLevelIndexProperty.value;
        hitArea.setPDOMAttribute("aria-checked", isSelected.toString());
      };

      // Set initial state
      updateAriaChecked();

      // Update when selection changes
      this.model.selectedEnergyLevelIndexProperty.link(updateAriaChecked);

      // Add click handler to hit area
      hitArea.addInputListener({
        down: () => {
          // Only set if index is within valid range
          if (index >= 0 && index < boundStates.energies.length) {
            this.model.selectedEnergyLevelIndexProperty.value = index;
          }
        },
        click: () => {
          // Also handle click events for mouse users
          if (index >= 0 && index < boundStates.energies.length) {
            this.model.selectedEnergyLevelIndexProperty.value = index;
          }
        },
        enter: () => {
          this.hoveredEnergyLevelIndex = index;
          this.updateEnergyLevelStyling();
        },
        exit: () => {
          this.hoveredEnergyLevelIndex = null;
          this.updateEnergyLevelStyling();
        },
      });

      this.energyLevelNodes.set(index, line);
      this.plotContentNode.addChild(line); // Add to clipped content
      this.plotContentNode.addChild(hitArea); // Add hit area on top

      // Add energy label (outside clipped area for visibility)
      // Use template string from i18n: "E{{level}} = {{value}} eV"
      const template = stringManager.energyLevelLabelStringProperty.value;
      const labelText = template
        .replace("{{level}}", (index + 1).toString())
        .replace("{{value}}", energy.toFixed(3));
      const label = new Text(labelText, {
        font: "10px sans-serif",
        fill: QPPWColors.labelFillProperty,
        left: x2 + 5,
        centerY: y,
        visible: isHovered, // Only show when hovering
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
  protected updateZeroLine(): void {
    const y = this.dataToViewY(0);
    this.zeroLine.x1 = this.chartMargins.left;
    this.zeroLine.y1 = y;
    this.zeroLine.x2 = this.chartWidth - this.chartMargins.right;
    this.zeroLine.y2 = y;

    // Hide zero line for Coulomb potentials to reduce clutter
    const potentialType = this.model.potentialTypeProperty.value;
    this.zeroLine.visible =
      potentialType !== PotentialType.COULOMB_1D &&
      potentialType !== PotentialType.COULOMB_3D;
  }

  /**
   * Updates the total energy line for the selected state.
   */
  private updateTotalEnergyLine(boundStates: BoundStateResult): void {
    const selectedIndex = this.model.selectedEnergyLevelIndexProperty.value;
    if (selectedIndex >= 0 && selectedIndex < boundStates.energies.length) {
      const energy =
        boundStates.energies[selectedIndex] * QuantumConstants.JOULES_TO_EV;
      const y = this.dataToViewY(energy);
      this.totalEnergyLine.x1 = this.chartMargins.left;
      this.totalEnergyLine.y1 = y;
      this.totalEnergyLine.x2 = this.chartWidth - this.chartMargins.right;
      this.totalEnergyLine.y2 = y;
      this.totalEnergyLine.visible =
        this.viewState.showTotalEnergyProperty.value;
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
      this.updateEnergyLevelStyling();
      this.updateTotalEnergyLine(boundStates);
      this.updateClassicalTurningPoints(boundStates);
    }
  }
}
