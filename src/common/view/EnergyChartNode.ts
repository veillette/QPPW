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
import { ChartTransform, ChartRectangle, AxisLine } from "scenerystack/bamboo";
import { OneWellModel } from "../../one-well/model/OneWellModel.js";
import { PotentialType, BoundStateResult } from "../model/PotentialFunction.js";
import QuantumConstants from "../model/QuantumConstants.js";
import QPPWColors from "../../QPPWColors.js";

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
  private readonly potentialPath: Path;
  private readonly energyLevelNodes: Map<number, Line>;
  private readonly totalEnergyLine: Line;
  private readonly zeroLine: Line;
  private readonly axesNode: Node;
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

    // Initialize view range (will be updated based on data)
    this.xMinProperty = new NumberProperty(0);
    this.xMaxProperty = new NumberProperty(1);
    this.yMinProperty = new NumberProperty(-1);
    this.yMaxProperty = new NumberProperty(10);

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

    // Create zero line
    this.zeroLine = new Line(0, 0, 0, 0, {
      stroke: "red",
      lineWidth: 1,
      lineDash: [5, 5],
    });
    this.addChild(this.zeroLine);

    // Create potential energy path
    this.potentialPath = new Path(null, {
      stroke: QPPWColors.potentialWellProperty,
      lineWidth: 3,
    });
    this.addChild(this.potentialPath);

    // Create energy level lines container
    this.energyLevelNodes = new Map();

    // Create total energy line
    this.totalEnergyLine = new Line(0, 0, 0, 0, {
      stroke: QPPWColors.energyLevelProperty,
      lineWidth: 2,
      lineDash: [10, 5],
    });
    this.addChild(this.totalEnergyLine);

    // Create legend
    this.legendNode = this.createLegend();
    this.addChild(this.legendNode);

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

    // Axis labels
    const yLabel = new Text("Energy (eV)", {
      font: "14px sans-serif",
      fill: QPPWColors.labelFillProperty,
      rotation: -Math.PI / 2,
      centerX: this.chartMargins.left - 40,
      centerY: this.chartHeight / 2,
    });
    axesNode.addChild(yLabel);

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
   * Creates the legend for toggling visibility of elements.
   */
  private createLegend(): Node {
    const legendContent = new VBox({
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

    const legendPanel = new Rectangle(0, 0, legendContent.width + 20, legendContent.height + 15, 5, 5, {
      fill: "rgba(255, 255, 255, 0.85)",
      stroke: QPPWColors.gridLineProperty,
      lineWidth: 1,
    });

    const legendNode = new Node({
      children: [legendPanel, legendContent],
      left: this.chartMargins.left + 10,
      top: this.chartMargins.top + 10,
    });

    legendContent.left = 10;
    legendContent.top = 8;

    return legendNode;
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
    this.model.selectedEnergyLevelIndexProperty.link(() => this.updateSelection());
    this.model.showTotalEnergyProperty.link((show) => {
      this.totalEnergyLine.visible = show;
    });
    this.model.showPotentialEnergyProperty.link((show) => {
      this.potentialPath.visible = show;
    });
  }

  /**
   * Main update method - recalculates and redraws everything.
   */
  private update(): void {
    const boundStates = this.model.getBoundStates();
    if (!boundStates) {
      return;
    }

    this.updateViewRange(boundStates);
    this.updatePotentialCurve(boundStates);
    this.updateEnergyLevels(boundStates);
    this.updateZeroLine();
    this.updateTotalEnergyLine(boundStates);
  }

  /**
   * Updates the view range based on the data and updates the ChartTransform.
   */
  private updateViewRange(boundStates: BoundStateResult): void {
    const xGrid = boundStates.xGrid;
    const xMin = xGrid[0] * QuantumConstants.M_TO_NM;
    const xMax = xGrid[xGrid.length - 1] * QuantumConstants.M_TO_NM;

    // Calculate energy range
    const energies = boundStates.energies.map((e) => e * QuantumConstants.JOULES_TO_EV);
    const maxEnergy = Math.max(...energies);
    const yMin = -this.model.wellDepthProperty.value * 0.2; // Small margin below zero
    const yMax = maxEnergy * 1.2; // 20% margin above highest energy

    this.xMinProperty.value = xMin;
    this.xMaxProperty.value = xMax;
    this.yMinProperty.value = yMin;
    this.yMaxProperty.value = yMax;

    // Update ChartTransform with new ranges
    this.chartTransform.setModelXRange(new Range(this.xMinProperty.value, this.xMaxProperty.value));
    this.chartTransform.setModelYRange(new Range(this.yMinProperty.value, this.yMaxProperty.value));
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

    if (potentialType === PotentialType.INFINITE_WELL) {
      // Draw square well
      const x1 = this.dataToViewX(0);
      const x2 = this.dataToViewX(wellWidth);
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
      // Draw finite square well
      const x1 = this.dataToViewX(0);
      const x2 = this.dataToViewX(wellWidth);
      const y0 = this.dataToViewY(0);
      const yDepth = this.dataToViewY(-wellDepth);

      shape.moveTo(this.chartMargins.left, y0);
      shape.lineTo(x1, y0);
      shape.lineTo(x1, yDepth);
      shape.lineTo(x2, yDepth);
      shape.lineTo(x2, y0);
      shape.lineTo(this.chartWidth - this.chartMargins.right, y0);
    } else if (potentialType === PotentialType.HARMONIC_OSCILLATOR) {
      // Draw parabola
      const centerX = wellWidth / 2;
      const numPoints = 100;
      let firstPoint = true;

      for (let i = 0; i < numPoints; i++) {
        const x = (xGrid[0] + (xGrid[xGrid.length - 1] - xGrid[0]) * i / (numPoints - 1)) * QuantumConstants.M_TO_NM;
        const dx = x - centerX;
        const k = (2 * wellDepth) / (wellWidth * wellWidth / 4); // Spring constant
        const V = 0.5 * k * dx * dx;
        const viewX = this.dataToViewX(x);
        const viewY = this.dataToViewY(V);

        if (firstPoint) {
          shape.moveTo(viewX, viewY);
          firstPoint = false;
        } else {
          shape.lineTo(viewX, viewY);
        }
      }
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
      this.removeChild(line);
    });
    this.energyLevelNodes.clear();

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
      this.addChild(line);

      // Add energy label
      const label = new Text(`E${index + 1} = ${energy.toFixed(3)} eV`, {
        font: "10px sans-serif",
        fill: QPPWColors.labelFillProperty,
        left: x2 + 5,
        centerY: y,
      });
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
   * Note: ChartTransform already handles the Y-axis inversion.
   */
  private dataToViewY(y: number): number {
    return this.chartMargins.top + (this.plotHeight - this.chartTransform.modelToViewY(y));
  }
}
