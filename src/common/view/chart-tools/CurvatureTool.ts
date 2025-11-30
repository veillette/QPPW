/**
 * CurvatureTool provides interactive visualization of the second derivative (curvature)
 * of the wavefunction at a specific point. Shows a parabola based on Taylor expansion.
 */

import {
  Node,
  Line,
  Path,
  Text,
  Circle,
  DragListener,
  KeyboardDragListener,
  DerivedProperty,
} from "scenerystack/scenery";
import { Shape } from "scenerystack/kite";
import { NumberProperty, BooleanProperty } from "scenerystack/axon";
import { PhetFont } from "scenerystack/scenery-phet";
import type { ScreenModel } from "../../model/ScreenModels.js";
import QPPWColors from "../../../QPPWColors.js";
import stringManager from "../../../i18n/StringManager.js";
import { utteranceQueue } from "scenerystack/utterance-queue";

export type CurvatureToolOptions = {
  chartMargins: { left: number; right: number; top: number; bottom: number };
  plotHeight: number;
  xMinProperty: NumberProperty;
  xMaxProperty: NumberProperty;
  dataToViewX: (x: number) => number;
  dataToViewY: (y: number) => number;
  viewToDataX: (x: number) => number;
  parentNode: Node; // Reference to parent chart node for coordinate conversions
};

export class CurvatureTool extends Node {
  private readonly model: ScreenModel;
  private readonly options: CurvatureToolOptions;

  public readonly showProperty: BooleanProperty;
  private readonly markerXProperty: NumberProperty; // X position in nm

  private readonly container: Node;
  private readonly marker: Line;
  private readonly markerHandle: Circle;
  private readonly positionCircle: Circle; // Circle tracking wavefunction position
  private readonly parabola: Path;
  private readonly label: Text;

  private readonly showPropertyInternal: BooleanProperty;

  constructor(
    model: ScreenModel,
    getEffectiveDisplayMode: () => string,
    options: CurvatureToolOptions,
  ) {
    // Initialize properties first so they can be used in super()
    const showPropertyInternal = new BooleanProperty(false);

    super({
      // pdom - container for the curvature tool
      tagName: "div",
      labelTagName: "h3",
      labelContent: "Curvature Visualization",
      descriptionTagName: "p",
      descriptionContent: new DerivedProperty(
        [showPropertyInternal],
        (enabled) =>
          enabled
            ? "Showing second derivative d²ψ/dx². Curvature is proportional to (V(x) - E)ψ(x) " +
              "according to the Schrödinger equation. Positive curvature where V > E, negative where V < E."
            : "Curvature visualization disabled.",
      ),
    });

    this.model = model;
    this.options = options;

    // Store properties
    this.showPropertyInternal = showPropertyInternal;
    this.showProperty = showPropertyInternal;
    this.markerXProperty = new NumberProperty(0);

    // Create container
    this.container = new Node({
      visible: false,
    });
    this.addChild(this.container);

    // Create marker line
    this.marker = new Line(0, 0, 0, 0, {
      stroke: QPPWColors.curvatureToolStrokeProperty,
      lineWidth: 2,
      lineDash: [4, 3],
    });
    this.container.addChild(this.marker);

    // Create marker handle (draggable circle)
    this.markerHandle = new Circle(8, {
      fill: QPPWColors.curvatureToolFillLightProperty,
      stroke: QPPWColors.backgroundColorProperty,
      lineWidth: 2,
      cursor: "ew-resize",

      // pdom - keyboard accessible marker
      tagName: "div",
      ariaRole: "slider",
      focusable: true,
      accessibleName: "Curvature Measurement Position",
      labelContent: "Position for curvature and second derivative display",
      ariaValueText: new DerivedProperty(
        [this.markerXProperty],
        (position) => `Position: ${position.toFixed(2)} nanometers`,
      ),
      ariaValueMin: -5,
      ariaValueMax: 5,
      ariaValueNow: this.markerXProperty,
      helpText:
        "Use Left/Right arrow keys to move marker. " +
        "Shift+Arrow for fine control. " +
        "Page Up/Down for large steps. " +
        "Shows second derivative (curvature) at selected position.",
    });
    this.container.addChild(this.markerHandle);

    // Create position tracking circle (shows position on wavefunction)
    this.positionCircle = new Circle(5, {
      fill: QPPWColors.curvatureToolFillDarkProperty,
      stroke: QPPWColors.backgroundColorProperty,
      lineWidth: 2,
    });
    this.container.addChild(this.positionCircle);

    // Create parabola path
    this.parabola = new Path(null, {
      stroke: QPPWColors.curvatureToolFillLightProperty,
      lineWidth: 3,
      fill: null,
    });
    this.container.addChild(this.parabola);

    // Create curvature label
    this.label = new Text("", {
      font: new PhetFont({ size: 14, weight: "bold" }),
      fill: QPPWColors.curvatureToolFillDarkProperty,
    });
    this.container.addChild(this.label);

    // Setup drag listener for marker
    this.setupDragListener();

    // Store display mode getter
    const getDisplayMode = getEffectiveDisplayMode;

    // Link visibility to property
    this.showProperty.link((show: boolean) => {
      this.container.visible = show;
      if (show) {
        this.update(getDisplayMode());
      }
    });

    // Update when marker position changes
    this.markerXProperty.link(() => {
      if (this.showProperty.value) {
        this.update(getDisplayMode());
      }
    });
  }

  /**
   * Setup drag listeners for marker handle (both mouse and keyboard)
   */
  private setupDragListener(): void {
    const {
      parentNode,
      chartMargins,
      viewToDataX,
      xMinProperty,
      xMaxProperty,
    } = this.options;

    // Mouse drag listener
    const dragListener = new DragListener({
      drag: (event) => {
        const parentPoint = parentNode.globalToParentPoint(event.pointer.point);
        const localX = parentPoint.x - chartMargins.left;
        const dataX = viewToDataX(localX);

        // Clamp to chart bounds
        const clampedX = Math.max(
          xMinProperty.value,
          Math.min(xMaxProperty.value, dataX),
        );
        this.markerXProperty.value = clampedX;
      },
    });

    this.markerHandle.addInputListener(dragListener);

    // Keyboard drag listener
    const keyboardDragListener = new KeyboardDragListener({
      drag: (vectorDelta) => {
        let newX = this.markerXProperty.value + vectorDelta.x * 0.1;

        // Clamp to chart bounds
        newX = Math.max(xMinProperty.value, Math.min(xMaxProperty.value, newX));

        this.markerXProperty.value = newX;
      },
      dragDelta: 1, // Regular arrow key step
      shiftDragDelta: 0.1, // Fine control with shift
      end: () => {
        // Announce position and curvature value on drag end
        const position = this.markerXProperty.value;
        const derivatives = this.calculateDerivatives(position, "waveFunction");

        if (derivatives !== null) {
          const message =
            `Marker at ${position.toFixed(2)} nanometers. ` +
            `Second derivative: ${derivatives.secondDerivative.toExponential(2)}.`;
          utteranceQueue.addToBack(message);
        }
      },
    });

    this.markerHandle.addInputListener(keyboardDragListener);
  }

  /**
   * Updates the curvature tool visualization.
   */
  public update(displayMode: string): void {
    if (!this.showProperty.value) {
      return;
    }

    const { chartMargins, plotHeight, dataToViewX } = this.options;
    const x = this.markerXProperty.value;
    const viewX = dataToViewX(x);
    const yTop = chartMargins.top;
    const yBottom = chartMargins.top + plotHeight;

    // Update marker line
    this.marker.setLine(viewX, yTop, viewX, yBottom);

    // Update marker handle
    this.markerHandle.centerX = viewX;
    this.markerHandle.centerY = yTop + 20;

    // Calculate and display derivatives
    const derivatives = this.calculateDerivatives(x, displayMode);

    if (derivatives !== null) {
      // Create parabola shape with first and second derivatives
      const parabolaShape = this.createCurvatureParabola(
        x,
        derivatives.firstDerivative,
        derivatives.secondDerivative,
        derivatives.wavefunctionValue,
      );
      this.parabola.shape = parabolaShape;

      // Update position tracking circle
      const { dataToViewY } = this.options;
      this.positionCircle.centerX = viewX;
      this.positionCircle.centerY = dataToViewY(derivatives.wavefunctionValue);

      // Update label with proper units (nm^-5/2)
      this.label.string =
        stringManager.secondDerivativeLabelStringProperty.value.replace(
          "{{value}}",
          derivatives.secondDerivative.toExponential(2),
        );
      this.label.centerX = viewX;
      this.label.top = yTop + 5; // Position just inside the chart area
    } else {
      this.parabola.shape = null;
      this.label.string = stringManager.notAvailableStringProperty.value;
    }
  }

  /**
   * Calculates the first and second derivatives of the wavefunction at a given x position.
   */
  private calculateDerivatives(
    xData: number,
    displayMode: string,
  ): {
    firstDerivative: number;
    secondDerivative: number;
    wavefunctionValue: number;
  } | null {
    // Only calculate for wavefunction mode
    if (displayMode !== "waveFunction") {
      return null;
    }

    const selectedIndex = this.model.selectedEnergyLevelIndexProperty.value;

    // Use model's calculation method
    const result = this.model.getWavefunctionAtPosition(selectedIndex, xData);
    if (!result) {
      return null;
    }

    return {
      firstDerivative: result.firstDerivative,
      secondDerivative: result.secondDerivative,
      wavefunctionValue: result.value,
    };
  }

  /**
   * Creates a parabola shape that matches the curvature at a point.
   * Uses Taylor expansion: y = y₀ + y'₀(x - x₀) + (1/2)y''₀(x - x₀)²
   */
  private createCurvatureParabola(
    xData: number,
    firstDerivative: number,
    secondDerivative: number,
    wavefunctionValue: number,
  ): Shape {
    const shape = new Shape();
    const { dataToViewX, dataToViewY } = this.options;

    // Create a small parabola centered at (xData, wavefunctionValue)
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

      const viewX = dataToViewX(x);
      const viewY = dataToViewY(y);

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
