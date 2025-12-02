/**
 * DerivativeTool provides interactive visualization of the first derivative (slope)
 * of the wavefunction at a specific point. Shows a tangent line.
 */

import { Line, Path, Text, Circle } from "scenerystack/scenery";
import { Shape } from "scenerystack/kite";
import { NumberProperty } from "scenerystack/axon";
import { PhetFont } from "scenerystack/scenery-phet";
import type { ScreenModel } from "../../model/ScreenModels.js";
import QPPWColors from "../../../QPPWColors.js";
import stringManager from "../../../i18n/StringManager.js";
import { Utterance } from "scenerystack/utterance-queue";
import {
  BaseChartTool,
  utteranceQueue,
  type BaseChartToolOptions,
} from "./BaseChartTool.js";

export type DerivativeToolOptions = BaseChartToolOptions;

export class DerivativeTool extends BaseChartTool {
  private readonly markerXProperty: NumberProperty; // X position in nm
  private marker!: Line;
  private positionCircle!: Circle; // Circle tracking wavefunction position
  private tangentLine!: Path;
  private label!: Text;

  constructor(
    model: ScreenModel,
    getEffectiveDisplayMode: () => string,
    options: DerivativeToolOptions,
  ) {
    super(
      model,
      getEffectiveDisplayMode,
      options,
      "Derivative Visualization",
      "Showing first derivative dψ/dx (slope) of the wavefunction. " +
        "The tangent line shows the rate of change at the selected position.",
    );

    this.markerXProperty = new NumberProperty(1);

    // Update when marker position changes
    this.markerXProperty.link(() => {
      if (this.showProperty.value) {
        this.update(getEffectiveDisplayMode());
      }
    });
  }

  /**
   * Setup visual elements for the derivative tool
   */
  protected setupVisualElements(): void {
    // Create marker line
    this.marker = new Line(0, 0, 0, 0, {
      stroke: QPPWColors.derivativeToolStrokeProperty,
      lineWidth: 2,
      lineDash: [4, 3],
      cursor: "ew-resize",
    });
    this.container.addChild(this.marker);

    // Create position tracking circle (shows position on wavefunction)
    this.positionCircle = new Circle(5, {
      fill: QPPWColors.derivativeToolFillDarkProperty,
      stroke: QPPWColors.backgroundColorProperty,
      lineWidth: 2,
    });
    this.container.addChild(this.positionCircle);

    // Create tangent line path
    this.tangentLine = new Path(null, {
      stroke: QPPWColors.derivativeToolFillLightProperty,
      lineWidth: 3,
      fill: null,
    });
    this.container.addChild(this.tangentLine);

    // Create derivative label
    this.label = new Text("", {
      font: new PhetFont({ size: 14, weight: "bold" }),
      fill: QPPWColors.derivativeToolFillDarkProperty,
    });
    this.container.addChild(this.label);

    // Setup drag listener for marker using base class method
    this.setupMarkerDragListener({
      markerNode: this.marker,
      positionProperty: this.markerXProperty,
      accessibleName: "Derivative Measurement Position",
      labelContent: "Position for derivative and slope display",
      helpText:
        "Use Left/Right arrow keys to move marker. " +
        "Shift+Arrow for fine control. " +
        "Page Up/Down for large steps. " +
        "Shows first derivative (slope) at selected position.",
      onDragEnd: (position) => {
        const derivativeData = this.calculateFirstDerivative(
          position,
          "waveFunction",
        );
        if (derivativeData !== null) {
          const message =
            `Marker at ${position.toFixed(2)} nanometers. ` +
            `First derivative: ${derivativeData.firstDerivative.toFixed(3)}.`;
          utteranceQueue.addToBack(new Utterance({ alert: message }));
        }
      },
    });
  }

  /**
   * Updates the derivative tool visualization.
   */
  public update(displayMode: string): void {
    if (!this.showProperty.value) {
      return;
    }

    const x = this.markerXProperty.value;
    const viewX = this.dataToViewX(x);
    const yTop = this.getChartTop();
    const yBottom = this.getChartBottom();

    // Update marker line
    this.marker.setLine(viewX, yTop, viewX, yBottom);

    // Calculate and display first derivative
    const derivativeData = this.calculateFirstDerivative(x, displayMode);

    if (derivativeData !== null) {
      // Create tangent line shape with first derivative
      const tangentLineShape = this.createTangentLine(
        x,
        derivativeData.firstDerivative,
        derivativeData.wavefunctionValue,
      );
      this.tangentLine.shape = tangentLineShape;

      // Update position tracking circle
      this.positionCircle.centerX = viewX;
      this.positionCircle.centerY = this.dataToViewY(
        derivativeData.wavefunctionValue,
      );

      // Update label with proper units (nm^-3/2)
      this.label.string =
        stringManager.firstDerivativeLabelStringProperty.value.replace(
          "{{value}}",
          derivativeData.firstDerivative.toFixed(3),
        );
      this.label.centerX = viewX;
      this.label.top = yTop + 5; // Position just inside the chart area
    } else {
      this.tangentLine.shape = null;
      this.label.string = stringManager.notAvailableStringProperty.value;
    }
  }

  /**
   * Calculates the first derivative of the wavefunction at a given x position.
   */
  private calculateFirstDerivative(
    xData: number,
    displayMode: string,
  ): {
    firstDerivative: number;
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
      wavefunctionValue: result.value,
    };
  }

  /**
   * Creates a tangent line shape that shows the slope at a point.
   * Uses the linear approximation: y = y₀ + y'₀(x - x₀)
   */
  private createTangentLine(
    xData: number,
    firstDerivative: number,
    wavefunctionValue: number,
  ): Shape {
    const shape = new Shape();

    // Create a tangent line centered at (xData, wavefunctionValue)
    const centerX = xData;
    const centerY = wavefunctionValue;

    // Width of tangent line in nm (±0.5 nm from center)
    const tangentHalfWidth = 0.5;

    // Calculate endpoints of the tangent line
    const x1 = centerX - tangentHalfWidth;
    const y1 = centerY + firstDerivative * -tangentHalfWidth;
    const x2 = centerX + tangentHalfWidth;
    const y2 = centerY + firstDerivative * tangentHalfWidth;

    // Convert to view coordinates
    const viewX1 = this.dataToViewX(x1);
    const viewY1 = this.dataToViewY(y1);
    const viewX2 = this.dataToViewX(x2);
    const viewY2 = this.dataToViewY(y2);

    // Draw the tangent line
    shape.moveTo(viewX1, viewY1);
    shape.lineTo(viewX2, viewY2);

    return shape;
  }
}
