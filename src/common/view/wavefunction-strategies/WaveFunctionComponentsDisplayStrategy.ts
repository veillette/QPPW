/**
 * WaveFunctionComponentsDisplayStrategy handles the wavefunction components display mode.
 * This mode shows the real part Re(ψ), imaginary part Im(ψ), and magnitude |ψ|
 * as separate line plots with optional visibility toggles.
 */

import { Range } from "scenerystack/dot";
import { Shape } from "scenerystack/kite";
import {
  WaveFunctionDisplayStrategy,
  type WavefunctionData,
  type RenderContext,
  type RMSIndicatorContext,
} from "./WaveFunctionDisplayStrategy.js";
import QuantumConstants from "../../model/QuantumConstants.js";

export class WaveFunctionComponentsDisplayStrategy extends WaveFunctionDisplayStrategy {
  /**
   * Calculate Y-axis range for wavefunction components mode.
   * Range is symmetric around zero (wavefunction can be negative).
   */
  calculateYRange(data: WavefunctionData): Range {
    // Find maximum absolute value across all components
    const maxReal = Math.max(...data.realPart.map(Math.abs));
    const maxImag = Math.max(...data.imagPart.map(Math.abs));
    const maxMagnitude = data.maxMagnitude;
    const maxAbs = Math.max(maxReal, maxImag, maxMagnitude);

    // Symmetric range around zero
    let range = new Range(-maxAbs, maxAbs);

    // Add 10% padding
    range = this.addPadding(range.min, range.max);

    // Ensure minimum range
    range = this.ensureMinimumRange(range);

    // Validate range
    if (!this.isValidRange(range)) {
      console.warn(
        "[WaveFunctionComponentsDisplayStrategy] Invalid range detected:",
        range,
        "Using defaults",
      );
      return this.getDefaultRange();
    }

    return range;
  }

  /**
   * Render wavefunction components (real, imaginary, magnitude) as separate lines.
   */
  render(data: WavefunctionData, context: RenderContext): void {
    const { xGrid, realPart, imagPart } = data;
    const {
      realPartPath,
      imaginaryPartPath,
      magnitudePath,
      probabilityDensityPath,
      dataToViewX,
      dataToViewY,
      showRealPart,
      showImaginaryPart,
      showMagnitude,
    } = context;

    // Hide probability density
    probabilityDensityPath.shape = null;

    // Build points for each component
    const realPoints: { x: number; y: number }[] = [];
    const imagPoints: { x: number; y: number }[] = [];
    const magnitudePoints: { x: number; y: number }[] = [];

    for (let i = 0; i < xGrid.length; i++) {
      const x = dataToViewX(xGrid[i] * QuantumConstants.M_TO_NM);

      // Calculate magnitude from real and imaginary parts
      const magnitude = Math.sqrt(
        realPart[i] * realPart[i] + imagPart[i] * imagPart[i],
      );

      realPoints.push({ x, y: dataToViewY(realPart[i]) });
      imagPoints.push({ x, y: dataToViewY(imagPart[i]) });
      magnitudePoints.push({ x, y: dataToViewY(magnitude) });
    }

    // Plot real part
    realPartPath.shape = this.createLineShape(realPoints);
    realPartPath.visible = showRealPart;

    // Plot imaginary part
    imaginaryPartPath.shape = this.createLineShape(imagPoints);
    imaginaryPartPath.visible = showImaginaryPart;

    // Plot magnitude
    magnitudePath.shape = this.createLineShape(magnitudePoints);
    magnitudePath.visible = showMagnitude;
  }

  /**
   * Create a line shape from an array of points.
   */
  private createLineShape(points: { x: number; y: number }[]): Shape {
    const shape = new Shape();

    if (points.length === 0) {
      return shape;
    }

    // Move to first point
    shape.moveTo(points[0].x, points[0].y);

    // Draw lines to remaining points
    for (let i = 1; i < points.length; i++) {
      shape.lineTo(points[i].x, points[i].y);
    }

    return shape;
  }

  /**
   * Update RMS indicators for wavefunction components mode.
   * This mode does not show RMS indicators, so they are hidden.
   */
  updateRMSIndicators(
    _data: WavefunctionData,
    context: RMSIndicatorContext,
  ): void {
    // Hide all RMS indicators in wavefunction components mode
    context.avgPositionLabel.string = "";
    context.rmsPositionLabel.string = "";
    context.avgPositionIndicator.setLine(0, 0, 0, 0);
    context.rmsPositionIndicator.shape = null;
  }

  /**
   * Get Y-axis label for wavefunction components mode.
   */
  getYAxisLabel(): string {
    return "Wave Function (nm⁻¹ᐟ²)";
  }

  /**
   * Format state label for wavefunction components mode.
   * Returns the label unchanged (e.g., ψ₁ stays as ψ₁).
   */
  getStateLabel(stateLabel: string): string {
    return stateLabel;
  }

  /**
   * Wavefunction components mode shows zeros visualization.
   */
  shouldShowZeros(): boolean {
    return true;
  }

  /**
   * Get display mode identifier.
   */
  getDisplayMode(): string {
    return "waveFunction";
  }
}
