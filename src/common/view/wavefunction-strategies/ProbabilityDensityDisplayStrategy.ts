/**
 * ProbabilityDensityDisplayStrategy handles the probability density display mode.
 * This mode shows |ψ(x)|² as a filled area under the curve, along with
 * RMS position indicators (average position and uncertainty).
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
import { calculateRMSStatistics, createDoubleArrowShape } from "../RMSIndicatorUtils.js";
import stringManager from "../../../i18n/StringManager.js";

export class ProbabilityDensityDisplayStrategy extends WaveFunctionDisplayStrategy {
  /**
   * Calculate Y-axis range for probability density mode.
   * Range is from 0 (probability density is always non-negative) to maximum value.
   */
  calculateYRange(data: WavefunctionData): Range {
    // Find maximum probability density value
    const maxValue = Math.max(...data.probabilityDensity);

    // Range from 0 to max (probability density is always >= 0)
    let range = new Range(0, maxValue);

    // Add 10% padding at top
    range = this.addPadding(range.min, range.max);

    // Ensure minimum range
    range = this.ensureMinimumRange(range);

    // Validate range
    if (!this.isValidRange(range)) {
      console.warn(
        "[ProbabilityDensityDisplayStrategy] Invalid range detected:",
        range,
        "Using defaults",
      );
      return new Range(0, 1);
    }

    return range;
  }

  /**
   * Render probability density as a filled area under the curve.
   */
  render(data: WavefunctionData, context: RenderContext): void {
    const { xGrid, probabilityDensity } = data;
    const { probabilityDensityPath, dataToViewX, dataToViewY } = context;

    // Hide wavefunction component paths
    context.realPartPath.visible = false;
    context.imaginaryPartPath.visible = false;
    context.magnitudePath.visible = false;

    // Plot probability density with filled area
    const shape = this.createFilledCurveShape(
      xGrid,
      probabilityDensity,
      dataToViewX,
      dataToViewY,
    );
    probabilityDensityPath.shape = shape;
  }

  /**
   * Create a filled area shape under the probability density curve.
   */
  private createFilledCurveShape(
    xGrid: number[],
    probabilityDensity: number[],
    dataToViewX: (x: number) => number,
    dataToViewY: (y: number) => number,
  ): Shape {
    const shape = new Shape();

    // Build points array
    const points: { x: number; y: number }[] = [];
    for (let i = 0; i < xGrid.length; i++) {
      const x = dataToViewX(xGrid[i] * QuantumConstants.M_TO_NM);
      const y = dataToViewY(probabilityDensity[i]);
      points.push({ x, y });
    }

    if (points.length === 0) {
      return shape;
    }

    // Create filled area under the curve
    const y0 = dataToViewY(0); // baseline

    // Start at bottom-left
    shape.moveTo(points[0].x, y0);
    shape.lineTo(points[0].x, points[0].y);

    // Trace the curve
    for (let i = 0; i < points.length - 1; i++) {
      shape.lineTo(points[i].x, points[i].y);
    }

    // Close the shape back to baseline
    shape.lineTo(points[points.length - 1].x, points[points.length - 1].y);
    shape.lineTo(points[points.length - 1].x, y0);
    shape.close();

    return shape;
  }

  /**
   * Update RMS indicators for probability density mode.
   * Shows average position <x> and position uncertainty Δx.
   */
  updateRMSIndicators(
    data: WavefunctionData,
    context: RMSIndicatorContext,
  ): void {
    if (!context.shouldShow) {
      // Hide indicators
      context.avgPositionLabel.string = "";
      context.rmsPositionLabel.string = "";
      context.avgPositionIndicator.setLine(0, 0, 0, 0);
      context.rmsPositionIndicator.shape = null;
      return;
    }

    const { xGrid, probabilityDensity } = data;

    // Convert xGrid from meters to nanometers for RMS calculation
    const xGridNm = xGrid.map((x) => x * 1e9);

    // Calculate average position and RMS uncertainty
    const { avg, rms } = calculateRMSStatistics(xGridNm, probabilityDensity);

    // Update labels
    context.avgPositionLabel.string =
      stringManager.averagePositionLabelStringProperty.value.replace(
        "{{value}}",
        avg.toFixed(2),
      );
    context.rmsPositionLabel.string =
      stringManager.rmsPositionLabelStringProperty.value.replace(
        "{{value}}",
        rms.toFixed(2),
      );

    // Update average position indicator: vertical line at ⟨x⟩
    const avgX = context.dataToViewX(avg);
    const yTop = context.dataToViewY(context.yMaxValue);
    const yBottom = context.dataToViewY(0);
    context.avgPositionIndicator.setLine(avgX, yTop, avgX, yBottom);

    // Update RMS indicator: horizontal double arrow from (avg - rms) to (avg + rms)
    const leftX = avg - rms;
    const rightX = avg + rms;
    const x1 = context.dataToViewX(leftX);
    const x2 = context.dataToViewX(rightX);
    // Position the indicator at 80% of the visible range
    const indicatorY = context.dataToViewY(context.yMaxValue * 0.8);
    context.rmsPositionIndicator.shape = createDoubleArrowShape(x1, x2, indicatorY);
  }

  /**
   * Get Y-axis label for probability density mode.
   */
  getYAxisLabel(): string {
    return "Probability Density (nm⁻¹)";
  }

  /**
   * Format state label for probability density mode.
   * Example: ψ₁ → |ψ₁|²
   */
  getStateLabel(stateLabel: string): string {
    return stringManager.stateLabelProbabilityStringProperty.value.replace(
      "{{label}}",
      stateLabel,
    );
  }

  /**
   * Probability density mode shows zeros visualization.
   */
  shouldShowZeros(): boolean {
    return true;
  }

  /**
   * Get display mode identifier.
   */
  getDisplayMode(): string {
    return "probabilityDensity";
  }
}
