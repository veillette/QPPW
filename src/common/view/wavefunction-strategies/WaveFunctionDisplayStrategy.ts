/**
 * WaveFunctionDisplayStrategy defines the interface for different display modes
 * in the WaveFunctionChartNode. This follows the Strategy pattern, similar to
 * how PotentialRenderer works for EnergyChartNode.
 *
 * Each display mode (probability density, wavefunction components, phase color)
 * has its own concrete implementation of this interface.
 */

import { Range } from "scenerystack/dot";
import type { BoundStateResult } from "../../model/PotentialFunction.js";
import type { ScreenModel } from "../../model/ScreenModels.js";
import { Path, Text, Line } from "scenerystack/scenery";
import { ChartTransform } from "scenerystack/bamboo";
import type { NumberProperty } from "scenerystack/axon";

/**
 * Data structure containing wavefunction information in nm units.
 */
export interface WavefunctionData {
  xGrid: number[]; // Position grid in meters
  realPart: number[]; // Real part in nm^-1/2
  imagPart: number[]; // Imaginary part in nm^-1/2
  magnitude: number[]; // Magnitude in nm^-1/2
  probabilityDensity: number[]; // |ψ|² in nm^-1
  maxMagnitude: number; // Maximum magnitude value
}

/**
 * Context object containing chart rendering resources and coordinate transforms.
 */
export interface RenderContext {
  // Visual elements (paths for rendering)
  realPartPath: Path;
  imaginaryPartPath: Path;
  magnitudePath: Path;
  probabilityDensityPath: Path;

  // Coordinate transformation functions
  dataToViewX: (x: number) => number;
  dataToViewY: (y: number) => number;

  // Chart dimensions and properties
  chartMargins: { left: number; right: number; top: number; bottom: number };
  yMinProperty: NumberProperty;
  yMaxProperty: NumberProperty;

  // Visibility control properties from view state
  showRealPart: boolean;
  showImaginaryPart: boolean;
  showMagnitude: boolean;
}

/**
 * Context for RMS indicator display.
 */
export interface RMSIndicatorContext {
  avgPositionIndicator: Line;
  rmsPositionIndicator: Path;
  avgPositionLabel: Text;
  rmsPositionLabel: Text;
  shouldShow: boolean;
  dataToViewX: (x: number) => number;
  dataToViewY: (y: number) => number;
  yMaxValue: number;
}

/**
 * Abstract base class for wavefunction display strategies.
 * Each display mode extends this class and implements mode-specific rendering logic.
 */
export abstract class WaveFunctionDisplayStrategy {
  /**
   * Calculate the appropriate Y-axis range for this display mode.
   * Different modes have different range requirements:
   * - Probability density: 0 to max (always positive)
   * - Wavefunction components: symmetric around 0 (can be negative)
   * - Phase color: 0 to max magnitude (always positive)
   *
   * @param data - Wavefunction data in nm units
   * @returns Range object with min and max values
   */
  abstract calculateYRange(data: WavefunctionData): Range;

  /**
   * Render the wavefunction visualization for this display mode.
   *
   * @param data - Wavefunction data in nm units
   * @param context - Rendering context with paths and coordinate transforms
   */
  abstract render(data: WavefunctionData, context: RenderContext): void;

  /**
   * Update RMS indicators (average position and uncertainty).
   * Only probability density mode shows these indicators.
   *
   * @param data - Wavefunction data in nm units
   * @param context - RMS indicator context
   */
  abstract updateRMSIndicators(
    data: WavefunctionData,
    context: RMSIndicatorContext,
  ): void;

  /**
   * Get the Y-axis label text for this display mode.
   *
   * @returns Axis label string (e.g., "Probability Density (nm⁻¹)")
   */
  abstract getYAxisLabel(): string;

  /**
   * Get the state label prefix for this display mode.
   * Used to format labels like "|ψ₁|²" for probability density or "ψ₁" for wavefunction.
   *
   * @param stateLabel - Base state label (e.g., "ψ₁")
   * @returns Formatted label for this display mode
   */
  abstract getStateLabel(stateLabel: string): string;

  /**
   * Check if zeros visualization should be shown in this mode.
   * Only wavefunction and probability density modes show zeros.
   *
   * @returns true if zeros should be displayed
   */
  abstract shouldShowZeros(): boolean;

  /**
   * Get the display mode name for this strategy.
   *
   * @returns Display mode identifier ("probabilityDensity", "waveFunction", or "phaseColor")
   */
  abstract getDisplayMode(): string;

  /**
   * Helper method: Add 10% padding to a range.
   */
  protected addPadding(min: number, max: number, paddingFraction = 0.1): Range {
    const range = max - min;
    const padding = range * paddingFraction;
    return new Range(min - padding, max + padding);
  }

  /**
   * Helper method: Ensure range is not too small (minimum 0.01).
   */
  protected ensureMinimumRange(range: Range): Range {
    if (range.max - range.min < 0.01) {
      return new Range(-0.01, 0.01);
    }
    return range;
  }

  /**
   * Helper method: Check if range values are valid (not too large).
   */
  protected isValidRange(range: Range): boolean {
    return Math.abs(range.min) < 1000 && Math.abs(range.max) < 1000;
  }

  /**
   * Helper method: Get default range if calculated range is invalid.
   */
  protected getDefaultRange(): Range {
    return new Range(-1, 1);
  }
}
