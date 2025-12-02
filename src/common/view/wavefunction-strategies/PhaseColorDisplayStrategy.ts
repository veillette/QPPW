/**
 * PhaseColorDisplayStrategy handles the phase color display mode.
 * This mode shows the wavefunction magnitude with colors representing the complex phase.
 * It delegates the actual rendering to the PhaseColorVisualization component.
 */

import { Range } from "scenerystack/dot";
import {
  WaveFunctionDisplayStrategy,
  type WavefunctionData,
  type RenderContext,
  type RMSIndicatorContext,
} from "./WaveFunctionDisplayStrategy.js";
import { PhaseColorVisualization } from "../chart-tools/PhaseColorVisualization.js";
import stringManager from "../../../i18n/StringManager.js";

export class PhaseColorDisplayStrategy extends WaveFunctionDisplayStrategy {
  private readonly phaseColorVisualization: PhaseColorVisualization;

  constructor(phaseColorVisualization: PhaseColorVisualization) {
    super();
    this.phaseColorVisualization = phaseColorVisualization;
  }

  /**
   * Calculate Y-axis range for phase color mode.
   * Range is from 0 to maximum magnitude (magnitude is always non-negative).
   */
  calculateYRange(data: WavefunctionData): Range {
    // Find maximum of probability density (which is |ψ|² in nm^-1)
    // For phase color, we want the range of |ψ| in nm^-1/2
    const maxValue = Math.max(...data.probabilityDensity);

    // Range from 0 to max (magnitude is always >= 0)
    let range = new Range(0, maxValue);

    // Add 10% padding at top
    range = this.addPadding(range.min, range.max);

    // Ensure minimum range
    range = this.ensureMinimumRange(range);

    // Validate range
    if (!this.isValidRange(range)) {
      console.warn(
        "[PhaseColorDisplayStrategy] Invalid range detected:",
        range,
        "Using defaults",
      );
      return new Range(0, 1);
    }

    return range;
  }

  /**
   * Render phase-colored wavefunction visualization.
   * Delegates to PhaseColorVisualization component.
   */
  render(data: WavefunctionData, context: RenderContext): void {
    const { xGrid, realPart, imagPart } = data;
    const {
      realPartPath,
      imaginaryPartPath,
      magnitudePath,
      probabilityDensityPath,
    } = context;

    // Hide other visualization paths
    probabilityDensityPath.shape = null;
    realPartPath.visible = false;
    imaginaryPartPath.visible = false;
    magnitudePath.visible = false;

    // Show and update phase color visualization
    this.phaseColorVisualization.show();
    this.phaseColorVisualization.plotSuperposition(xGrid, realPart, imagPart);
  }

  /**
   * Render phase-colored wavefunction for a single eigenstate.
   * This variant handles time evolution phase for single states.
   */
  renderSingleState(
    data: WavefunctionData,
    context: RenderContext,
    globalPhase: number,
  ): void {
    const { xGrid } = data;
    const {
      realPartPath,
      imaginaryPartPath,
      magnitudePath,
      probabilityDensityPath,
    } = context;

    // Hide other visualization paths
    probabilityDensityPath.shape = null;
    realPartPath.visible = false;
    imaginaryPartPath.visible = false;
    magnitudePath.visible = false;

    // Show and update phase color visualization with global phase
    this.phaseColorVisualization.show();

    // For single eigenstate, the wavefunction is real, so we apply the global phase
    // to get the time-evolved complex wavefunction
    const realTimeDep: number[] = [];
    const imagTimeDep: number[] = [];

    for (let i = 0; i < data.realPart.length; i++) {
      // Since eigenstate wavefunction is real, apply time evolution phase
      // ψ(x,t) = ψ(x) * exp(-iEt/ℏ) = ψ(x) * (cos(Et/ℏ) - i*sin(Et/ℏ))
      const psi = data.realPart[i];
      realTimeDep.push(psi * Math.cos(globalPhase));
      imagTimeDep.push(-psi * Math.sin(globalPhase));
    }

    this.phaseColorVisualization.plotSuperposition(
      xGrid,
      realTimeDep,
      imagTimeDep,
    );
  }

  /**
   * Hide phase color visualization.
   */
  hide(): void {
    this.phaseColorVisualization.hide();
  }

  /**
   * Update RMS indicators for phase color mode.
   * This mode does not show RMS indicators, so they are hidden.
   */
  updateRMSIndicators(
    _data: WavefunctionData,
    context: RMSIndicatorContext,
  ): void {
    // Hide all RMS indicators in phase color mode
    context.avgPositionLabel.string = "";
    context.rmsPositionLabel.string = "";
    context.avgPositionIndicator.setLine(0, 0, 0, 0);
    context.rmsPositionIndicator.shape = null;
  }

  /**
   * Get Y-axis label for phase color mode.
   */
  getYAxisLabel(): string {
    return "Wave Function Magnitude (nm⁻¹ᐟ²)";
  }

  /**
   * Format state label for phase color mode.
   * Example: ψ₁ → ψ₁ (wrapped in wavefunction label template)
   */
  getStateLabel(stateLabel: string): string {
    return stringManager.stateLabelWavefunctionStringProperty.value.replace(
      "{{label}}",
      stateLabel,
    );
  }

  /**
   * Phase color mode does NOT show zeros visualization.
   * Zeros are not meaningful in phase color representation.
   */
  shouldShowZeros(): boolean {
    return false;
  }

  /**
   * Get display mode identifier.
   */
  getDisplayMode(): string {
    return "phaseColor";
  }
}
