/**
 * PhaseColorVisualization displays the wavefunction with colors representing the phase
 * of the complex wavefunction. The height represents magnitude, and the hue represents phase.
 */

import { Rectangle } from "scenerystack/scenery";
import QuantumConstants from "../../model/QuantumConstants.js";
import QPPWColors from "../../../QPPWColors.js";
import {
  BaseVisualization,
  type BaseVisualizationOptions,
} from "./BaseVisualization.js";

export type PhaseColorVisualizationOptions = BaseVisualizationOptions;

export class PhaseColorVisualization extends BaseVisualization {
  private phaseColorStrips: Rectangle[] = [];

  constructor(options: PhaseColorVisualizationOptions) {
    super(options);
  }

  /**
   * Plot phase-colored wavefunction for a single eigenstate
   */
  public plotWavefunction(
    xGrid: number[],
    wavefunction: number[],
    globalPhase: number,
  ): void {
    const y0 = this.dataToViewY(0);
    const numStrips = xGrid.length - 1;

    // Ensure we have enough rectangles in the pool
    while (this.phaseColorStrips.length < numStrips) {
      const rect = new Rectangle(0, 0, 1, 1, {
        fill: QPPWColors.phaseIndicatorProperty,
        stroke: null,
      });
      this.phaseColorStrips.push(rect);
      this.container.addChild(rect);
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
      const normalizedPhase = (localPhase + Math.PI) / (2 * Math.PI);
      const hue = Math.round(normalizedPhase * 360);

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
   * Plot phase-colored superposition from real and imaginary parts
   */
  public plotSuperposition(
    xGrid: number[],
    realPart: number[],
    imagPart: number[],
  ): void {
    const y0 = this.dataToViewY(0);
    const numStrips = xGrid.length - 1;

    // Ensure we have enough rectangles in the pool
    while (this.phaseColorStrips.length < numStrips) {
      const rect = new Rectangle(0, 0, 1, 1, {
        fill: QPPWColors.phaseIndicatorProperty,
        stroke: null,
      });
      this.phaseColorStrips.push(rect);
      this.container.addChild(rect);
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
      const normalizedPhase = (localPhase + Math.PI) / (2 * Math.PI);
      const hue = Math.round(normalizedPhase * 360);

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
}
