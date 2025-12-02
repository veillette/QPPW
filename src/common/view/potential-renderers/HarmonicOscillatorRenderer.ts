/**
 * Renderer for harmonic oscillator potential.
 * V(x) = (1/2) * k * (x - xCenter)²
 * where k = 2 * wellDepth / (wellWidth²/4)
 */

import { Shape } from "scenerystack/kite";
import QuantumConstants from "../../model/QuantumConstants.js";
import { PotentialRenderer, RenderContext } from "./PotentialRenderer.js";

export class HarmonicOscillatorRenderer extends PotentialRenderer {
  public render(context: RenderContext): Shape {
    const { xGrid, xCenter, wellWidth, wellDepth, dataToViewX, dataToViewY } = context;
    const shape = new Shape();

    // Draw parabola centered at xCenter
    const centerX = xCenter;
    const numPoints = 100;
    const k = (2 * wellDepth) / ((wellWidth * wellWidth) / 4); // Spring constant
    let firstPoint = true;

    this.iterateXGrid(xGrid, numPoints, (xMeters) => {
      const x = xMeters * QuantumConstants.M_TO_NM;
      const dx = x - centerX;
      const V = 0.5 * k * dx * dx;
      const viewX = dataToViewX(x);
      const viewY = dataToViewY(V);

      if (firstPoint) {
        shape.moveTo(viewX, viewY);
        firstPoint = false;
      } else {
        shape.lineTo(viewX, viewY);
      }
    });

    return shape;
  }
}
