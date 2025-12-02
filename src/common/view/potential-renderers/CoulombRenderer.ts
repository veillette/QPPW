/**
 * Renderer for Coulomb potential (both 1D and 3D).
 * V(x) = -k/|x| where k is determined by wellDepth
 * V→0 as x→∞, V→-wellDepth at distance wellWidth/2
 * At distance wellWidth/2, V = -wellDepth
 */

import { Shape } from "scenerystack/kite";
import QuantumConstants from "../../model/QuantumConstants.js";
import { PotentialRenderer, RenderContext } from "./PotentialRenderer.js";

export class CoulombRenderer extends PotentialRenderer {
  public render(context: RenderContext): Shape {
    const { xGrid, xCenter, wellWidth, wellDepth, dataToViewX, dataToViewY, model } = context;
    const shape = new Shape();

    const centerX = xCenter;
    const numPoints = 200;
    let firstPoint = true;

    // Scale factor: at distance wellWidth/2, V = -wellDepth
    const k = wellDepth * (wellWidth / 2);

    this.iterateXGrid(xGrid, numPoints, (xMeters) => {
      const x = xMeters * QuantumConstants.M_TO_NM;
      const dx = x - centerX;

      // Avoid singularity at x=0
      const distance = Math.max(Math.abs(dx), 0.01);
      const V = -k / distance;

      // Clamp V to reasonable range for display
      // Access yMinProperty from the model's view (not ideal, but needed for clamping)
      const VClamped = Math.max(V, -20); // Reasonable default clamp value

      const viewX = dataToViewX(x);
      const viewY = dataToViewY(VClamped);

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
