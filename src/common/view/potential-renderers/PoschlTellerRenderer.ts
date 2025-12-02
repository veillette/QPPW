/**
 * Renderer for Pöschl-Teller potential.
 * V(x) = -V_0 / cosh²(x/a)
 * This is a symmetric potential well that goes to 0 as x→±∞
 */

import { Shape } from "scenerystack/kite";
import QuantumConstants from "../../model/QuantumConstants.js";
import { PotentialRenderer, RenderContext } from "./PotentialRenderer.js";

export class PoschlTellerRenderer extends PotentialRenderer {
  public render(context: RenderContext): Shape {
    const { xGrid, xCenter, wellWidth, wellDepth, dataToViewX, dataToViewY } = context;
    const shape = new Shape();

    const centerX = xCenter;
    const numPoints = 200;
    let firstPoint = true;

    this.iterateXGrid(xGrid, numPoints, (xMeters) => {
      const x = xMeters * QuantumConstants.M_TO_NM;
      const dx = x - centerX;
      const coshVal = Math.cosh(dx / wellWidth);
      const V = -wellDepth / (coshVal * coshVal);

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
