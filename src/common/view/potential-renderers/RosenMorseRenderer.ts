/**
 * Renderer for Rosen-Morse potential.
 * V(x) = -V_0 / cosh²(x/a) + V_1 * tanh(x/a)
 * This adds an asymmetric term to the Pöschl-Teller potential
 */

import { Shape } from "scenerystack/kite";
import QuantumConstants from "../../model/QuantumConstants.js";
import { hasBarrierHeight } from "../../model/ModelTypeGuards.js";
import { PotentialRenderer, RenderContext } from "./PotentialRenderer.js";

export class RosenMorseRenderer extends PotentialRenderer {
  public render(context: RenderContext): Shape {
    const { xGrid, xCenter, wellWidth, wellDepth, dataToViewX, dataToViewY, model } = context;
    const shape = new Shape();

    const centerX = xCenter;
    const barrierHeight = hasBarrierHeight(model)
      ? model.barrierHeightProperty.value
      : 0;
    const numPoints = 200;
    let firstPoint = true;

    this.iterateXGrid(xGrid, numPoints, (xMeters) => {
      const x = xMeters * QuantumConstants.M_TO_NM;
      const dx = x - centerX;
      const coshVal = Math.cosh(dx / wellWidth);
      const tanhVal = Math.tanh(dx / wellWidth);
      const V = -wellDepth / (coshVal * coshVal) + barrierHeight * tanhVal;

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
