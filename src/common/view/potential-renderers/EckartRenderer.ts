/**
 * Renderer for Eckart potential.
 * V(x) = V_0 / (1 + exp(x/a))² - V_1 / (1 + exp(x/a))
 * This potential has an asymmetric shape and is used to model molecular barriers
 */

import { Shape } from "scenerystack/kite";
import QuantumConstants from "../../model/QuantumConstants.js";
import { hasBarrierHeight } from "../../model/ModelTypeGuards.js";
import { PotentialRenderer, RenderContext } from "./PotentialRenderer.js";

export class EckartRenderer extends PotentialRenderer {
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
      const expVal = Math.exp(dx / wellWidth);
      const denom = 1 + expVal;
      const V = wellDepth / (denom * denom) - barrierHeight / denom;

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
