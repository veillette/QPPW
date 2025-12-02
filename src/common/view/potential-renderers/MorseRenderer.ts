/**
 * Renderer for Morse potential.
 * V(x) = D_e * (1 - exp(-(x-x_e)/a))² - D_e
 * With x_e = 0 (centered), D_e = wellDepth, a = wellWidth
 * V=0 at dissociation (x→∞), V=-D_e at bottom of well (x=x_e)
 */

import { Shape } from "scenerystack/kite";
import QuantumConstants from "../../model/QuantumConstants.js";
import { PotentialRenderer, RenderContext } from "./PotentialRenderer.js";

export class MorseRenderer extends PotentialRenderer {
  public render(context: RenderContext): Shape {
    const { xGrid, xCenter, wellWidth, wellDepth, dataToViewX, dataToViewY } = context;
    const shape = new Shape();

    const centerX = xCenter;
    const numPoints = 200;
    let firstPoint = true;

    this.iterateXGrid(xGrid, numPoints, (xMeters) => {
      const x = xMeters * QuantumConstants.M_TO_NM;
      const dx = x - centerX;
      const exponent = Math.exp(-dx / wellWidth);
      // V(x) = D_e * (1 - e^(-dx/a))^2 - D_e
      // At x=x_e (dx=0): V = D_e * (1-1)^2 - D_e = -D_e (bottom of well)
      // At x→∞: V = D_e * (1-0)^2 - D_e = 0 (dissociation limit)
      const V = wellDepth * Math.pow(1 - exponent, 2) - wellDepth;

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
