/**
 * Renderer for finite square well potential.
 * V(x) = 0 at infinity
 * V(x) = -wellDepth inside the well [xCenter - wellWidth/2, xCenter + wellWidth/2]
 */

import { Shape } from "scenerystack/kite";
import { PotentialRenderer, RenderContext } from "./PotentialRenderer.js";

export class FiniteWellRenderer extends PotentialRenderer {
  public render(context: RenderContext): Shape {
    const { xCenter, wellWidth, wellDepth, dataToViewX, dataToViewY, chartMargins, chartWidth } = context;
    const shape = new Shape();

    // Draw finite square well centered at xCenter
    // V=0 at infinity, V=-wellDepth inside the well
    const x1 = dataToViewX(xCenter - wellWidth / 2);
    const x2 = dataToViewX(xCenter + wellWidth / 2);
    const y0 = dataToViewY(0);
    const yDepth = dataToViewY(-wellDepth);

    shape.moveTo(chartMargins.left, y0);
    shape.lineTo(x1, y0);
    shape.lineTo(x1, yDepth);
    shape.lineTo(x2, yDepth);
    shape.lineTo(x2, y0);
    shape.lineTo(chartWidth - chartMargins.right, y0);

    return shape;
  }
}
