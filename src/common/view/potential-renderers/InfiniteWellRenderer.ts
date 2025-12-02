/**
 * Renderer for infinite square well potential.
 * V(x) = 0 inside the well [-wellWidth/2, wellWidth/2]
 * V(x) = ∞ outside (displayed as 15 eV)
 */

import { Shape } from "scenerystack/kite";
import { PotentialRenderer, RenderContext } from "./PotentialRenderer.js";

export class InfiniteWellRenderer extends PotentialRenderer {
  public render(context: RenderContext): Shape {
    const { wellWidth, dataToViewX, dataToViewY, chartMargins, chartWidth } = context;
    const shape = new Shape();

    // Draw square well centered at x=0
    // Well extends from -wellWidth/2 to +wellWidth/2
    // V=0 inside, V=15 eV outside (representing infinity)
    const x1 = dataToViewX(-wellWidth / 2);
    const x2 = dataToViewX(wellWidth / 2);
    const y0 = dataToViewY(0);
    const y15 = dataToViewY(15); // Display infinity as 15 eV

    // Left region (x < -wellWidth/2): horizontal line at 15 eV
    shape.moveTo(chartMargins.left, y15);
    shape.lineTo(x1, y15);

    // Left wall: vertical line from 15 eV down to 0 eV
    shape.lineTo(x1, y0);

    // Bottom (inside well): horizontal line at 0 eV
    shape.lineTo(x2, y0);

    // Right wall: vertical line from 0 eV up to 15 eV
    shape.lineTo(x2, y15);

    // Right region (x > wellWidth/2): horizontal line at 15 eV
    shape.lineTo(chartWidth - chartMargins.right, y15);

    return shape;
  }
}
