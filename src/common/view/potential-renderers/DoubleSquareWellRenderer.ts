/**
 * Renderer for double square well potential.
 * Convention: V=0 in wells, V=wellDepth in barrier
 * Wells are centered at ±(separation/2 + wellWidth/2)
 */

import { Shape } from "scenerystack/kite";
import { isTwoWellsModel } from "../../model/ModelTypeGuards.js";
import { PotentialRenderer, RenderContext } from "./PotentialRenderer.js";

export class DoubleSquareWellRenderer extends PotentialRenderer {
  public render(context: RenderContext): Shape {
    const { wellWidth, wellDepth, dataToViewX, dataToViewY, chartMargins, chartWidth, model } = context;
    const shape = new Shape();

    // DOUBLE_SQUARE_WELL requires TwoWellsModel
    if (!isTwoWellsModel(model)) {
      // Return empty shape if model doesn't support this potential
      return shape;
    }

    const separation = model.wellSeparationProperty.value;

    const leftWellCenter = -(separation / 2 + wellWidth / 2);
    const rightWellCenter = separation / 2 + wellWidth / 2;
    const halfWidth = wellWidth / 2;

    // Well boundaries
    const leftWellLeft = leftWellCenter - halfWidth;
    const leftWellRight = leftWellCenter + halfWidth;
    const rightWellLeft = rightWellCenter - halfWidth;
    const rightWellRight = rightWellCenter + halfWidth;

    const y0 = dataToViewY(0); // Well energy
    const yBarrier = dataToViewY(wellDepth); // Barrier energy

    // Draw from left to right
    // Left outside region (at barrier height)
    shape.moveTo(chartMargins.left, yBarrier);
    shape.lineTo(dataToViewX(leftWellLeft), yBarrier);

    // Left well (drop down to V=0)
    shape.lineTo(dataToViewX(leftWellLeft), y0);
    shape.lineTo(dataToViewX(leftWellRight), y0);
    shape.lineTo(dataToViewX(leftWellRight), yBarrier);

    // Barrier between wells
    shape.lineTo(dataToViewX(rightWellLeft), yBarrier);

    // Right well (drop down to V=0)
    shape.lineTo(dataToViewX(rightWellLeft), y0);
    shape.lineTo(dataToViewX(rightWellRight), y0);
    shape.lineTo(dataToViewX(rightWellRight), yBarrier);

    // Right outside region (at barrier height)
    shape.lineTo(chartWidth - chartMargins.right, yBarrier);

    return shape;
  }
}
