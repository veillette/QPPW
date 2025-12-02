/**
 * Renderer for triangular potential (finite well version).
 * V(x) = height + offset for x < 0
 * V(x) = offset at x = 0
 * V(x) = offset + (height/width) * x for 0 < x < width
 * V(x) = height + offset for x > width
 */

import { Shape } from "scenerystack/kite";
import { hasPotentialOffset } from "../../model/ModelTypeGuards.js";
import { PotentialRenderer, RenderContext } from "./PotentialRenderer.js";

export class TriangularRenderer extends PotentialRenderer {
  public render(context: RenderContext): Shape {
    const { wellWidth, wellDepth, dataToViewX, dataToViewY, chartMargins, chartWidth, model } = context;
    const shape = new Shape();

    // Get the offset from potentialOffsetProperty (OneWellModel or IntroModel only)
    const offset = hasPotentialOffset(model)
      ? model.potentialOffsetProperty.value
      : 0;
    const height = wellDepth;
    const barrierTop = height + offset;
    const slope = height / wellWidth; // eV/nm

    const yBarrier = dataToViewY(barrierTop);
    const yOffset = dataToViewY(offset);

    // Left region (x < 0): horizontal line at height + offset
    shape.moveTo(chartMargins.left, yBarrier);
    shape.lineTo(dataToViewX(0), yBarrier);

    // Drop to offset at x = 0
    shape.lineTo(dataToViewX(0), yOffset);

    // Linear region (0 < x < width): V = offset + slope * x
    const numPoints = 50;
    for (let i = 1; i <= numPoints; i++) {
      const x = (wellWidth * i) / numPoints;
      const V = offset + slope * x;
      const viewX = dataToViewX(x);
      const viewY = dataToViewY(V);
      shape.lineTo(viewX, viewY);
    }

    // At x = width, we should be back at height + offset
    // Then continue as horizontal line to the right
    shape.lineTo(chartWidth - chartMargins.right, yBarrier);

    return shape;
  }
}
