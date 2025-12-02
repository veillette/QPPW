/**
 * Renderer for asymmetric triangle potential (infinite wall version).
 * V(x) = ∞ for x < 0 (displayed as 15 eV)
 * V(x) = F·x for x ≥ 0 (linear increasing)
 * where F = wellDepth/wellWidth (slope)
 */

import { Shape } from "scenerystack/kite";
import QuantumConstants from "../../model/QuantumConstants.js";
import { PotentialRenderer, RenderContext } from "./PotentialRenderer.js";

export class AsymmetricTriangleRenderer extends PotentialRenderer {
  public render(context: RenderContext): Shape {
    const { xGrid, wellWidth, wellDepth, dataToViewX, dataToViewY, chartMargins } = context;
    const shape = new Shape();

    const F_eV_per_nm = wellDepth / wellWidth; // slope in eV/nm
    const y15eV = dataToViewY(15); // Display infinity as 15 eV
    const y0eV = dataToViewY(0);

    // Left region (x < 0): vertical line at x=0 representing infinite wall
    shape.moveTo(chartMargins.left, y15eV);
    shape.lineTo(dataToViewX(0), y15eV);
    shape.lineTo(dataToViewX(0), y0eV);

    // Right region (x ≥ 0): linear increasing potential V = F·x
    const numPoints = 100;
    for (let i = 0; i <= numPoints; i++) {
      const x =
        ((xGrid[xGrid.length - 1] * i) / numPoints) * QuantumConstants.M_TO_NM;
      if (x >= 0) {
        const V = F_eV_per_nm * x;
        const viewX = dataToViewX(x);
        const viewY = dataToViewY(Math.min(V, 15)); // Clamp to chart range
        shape.lineTo(viewX, viewY);
      }
    }

    return shape;
  }
}
