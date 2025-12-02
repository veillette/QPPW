/**
 * Renderer for multi-square well potential (generalization of double square well).
 * Convention: V=0 in wells, V=wellDepth in barrier
 * Supports electric field tilt: V_field(x) = -E * x
 */

import { Shape } from "scenerystack/kite";
import { PotentialRenderer, RenderContext } from "./PotentialRenderer.js";

export class MultiSquareWellRenderer extends PotentialRenderer {
  public render(context: RenderContext): Shape {
    const { wellWidth, wellDepth, dataToViewX, dataToViewY, chartMargins, chartWidth, model } = context;
    const shape = new Shape();

    // Access ManyWellsModel properties
    const manyWellsModel = model as any; // Type cast to access dynamic properties
    const numberOfWells =
      "numberOfWellsProperty" in manyWellsModel
        ? manyWellsModel.numberOfWellsProperty.value
        : 3;
    const separationParam =
      "wellSeparationProperty" in manyWellsModel
        ? manyWellsModel.wellSeparationProperty.value
        : 0.2;
    const electricField =
      "electricFieldProperty" in manyWellsModel
        ? manyWellsModel.electricFieldProperty.value
        : 0.0;

    // Calculate the natural total span needed
    const maxVisibleRange = 7.5; // Use 7.5nm of the 8nm available range (leave small margins)

    // Start with a base scaling for separation
    let effectiveSeparation = separationParam * 10;
    let effectiveWellWidth = wellWidth;

    // Calculate what the total span would be
    const naturalSpan =
      numberOfWells * effectiveWellWidth +
      (numberOfWells - 1) * effectiveSeparation;

    // If it exceeds the visible range, scale everything down proportionally
    if (naturalSpan > maxVisibleRange) {
      const scaleFactor = maxVisibleRange / naturalSpan;
      effectiveSeparation *= scaleFactor;
      effectiveWellWidth *= scaleFactor;
    }

    // Calculate total span and center the wells
    const totalSpan =
      numberOfWells * effectiveWellWidth +
      (numberOfWells - 1) * effectiveSeparation;
    const startX = -totalSpan / 2;

    // Helper function to get potential with electric field tilt
    const getV = (x: number, baseV: number) => baseV - electricField * x;

    // Access xMinProperty and xMaxProperty from context (we need to add these to RenderContext)
    // For now, use reasonable defaults
    const leftmostX = -4; // Reasonable default for xMin
    const rightmostX = 4; // Reasonable default for xMax

    // Start from left at barrier height (with field tilt)
    shape.moveTo(
      chartMargins.left,
      dataToViewY(getV(leftmostX, wellDepth)),
    );

    // Draw each well
    for (let i = 0; i < numberOfWells; i++) {
      const wellLeft =
        startX + i * (effectiveWellWidth + effectiveSeparation);
      const wellRight = wellLeft + effectiveWellWidth;

      // Drop into well (with field tilt)
      shape.lineTo(
        dataToViewX(wellLeft),
        dataToViewY(getV(wellLeft, wellDepth)),
      );
      shape.lineTo(
        dataToViewX(wellLeft),
        dataToViewY(getV(wellLeft, 0)),
      );
      shape.lineTo(
        dataToViewX(wellRight),
        dataToViewY(getV(wellRight, 0)),
      );
      shape.lineTo(
        dataToViewX(wellRight),
        dataToViewY(getV(wellRight, wellDepth)),
      );

      // Barrier between wells (except after last well)
      if (i < numberOfWells - 1) {
        const barrierRight = wellRight + effectiveSeparation;
        shape.lineTo(
          dataToViewX(barrierRight),
          dataToViewY(getV(barrierRight, wellDepth)),
        );
      }
    }

    // Right outside region (with field tilt)
    shape.lineTo(
      chartWidth - chartMargins.right,
      dataToViewY(getV(rightmostX, wellDepth)),
    );

    return shape;
  }
}
