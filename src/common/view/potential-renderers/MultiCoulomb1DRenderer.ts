/**
 * Renderer for multi-Coulomb 1D potential (multiple Coulomb centers).
 * V(x) = Σ(-k/|x - x_i|) where x_i are the center positions
 * Supports electric field tilt: V_field(x) = -E * x
 */

import { Shape } from "scenerystack/kite";
import QuantumConstants from "../../model/QuantumConstants.js";
import { PotentialRenderer, RenderContext } from "./PotentialRenderer.js";

export class MultiCoulomb1DRenderer extends PotentialRenderer {
  public render(context: RenderContext): Shape {
    const { xGrid, wellWidth, wellDepth, dataToViewX, dataToViewY, model } = context;
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

    const numPoints = 200;
    let firstPoint = true;

    // Calculate positions of Coulomb centers
    // Ensure they fit within the visible range of -4nm to 4nm
    const maxVisibleRange = 7.5; // Use 7.5nm of the 8nm available range

    // Start with base scaling
    let effectiveSeparation = separationParam * 10;

    // Calculate natural span
    const naturalSpan = (numberOfWells - 1) * effectiveSeparation;

    // Scale down if needed to fit within visible range
    if (naturalSpan > maxVisibleRange) {
      effectiveSeparation = maxVisibleRange / (numberOfWells - 1);
    }

    const totalSpan = (numberOfWells - 1) * effectiveSeparation;
    const centers: number[] = [];
    for (let i = 0; i < numberOfWells; i++) {
      centers.push(-totalSpan / 2 + i * effectiveSeparation);
    }

    // Scale factor for Coulomb strength
    const k = wellDepth * (wellWidth / 2);

    this.iterateXGrid(xGrid, numPoints, (xMeters) => {
      const x = xMeters * QuantumConstants.M_TO_NM;

      // Sum contributions from all Coulomb centers
      let V = 0;
      for (const center of centers) {
        const dx = x - center;
        const distance = Math.max(Math.abs(dx), 0.01); // Avoid singularity
        V += -k / distance;
      }

      // Add electric field tilt: V_field(x) = -E * x
      V += -electricField * x;

      // Clamp V to reasonable range for display
      const VClamped = Math.max(V, -20); // Reasonable default clamp value

      const viewX = dataToViewX(x);
      const viewY = dataToViewY(VClamped);

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
