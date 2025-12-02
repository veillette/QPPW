/**
 * Base class for potential energy curve renderers.
 * Each potential type has its own renderer that knows how to draw that specific potential.
 * This follows the Strategy pattern to decouple rendering logic from the chart.
 */

import { Shape } from "scenerystack/kite";
import type { ScreenModel } from "../../model/ScreenModels.js";

/**
 * Context object containing all the information needed to render a potential curve.
 */
export interface RenderContext {
  // Grid and model data
  xGrid: number[]; // Position grid in meters
  xCenter: number; // Center of the grid in nanometers
  wellWidth: number; // Well width in nanometers
  wellDepth: number; // Well depth in electron volts

  // Coordinate transformation functions
  dataToViewX: (x: number) => number; // Convert data x (nm) to view x (pixels)
  dataToViewY: (y: number) => number; // Convert data y (eV) to view y (pixels)

  // Chart dimensions (for drawing outside plot area)
  chartMargins: { left: number; right: number; top: number; bottom: number };
  chartWidth: number;
  chartHeight: number;

  // Model reference for accessing additional properties
  model: ScreenModel;
}

/**
 * Abstract base class for potential renderers.
 * Each concrete renderer implements the render() method for a specific potential type.
 */
export abstract class PotentialRenderer {
  /**
   * Render the potential curve and return the Shape to be drawn.
   * @param context - All the information needed to render the potential
   * @returns A Shape object representing the potential curve
   */
  public abstract render(context: RenderContext): Shape;

  /**
   * Helper method to iterate over x grid with a specified number of points.
   * Useful for smooth curves like harmonic oscillator, Morse, etc.
   */
  protected iterateXGrid(
    xGrid: number[],
    numPoints: number,
    callback: (x: number, i: number) => void,
  ): void {
    for (let i = 0; i < numPoints; i++) {
      const x =
        (xGrid[0] + ((xGrid[xGrid.length - 1] - xGrid[0]) * i) / (numPoints - 1));
      callback(x, i);
    }
  }
}
