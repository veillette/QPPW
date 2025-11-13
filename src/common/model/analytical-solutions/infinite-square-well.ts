/**
 * Analytical solution for an infinite square well (particle in a box).
 * V(x) = 0 for 0 < x < L, V(x) = ∞ otherwise
 */

import QuantumConstants from "../QuantumConstants.js";
import { BoundStateResult, GridConfig } from "../PotentialFunction.js";

/**
 * Analytical solution for an infinite square well (particle in a box).
 * V(x) = 0 for 0 < x < L, V(x) = ∞ otherwise
 *
 * @param wellWidth - Width of the well (L) in meters
 * @param mass - Particle mass in kg
 * @param numStates - Number of energy levels to calculate
 * @param gridConfig - Grid configuration for wavefunction evaluation
 * @returns Bound state results with exact energies and wavefunctions
 */
export function solveInfiniteWell(
  wellWidth: number,
  mass: number,
  numStates: number,
  gridConfig: GridConfig,
): BoundStateResult {
  const { HBAR } = QuantumConstants;
  const L = wellWidth;

  // Calculate energies: E_n = (n^2 * π^2 * ℏ^2) / (2 * m * L^2) for n = 1, 2, 3, ...
  const energies: number[] = [];
  for (let n = 1; n <= numStates; n++) {
    const energy =
      (n * n * Math.PI * Math.PI * HBAR * HBAR) / (2 * mass * L * L);
    energies.push(energy);
  }

  // Generate grid
  const numPoints = gridConfig.numPoints;
  const xGrid: number[] = [];
  const dx = (gridConfig.xMax - gridConfig.xMin) / (numPoints - 1);
  for (let i = 0; i < numPoints; i++) {
    xGrid.push(gridConfig.xMin + i * dx);
  }

  // Calculate wavefunctions: ψ_n(x) = sqrt(2/L) * sin(n * π * x / L)
  const wavefunctions: number[][] = [];
  for (let n = 1; n <= numStates; n++) {
    const wavefunction: number[] = [];
    const normalization = Math.sqrt(2 / L);

    for (const x of xGrid) {
      // Well is centered, so shift x coordinate
      const xShifted = x - gridConfig.xMin;
      if (xShifted >= 0 && xShifted <= L) {
        const value = normalization * Math.sin((n * Math.PI * xShifted) / L);
        wavefunction.push(value);
      } else {
        wavefunction.push(0);
      }
    }
    wavefunctions.push(wavefunction);
  }

  return {
    energies,
    wavefunctions,
    xGrid,
    method: "analytical",
  };
}
