/**
 * Analytical solution for the asymmetric triangle potential.
 * V(x) = ∞ for x < 0 (infinite wall)
 * V(x) = F·x for x ≥ 0 (linear increasing potential)
 *
 * This is the standard triangular well problem with an infinite wall at x=0.
 * The eigenvalues are related to the zeros of the Airy function Ai(z).
 */

import { BoundStateResult, GridConfig } from "../PotentialFunction.js";
import { airyAi } from "./math-utilities.js";
import {
  calculateAiryAlpha,
  getAiryZero,
  calculateTriangularWellEnergy,
  generateGrid,
  normalizeWavefunction,
} from "./airy-utilities.js";

/**
 * Analytical solution for the asymmetric triangle potential with infinite wall.
 *
 * @param slope - Slope parameter F in Joules/meter (field strength)
 * @param _wellWidth - Width parameter (not used for infinite well, kept for API compatibility)
 * @param mass - Particle mass in kg
 * @param numStates - Number of energy levels to calculate
 * @param gridConfig - Grid configuration for wavefunction evaluation
 * @returns Bound state results with exact energies and wavefunctions
 */
export function solveAsymmetricTrianglePotential(
  slope: number,
  _wellWidth: number,
  mass: number,
  numStates: number,
  gridConfig: GridConfig,
): BoundStateResult {
  const F = slope;
  const alpha = calculateAiryAlpha(mass, F);

  // Calculate energies using Airy zeros
  const energies: number[] = [];
  for (let n = 0; n < numStates; n++) {
    const z_n = getAiryZero(n);
    const energy = calculateTriangularWellEnergy(z_n, mass, F);
    energies.push(energy);
  }

  const actualNumStates = energies.length;

  // Generate grid
  const { xGrid, dx } = generateGrid(gridConfig);

  // Calculate wavefunctions using Airy functions
  const wavefunctions: number[][] = [];

  for (let n = 0; n < actualNumStates; n++) {
    const E = energies[n];

    // Classical turning point: x_0 = E/F (where V(x_0) = F·x_0 = E)
    const x0 = E / F;

    // Calculate unnormalized wavefunction
    const psiRaw: number[] = [];
    for (const x of xGrid) {
      if (x < 0) {
        // Region x < 0: infinite wall, ψ = 0
        psiRaw.push(0);
      } else {
        // Region x ≥ 0: Airy function solution
        // ψ(x) = N · Ai(α(x - x_0))
        const z = alpha * (x - x0);
        psiRaw.push(airyAi(z));
      }
    }

    // Normalize wavefunction
    const wavefunction = normalizeWavefunction(psiRaw, dx);
    wavefunctions.push(wavefunction);
  }

  return {
    energies,
    wavefunctions,
    xGrid,
    method: "analytical",
  };
}
