/**
 * Analytical solution for the asymmetric triangle potential.
 * V(x) = ∞ for x < 0 (infinite wall)
 * V(x) = F·x for x ≥ 0 (linear increasing potential)
 *
 * This is the standard triangular well problem with an infinite wall at x=0.
 * The eigenvalues are related to the zeros of the Airy function Ai(z).
 *
 * REFERENCES:
 * - Griffiths, D. J., & Schroeter, D. F. (2018). "Introduction to Quantum Mechanics" (3rd ed.).
 *   Cambridge University Press. Problem 2.43, p. 89.
 *   https://doi.org/10.1017/9781316995433
 *   Linear potential with hard wall boundary condition.
 *
 * - Landau, L. D., & Lifshitz, E. M. (1977). "Quantum Mechanics: Non-Relativistic Theory" (3rd ed.).
 *   Pergamon Press. Section 25, pp. 78-81.
 *   Linear potential and WKB approximation.
 *
 * - Schiff, L. I. (1968). "Quantum Mechanics" (3rd ed.). McGraw-Hill.
 *   Problem 14, pp. 269-270.
 *   Exact solution using Airy functions.
 *
 * - Vallée, O., & Soares, M. (2004). "Airy Functions and Applications to Physics".
 *   Imperial College Press. Chapter 5, pp. 115-145.
 *   https://doi.org/10.1142/p345
 *   Comprehensive treatment of Airy functions in quantum mechanics.
 *
 * - Abramowitz, M., & Stegun, I. A. (1964). "Handbook of Mathematical Functions".
 *   National Bureau of Standards. Section 10.4, pp. 446-452; Table 10.13, p. 478.
 *   https://doi.org/10.1119/1.15378
 *   Zeros of Airy function Ai(z): z_n for n = 1, 2, 3, ...
 *
 * - Fröman, N., & Fröman, P. O. (1965). "JWKB Approximation: Contributions to the Theory".
 *   North-Holland Publishing Company.
 *   Connection formulas for linear turning points.
 *
 * ENERGY EIGENVALUES (exact):
 *   E_n = (ℏ²/2m)^(1/3) · F^(2/3) · |z_n|
 * where z_n is the n-th zero of the Airy function Ai(z) (all negative):
 *   z_1 ≈ -2.338107, z_2 ≈ -4.087949, z_3 ≈ -5.520560, ...
 *
 * WAVEFUNCTIONS (exact):
 *   ψ_n(x) = N_n · Ai(α(x - x_n))
 * where α = (2mF/ℏ²)^(1/3), x_n = E_n/F is the classical turning point,
 * and N_n is the normalization constant.
 *
 * Boundary condition: ψ(0) = 0 leads to Ai(-αx_n) = 0, giving αx_n = -z_n.
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
