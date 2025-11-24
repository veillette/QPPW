/**
 * Analytical solution for the Morse potential.
 * V(x) = D_e * (1 - exp(-(x - x_e)/a))^2
 *
 * The Morse potential describes molecular vibrations more accurately than the harmonic oscillator
 * by including anharmonic effects and bond dissociation.
 *
 * REFERENCES:
 * - Morse, P. M. (1929). "Diatomic Molecules According to the Wave Mechanics. II. Vibrational Levels".
 *   Physical Review, 34(1), 57-64.
 *   https://doi.org/10.1103/PhysRev.34.57
 *   ORIGINAL PAPER: Introduced the Morse potential and derived exact solutions.
 *
 * - Flügge, S. (1999). "Practical Quantum Mechanics". Springer.
 *   Problem 38, pp. 94-95. https://doi.org/10.1007/978-3-642-61995-3
 *   Detailed solution procedure using associated Laguerre polynomials.
 *
 * - Cooper, I. L. (1993). "An Accurate Analytic Solution of the Morse Oscillator Problem".
 *   Journal of Chemical Education, 70(11), 887.
 *   https://doi.org/10.1021/ed070p887
 *   Pedagogical treatment with explicit wavefunctions.
 *
 * - Dahl, J. P., & Springborg, M. (1988). "The Morse oscillator in position space, momentum space,
 *   and phase space". Journal of Chemical Physics, 88(7), 4535-4547.
 *   https://doi.org/10.1063/1.453761
 *   Complete phase-space analysis of Morse oscillator.
 *
 * ENERGY EIGENVALUES:
 *   E_n = ℏω(n + 1/2) - (ℏω)²(n + 1/2)²/(4D_e) - D_e,  n = 0, 1, 2, ..., n_max
 *   where ω = (1/a)√(2D_e/m) and n_max = floor(a√(2mD_e)/ℏ - 1/2)
 *
 * WAVEFUNCTIONS:
 *   ψ_n(z) = N_n · z^(λ-n-1/2) · exp(-z/2) · L_n^(2λ-2n-1)(z)
 *   where z = 2λ exp(-(x-x_e)/a), λ = a√(2mD_e)/ℏ
 */

import QuantumConstants from "../QuantumConstants.js";
import { BoundStateResult, GridConfig } from "../PotentialFunction.js";
import { associatedLaguerre, factorial, gamma } from "./math-utilities.js";

/**
 * Analytical solution for the Morse potential.
 * V(x) = D_e * (1 - exp(-(x - x_e)/a))^2
 *
 * The Morse potential describes molecular vibrations more accurately than the harmonic oscillator
 * by including anharmonic effects and bond dissociation.
 *
 * @param dissociationEnergy - Dissociation energy D_e in Joules
 * @param wellWidth - Width parameter a in meters
 * @param equilibriumPosition - Equilibrium position x_e in meters
 * @param mass - Particle mass in kg
 * @param numStates - Number of energy levels to calculate
 * @param gridConfig - Grid configuration for wavefunction evaluation
 * @returns Bound state results with exact energies and wavefunctions
 */
export function solveMorsePotential(
  dissociationEnergy: number,
  wellWidth: number,
  equilibriumPosition: number,
  mass: number,
  numStates: number,
  gridConfig: GridConfig,
): BoundStateResult {
  const { HBAR } = QuantumConstants;
  const De = dissociationEnergy;
  const a = wellWidth;
  const xe = equilibriumPosition;

  // Calculate the maximum quantum number
  // With substitution a_old = 1/a_new:
  // n_max = floor(a * sqrt(2*m*D_e)/ℏ - 1/2)
  const lambda = (a * Math.sqrt(2 * mass * De)) / HBAR;
  const nMax = Math.floor(lambda - 0.5);
  const actualNumStates = Math.min(numStates, nMax + 1);

  if (actualNumStates <= 0) {
    throw new Error("Morse potential too shallow to support bound states");
  }

  // Calculate the characteristic frequency
  // ω = sqrt(2*D_e/m) / a
  const omega = Math.sqrt((2 * De) / mass) / a;

  // Calculate energies: E_n = ℏω(n + 1/2) - (ℏω)²(n + 1/2)² / (4*D_e)
  // Relative to the bottom of the well
  const energies: number[] = [];
  for (let n = 0; n < actualNumStates; n++) {
    const term1 = HBAR * omega * (n + 0.5);
    const term2 =
      (HBAR * HBAR * omega * omega * (n + 0.5) * (n + 0.5)) / (4 * De);
    const energy = term1 - term2 - De; // Energy relative to dissociation limit
    energies.push(energy);
  }

  // Generate grid
  const numPoints = gridConfig.numPoints;
  const xGrid: number[] = [];
  const dx = (gridConfig.xMax - gridConfig.xMin) / (numPoints - 1);
  for (let i = 0; i < numPoints; i++) {
    xGrid.push(gridConfig.xMin + i * dx);
  }

  // Calculate wavefunctions using associated Laguerre polynomials
  // ψ_n(z) = N_n * z^((λ-n-1/2)) * exp(-z/2) * L_n^(2λ-2n-1)(z)
  // where z = 2λ * exp(-(x-xe)/a), λ = a * sqrt(2*m*D_e)/ℏ
  const wavefunctions: number[][] = [];

  for (let n = 0; n < actualNumStates; n++) {
    const wavefunction: number[] = [];

    // Normalization constant (includes 1/a factor from variable change)
    const alpha = 2 * lambda - 2 * n - 1;
    const normalization = Math.sqrt(factorial(n) / a / gamma(2 * lambda - n));

    for (const x of xGrid) {
      const z = 2 * lambda * Math.exp(-(x - xe) / a);

      // Calculate wavefunction
      const exponent = lambda - n - 0.5;
      const laguerre = associatedLaguerre(n, alpha, z);
      const value =
        normalization * Math.pow(z, exponent) * Math.exp(-z / 2) * laguerre;

      wavefunction.push(value);
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
