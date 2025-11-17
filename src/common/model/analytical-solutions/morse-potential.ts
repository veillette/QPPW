/**
 * Analytical solution for the Morse potential.
 * V(x) = D_e * (1 - exp(-(x - x_e)/a))^2
 *
 * The Morse potential describes molecular vibrations more accurately than the harmonic oscillator
 * by including anharmonic effects and bond dissociation.
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
  const omega = Math.sqrt(2 * De / mass) / a;

  // Calculate energies: E_n = ℏω(n + 1/2) - (ℏω)²(n + 1/2)² / (4*D_e)
  // Relative to the bottom of the well
  const energies: number[] = [];
  for (let n = 0; n < actualNumStates; n++) {
    const term1 = HBAR * omega * (n + 0.5);
    const term2 = (HBAR * HBAR * omega * omega * (n + 0.5) * (n + 0.5)) / (4 * De);
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
    const normalization = Math.sqrt(
      (factorial(n) / a) / (gamma(2 * lambda - n))
    );

    for (const x of xGrid) {
      const z = 2 * lambda * Math.exp(-(x - xe) / a);

      // Calculate wavefunction
      const exponent = lambda - n - 0.5;
      const laguerre = associatedLaguerre(n, alpha, z);
      const value = normalization * Math.pow(z, exponent) * Math.exp(-z / 2) * laguerre;

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
