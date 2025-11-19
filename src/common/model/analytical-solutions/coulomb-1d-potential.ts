/**
 * Analytical solution for the 1D Coulomb potential.
 * V(x) = -α/|x|
 *
 * This potential has a singularity at x=0 and describes a 1D hydrogen-like atom.
 */

import QuantumConstants from "../QuantumConstants.js";
import { BoundStateResult, GridConfig } from "../PotentialFunction.js";
import { associatedLaguerre } from "./math-utilities.js";

/**
 * Analytical solution for the 1D Coulomb potential.
 * V(x) = -α/|x|
 *
 * This potential has a singularity at x=0 and describes a 1D hydrogen-like atom.
 * The energy eigenvalues are given by E_n = -mα²/(2ℏ²(n+1/2)²)
 *
 * @param coulombStrength - Coulomb strength parameter α in J·m
 * @param mass - Particle mass in kg
 * @param numStates - Number of energy levels to calculate
 * @param gridConfig - Grid configuration for wavefunction evaluation
 * @returns Bound state results with exact energies and wavefunctions
 */
export function solveCoulomb1DPotential(
  coulombStrength: number,
  mass: number,
  numStates: number,
  gridConfig: GridConfig,
): BoundStateResult {
  const { HBAR } = QuantumConstants;
  const alpha = coulombStrength;

  // Calculate energies: E_n = -mα²/(2ℏ²(n+1/2)²) for n = 0, 1, 2, ...
  const energies: number[] = [];
  for (let n = 0; n < numStates; n++) {
    const energy = -(mass * alpha * alpha) / (2 * HBAR * HBAR * (n + 0.5) * (n + 0.5));
    energies.push(energy);
  }

  // Generate grid
  const numPoints = gridConfig.numPoints;
  const xGrid: number[] = [];
  const dx = (gridConfig.xMax - gridConfig.xMin) / (numPoints - 1);
  for (let i = 0; i < numPoints; i++) {
    xGrid.push(gridConfig.xMin + i * dx);
  }

  // Calculate effective Bohr radius for 1D: a_0 = ℏ²/(mα)
  const a0 = (HBAR * HBAR) / (mass * alpha);

  // Calculate wavefunctions
  // For 1D Coulomb, we use a hydrogen-like form with modified quantum numbers
  // ψ_n(x) ∝ exp(-|x|/n*a_0) * L_n(2|x|/(n*a_0))
  // where the effective n is (n + 1/2) for the 1D case
  const wavefunctions: number[][] = [];

  for (let n = 0; n < numStates; n++) {
    const wavefunction: number[] = [];

    // Effective principal quantum number for 1D
    const nEff = n + 0.5;
    const a_n = nEff * a0;

    // Normalization constant for 1D Coulomb
    // For ψ_n(x) = N * exp(-|x|/a_n) * L_n^1(2|x|/a_n)
    // with ∫_{-∞}^{∞} |ψ|² dx = 1
    // Using ∫_0^∞ exp(-ρ) * |L_n^1(ρ)|² dρ = (n+1)
    // Full integral: 2 * N² * (a_n/2) * (n+1) = N² * a_n * (n+1) = 1
    const normalization = Math.sqrt(1.0 / (a_n * (n + 1)));

    for (const x of xGrid) {
      const absX = Math.abs(x);
      const rho = 2 * absX / a_n;

      // Wavefunction: ψ(x) = N * exp(-ρ/2) * L_n^1(ρ)
      const laguerre = associatedLaguerre(n, 1, rho);
      const value = normalization * Math.exp(-rho / 2) * laguerre;

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
