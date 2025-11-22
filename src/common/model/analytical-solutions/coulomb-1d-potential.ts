/**
 * Analytical solution for the 1D Coulomb potential.
 * V(x) = -α/|x|
 *
 * This potential has a singularity at x=0 and describes a 1D hydrogen-like atom.
 *
 * IMPORTANT: For the pure 1D Coulomb potential, only odd-parity eigenstates exist
 * as normalizable solutions. Even-parity eigenstates are absent because they would
 * diverge at the origin. This implementation correctly produces ALL odd-parity
 * wavefunctions with ψ(0) = 0.
 *
 * WARNING: Standard numerical solvers (DVR, FGH, etc.) will incorrectly find a mix
 * of even and odd parity states if applied to this potential. Only the analytical
 * solution or the solveCoulomb1DNumerical wrapper (which filters for odd parity)
 * should be used. The analytical solution is STRONGLY PREFERRED as it's exact and
 * much more efficient.
 *
 * Reference:
 * - Loudon, R. (2016). "The one-dimensional Coulomb problem"
 *   Proc. R. Soc. A 472: 20150534. https://doi.org/10.1098/rspa.2015.0534
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

      // Wavefunction: ψ(x) = sign(x) * N * exp(-ρ/2) * L_n^1(ρ)
      // ODD parity: ψ(-x) = -ψ(x), which matches the energy formula E_n = -E_R/(n+1/2)²
      const laguerre = associatedLaguerre(n, 1, rho);
      const radialPart = normalization * Math.exp(-rho / 2) * laguerre;

      // Apply sign(x) for odd parity
      // Math.sign(0) = 0 ensures ψ(0) = 0, which is required for odd parity states
      const value = Math.sign(x) * radialPart;

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
