/**
 * Analytical solution for the 3D Coulomb potential (hydrogen atom, radial equation with L=0).
 * V(r) = -α/r
 *
 * This solves the radial Schrödinger equation for the hydrogen atom with L=0 (s-waves).
 */

import QuantumConstants from "../QuantumConstants.js";
import { BoundStateResult, GridConfig } from "../PotentialFunction.js";
import { associatedLaguerre, factorial } from "./math-utilities.js";

/**
 * Analytical solution for the 3D Coulomb potential (hydrogen atom, radial equation with L=0).
 * V(r) = -α/r
 *
 * This solves the radial Schrödinger equation for the hydrogen atom with L=0 (s-waves).
 * The energy eigenvalues are given by E_n = -mα²/(2ℏ²n²) for n = 1, 2, 3, ...
 *
 * @param coulombStrength - Coulomb strength parameter α in J·m
 * @param mass - Particle mass in kg
 * @param numStates - Number of energy levels to calculate
 * @param gridConfig - Grid configuration for radial wavefunction evaluation (r > 0)
 * @returns Bound state results with exact energies and radial wavefunctions
 */
export function solveCoulomb3DPotential(
  coulombStrength: number,
  mass: number,
  numStates: number,
  gridConfig: GridConfig,
): BoundStateResult {
  const { HBAR } = QuantumConstants;
  const alpha = coulombStrength;

  // Calculate energies: E_n = -mα²/(2ℏ²n²) for n = 1, 2, 3, ...
  const energies: number[] = [];
  for (let n = 1; n <= numStates; n++) {
    const energy = -(mass * alpha * alpha) / (2 * HBAR * HBAR * n * n);
    energies.push(energy);
  }

  // Generate grid (must be r > 0 for radial equation)
  const numPoints = gridConfig.numPoints;
  const xGrid: number[] = [];
  const dx = (gridConfig.xMax - gridConfig.xMin) / (numPoints - 1);
  for (let i = 0; i < numPoints; i++) {
    const r = gridConfig.xMin + i * dx;
    xGrid.push(r);
  }

  // Calculate Bohr radius: a_0 = ℏ²/(mα)
  const a0 = (HBAR * HBAR) / (mass * alpha);

  // Calculate radial wavefunctions for L=0
  // R_n0(r) = N_n0 * (2/na_0)^(3/2) * (2r/na_0)^0 * exp(-r/na_0) * L^1_(n-1)(2r/na_0)
  const wavefunctions: number[][] = [];

  for (let n = 1; n <= numStates; n++) {
    const wavefunction: number[] = [];

    const a_n = n * a0;

    // Normalization constant for L=0
    // N_n0 = 2/(n²a_0^3) * sqrt(1/n)
    const normalization = 2.0 / (a_n * Math.sqrt(a_n)) * Math.sqrt(factorial(n - 1) / (2 * n * factorial(n)));

    for (const r of xGrid) {
      if (r <= 0) {
        // Radial wavefunction must be zero at r=0 for L>0, but for L=0 it's finite
        // However, we'll handle r <= 0 by setting to 0
        wavefunction.push(0);
        continue;
      }

      const rho = 2 * r / a_n;

      // Radial wavefunction: R_n0(r) = N * exp(-ρ/2) * L^1_(n-1)(ρ)
      const laguerre = associatedLaguerre(n - 1, 1, rho);
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
