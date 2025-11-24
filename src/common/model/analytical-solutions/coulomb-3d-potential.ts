/**
 * Analytical solution for the 3D Coulomb potential (hydrogen atom, radial equation with L=0).
 * V(r) = -α/r
 *
 * This solves the radial Schrödinger equation for the hydrogen atom with L=0 (s-waves).
 *
 * REFERENCES:
 * - Griffiths, D. J., & Schroeter, D. F. (2018). "Introduction to Quantum Mechanics" (3rd ed.).
 *   Cambridge University Press. Section 4.2, pp. 139-150.
 *   https://doi.org/10.1017/9781316995433
 *   Complete derivation of hydrogen atom wavefunctions using series solutions.
 *
 * - Bethe, H. A., & Salpeter, E. E. (1957). "Quantum Mechanics of One- and Two-Electron Atoms".
 *   Springer. Section 1-3, pp. 1-35.
 *   https://doi.org/10.1007/978-3-662-12869-5
 *   Classic comprehensive treatment of hydrogen-like atoms.
 *
 * - Schiff, L. I. (1968). "Quantum Mechanics" (3rd ed.). McGraw-Hill.
 *   Section 16, pp. 88-97.
 *   Derivation using separation of variables and associated Laguerre polynomials.
 *
 * - Pauli, W. (1926). "Über das Wasserstoffspektrum vom Standpunkt der neuen Quantenmechanik"
 *   Zeitschrift für Physik, 36(5), 336-363.
 *   https://doi.org/10.1007/BF01450175
 *   One of the first solutions of the hydrogen atom using Schrödinger's equation.
 *
 * - Abramowitz, M., & Stegun, I. A. (1964). "Handbook of Mathematical Functions".
 *   National Bureau of Standards. Section 13.4, pp. 509-515 and Section 22, pp. 773-802.
 *   https://doi.org/10.1119/1.15378
 *   Properties of associated Laguerre polynomials.
 *
 * ENERGY EIGENVALUES:
 *   E_n = -mα²/(2ℏ²n²),  n = 1, 2, 3, ...
 *   For electron in hydrogen: E_n = -13.6 eV/n² (Rydberg formula)
 *
 * RADIAL WAVEFUNCTIONS (L=0, s-states):
 *   R_{n0}(r) = N_{n0} · (2/na₀)^(3/2) · exp(-r/na₀) · L^1_{n-1}(2r/na₀)
 *   where a₀ = ℏ²/(mα) is the Bohr radius and L^1 are associated Laguerre polynomials
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
 * @param gridConfig - Grid configuration for wavefunction evaluation (uses |x| for radial symmetry)
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

  // Generate grid (can include negative values; we use |r| for radial symmetry)
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
    // N_n0 = sqrt{(2/(na₀))³ * (n-1)!/[2n(n!)³]}
    // Simplified form that ensures ∫|R_nl|²r²dr = 1
    const normalization =
      (2.0 / (a_n * Math.sqrt(a_n))) *
      Math.sqrt(factorial(n - 1) / (n * factorial(n)));

    for (const r of xGrid) {
      // For visualization on a symmetric x-axis, use |r| so the wavefunction appears on both sides
      const r_abs = Math.abs(r);

      if (r_abs === 0) {
        // At r=0, use a small offset to avoid singularity
        wavefunction.push(0);
        continue;
      }

      const rho = (2 * r_abs) / a_n;

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
