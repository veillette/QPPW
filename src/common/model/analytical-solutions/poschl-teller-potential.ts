/**
 * Analytical solution for the Pöschl-Teller potential.
 * V(x) = -V_0 / cosh²(x/a)
 *
 * This potential is useful for modeling quantum wells and has exact solutions.
 *
 * REFERENCES:
 * - Pöschl, G., & Teller, E. (1933). "Bemerkungen zur Quantenmechanik des anharmonischen Oszillators"
 *   Zeitschrift für Physik, 83(3-4), 143-151.
 *   https://doi.org/10.1007/BF01331132
 *   ORIGINAL PAPER: Introduced this potential and its exact solutions.
 *
 * - Flügge, S. (1999). "Practical Quantum Mechanics". Springer.
 *   Problem 39, pp. 95-97. https://doi.org/10.1007/978-3-642-61995-3
 *   Detailed solution using Jacobi polynomials.
 *
 * - Cooper, F., Khare, A., & Sukhatme, U. (1995). "Supersymmetry and quantum mechanics"
 *   Physics Reports, 251(5-6), 267-385.
 *   https://doi.org/10.1016/0370-1573(94)00080-M
 *   Section 3.2, pp. 281-283: Pöschl-Teller as a shape-invariant potential.
 *
 * - Natanzon, G. A. (1979). "General properties of potentials for which the Schrödinger equation
 *   can be solved by means of hypergeometric functions". Theoretical and Mathematical Physics, 38(2), 146-153.
 *   https://doi.org/10.1007/BF01016836
 *   Classification of exactly solvable potentials including Pöschl-Teller.
 *
 * ENERGY EIGENVALUES:
 *   E_n = -V_0 [(λ - n - 1/2)/λ]²,  n = 0, 1, 2, ..., n_max
 *   where λ = a√(2mV_0)/ℏ and n_max = floor(λ - 1/2)
 *
 * WAVEFUNCTIONS:
 *   ψ_n(x) = N_n · sech^(λ-n-1/2)(x/a) · P_n^(α,α)(tanh(x/a))
 *   where α = λ - n - 1/2 and P_n^(α,α) are Jacobi polynomials
 */

import QuantumConstants from "../QuantumConstants.js";
import { BoundStateResult, GridConfig } from "../PotentialFunction.js";
import { jacobiPolynomial, factorial } from "./math-utilities.js";

/**
 * Analytical solution for the Pöschl-Teller potential.
 * V(x) = -V_0 / cosh²(x/a)
 *
 * This potential is useful for modeling quantum wells and has exact solutions.
 *
 * @param potentialDepth - Potential depth V_0 in Joules (positive value)
 * @param wellWidth - Width parameter a in meters
 * @param mass - Particle mass in kg
 * @param numStates - Number of energy levels to calculate
 * @param gridConfig - Grid configuration for wavefunction evaluation
 * @returns Bound state results with exact energies and wavefunctions
 */
export function solvePoschlTellerPotential(
  potentialDepth: number,
  wellWidth: number,
  mass: number,
  numStates: number,
  gridConfig: GridConfig,
): BoundStateResult {
  const { HBAR } = QuantumConstants;
  const V0 = potentialDepth;
  const a = wellWidth;

  // Calculate λ = a * sqrt(2*m*V_0) / ℏ
  // (Note: with x/a substitution, a_old = 1/a_new, so λ_new = λ_old)
  const lambda = (a * Math.sqrt(2 * mass * V0)) / HBAR;

  // Maximum number of bound states
  const nMax = Math.floor(lambda - 0.5);
  const actualNumStates = Math.min(numStates, nMax + 1);

  if (actualNumStates <= 0) {
    throw new Error(
      "Pöschl-Teller potential too shallow to support bound states",
    );
  }

  // Calculate energies: E_n = -V_0 * [(λ - n - 1/2)/λ]²
  // Ground state (n=0) has lowest energy, excited states have higher energy
  const energies: number[] = [];
  for (let n = 0; n < actualNumStates; n++) {
    const term = lambda - n - 0.5;
    const energy = (-V0 * (term * term)) / (lambda * lambda);
    energies.push(energy);
  }

  // Generate grid
  const numPoints = gridConfig.numPoints;
  const xGrid: number[] = [];
  const dx = (gridConfig.xMax - gridConfig.xMin) / (numPoints - 1);
  for (let i = 0; i < numPoints; i++) {
    xGrid.push(gridConfig.xMin + i * dx);
  }

  // Calculate wavefunctions
  // ψ_n(x) = N_n * sech^(λ-n-1/2)(x/a) * P_n^(λ-n-1/2, λ-n-1/2)(tanh(x/a))
  // where P is the Jacobi polynomial
  const wavefunctions: number[][] = [];

  for (let n = 0; n < actualNumStates; n++) {
    const wavefunction: number[] = [];
    const alpha = lambda - n - 0.5;

    // Normalization (with 1/a factor from the variable change)
    const normalization =
      Math.sqrt(((1 / a) * (2 * alpha)) / factorial(n)) *
      Math.sqrt(factorial(n));

    for (const x of xGrid) {
      const tanhVal = Math.tanh(x / a);
      const sechVal = 1.0 / Math.cosh(x / a);

      // Use Legendre polynomials for Jacobi with α=β
      const jacobiPoly = jacobiPolynomial(n, alpha, alpha, tanhVal);
      const value = normalization * Math.pow(sechVal, alpha) * jacobiPoly;

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
