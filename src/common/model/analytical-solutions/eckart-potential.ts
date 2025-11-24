/**
 * Analytical solution for the Eckart potential.
 * V(x) = V_0 / (1 + exp(x/a))² - V_1 / (1 + exp(x/a))
 *
 * This potential is useful for modeling molecular barriers and chemical reactions.
 *
 * REFERENCES:
 * - Eckart, C. (1930). "The Penetration of a Potential Barrier by Electrons"
 *   Physical Review, 35(11), 1303-1309.
 *   https://doi.org/10.1103/PhysRev.35.1303
 *   ORIGINAL PAPER: Introduced this potential for modeling barrier penetration.
 *
 * - Flügge, S. (1999). "Practical Quantum Mechanics". Springer.
 *   Problem 40, pp. 97-99. https://doi.org/10.1007/978-3-642-61995-3
 *   Detailed solution using hypergeometric functions.
 *
 * - Cooper, F., Khare, A., & Sukhatme, U. (1995). "Supersymmetry and quantum mechanics"
 *   Physics Reports, 251(5-6), 267-385.
 *   https://doi.org/10.1016/0370-1573(94)00080-M
 *   Section 3.3, pp. 283-285: Eckart as a shape-invariant potential.
 *
 * - Natanzon, G. A. (1979). "General properties of potentials for which the Schrödinger equation
 *   can be solved by means of hypergeometric functions". Theoretical and Mathematical Physics, 38(2), 146-153.
 *   https://doi.org/10.1007/BF01016836
 *   Classification including the Eckart potential.
 *
 * - Wei, H., & Xia, X. (1991). "Algebraic approach to energy spectra of quantum systems"
 *   Journal of Physics A, 24(1), 151.
 *   https://doi.org/10.1088/0305-4470/24/1/024
 *   Algebraic solution methods.
 *
 * ENERGY EIGENVALUES:
 *   E_n = -(ℏ²/2ma²)(s_2 - n)²,  n = 0, 1, 2, ..., n_max
 *   where s_1 = -1/2 + √(1/4 + α), s_2 = -1/2 + √(1/4 + α - β)
 *   α = a²(2mV_0)/ℏ², β = V_1·a·√(2m)/(2ℏ√V_0)
 *
 * WAVEFUNCTIONS:
 *   ψ_n(ξ) = N_n · ξ^(s_2-n) · (1+ξ)^(-s_1-s_2+n) · P_n^(α,β)(1-2ξ/(1+ξ))
 *   where ξ = exp(x/a) and P_n^(α,β) are Jacobi polynomials
 */

import QuantumConstants from "../QuantumConstants.js";
import { BoundStateResult, GridConfig } from "../PotentialFunction.js";
import { jacobiPolynomial, factorial, logGamma } from "./math-utilities.js";

/**
 * Analytical solution for the Eckart potential.
 * V(x) = V_0 / (1 + exp(x/a))² - V_1 / (1 + exp(x/a))
 *
 * This potential is useful for modeling molecular barriers and chemical reactions.
 *
 * @param potentialDepth - Potential depth V_0 in Joules
 * @param barrierHeight - Barrier height V_1 in Joules
 * @param wellWidth - Width parameter a in meters
 * @param mass - Particle mass in kg
 * @param numStates - Number of energy levels to calculate
 * @param gridConfig - Grid configuration for wavefunction evaluation
 * @returns Bound state results with exact energies and wavefunctions
 */
export function solveEckartPotential(
  potentialDepth: number,
  barrierHeight: number,
  wellWidth: number,
  mass: number,
  numStates: number,
  gridConfig: GridConfig,
): BoundStateResult {
  const { HBAR } = QuantumConstants;
  const V0 = potentialDepth;
  const V1 = barrierHeight;
  const a = wellWidth;

  // Calculate dimensionless parameters (with x/a substitution)
  // α = a * sqrt(2*m*V_0) / ℏ
  // β = V_1 * a * sqrt(2*m) / (2*ℏ*sqrt(V_0))
  const alpha_param = (a * Math.sqrt(2 * mass * V0)) / HBAR;
  const beta_param =
    (V1 * a * Math.sqrt(2 * mass)) / (2 * HBAR * Math.sqrt(V0));

  // For bound states, we need certain conditions on α and β
  const s1 = -0.5 + Math.sqrt(0.25 + alpha_param);
  const s2 = -0.5 + Math.sqrt(0.25 + alpha_param - beta_param);

  if (s2 <= 0) {
    throw new Error("Eckart potential too shallow to support bound states");
  }

  // Maximum number of bound states
  const nMax = Math.floor(Math.min(s1, s2));
  const actualNumStates = Math.min(numStates, nMax);

  if (actualNumStates <= 0) {
    throw new Error("Eckart potential too shallow to support bound states");
  }

  // Calculate energies
  // E_n = -(ℏ²/2ma²) * (s_2 - n)²
  const energies: number[] = [];
  const energyFactor = (HBAR * HBAR) / (2 * mass * a * a);

  for (let n = 0; n < actualNumStates; n++) {
    const energy = -energyFactor * Math.pow(s2 - n, 2);
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
  // ψ_n(ξ) = N_n * ξ^(s_2-n) * (1+ξ)^(-s_1-s_2+n) * P_n^(α,β)(1-2ξ/(1+ξ))
  // where ξ = exp(x/a)
  const wavefunctions: number[][] = [];

  for (let n = 0; n < actualNumStates; n++) {
    const wavefunction: number[] = [];

    const alpha_jac = 2 * (s2 - n) - 1;
    const beta_jac = 2 * (s1 - s2 + n);

    // Normalization (with 1/a factor from variable change)
    const normalization = Math.sqrt(
      (1 / a) *
        factorial(n) *
        Math.exp(
          logGamma(2 * s2 - n) +
            logGamma(alpha_jac + beta_jac + n + 1) -
            logGamma(n + alpha_jac + 1) -
            logGamma(n + beta_jac + 1),
        ),
    );

    for (const x of xGrid) {
      const xi = Math.exp(x / a);
      const xiPlus1 = 1 + xi;

      // Jacobi polynomial argument
      const jacobiArg = 1 - (2 * xi) / xiPlus1;
      const jacobiPoly = jacobiPolynomial(n, alpha_jac, beta_jac, jacobiArg);

      // Calculate wavefunction
      const value =
        normalization *
        Math.pow(xi, s2 - n) *
        Math.pow(xiPlus1, -s1 - s2 + n) *
        jacobiPoly;

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
