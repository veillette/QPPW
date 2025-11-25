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

/**
 * Create the potential function for an Eckart potential.
 * V(x) = V_0 / (1 + exp(x/a))² - V_1 / (1 + exp(x/a))
 *
 * @param potentialDepth - Potential depth V_0 in Joules
 * @param barrierHeight - Barrier height V_1 in Joules
 * @param wellWidth - Width parameter a in meters
 * @returns Potential function V(x) in Joules
 */
export function createEckartPotential(
  potentialDepth: number,
  barrierHeight: number,
  wellWidth: number,
): (x: number) => number {
  const V0 = potentialDepth;
  const V1 = barrierHeight;
  const a = wellWidth;

  return (x: number) => {
    const expVal = Math.exp(x / a);
    const onePlusExp = 1 + expVal;
    return V0 / (onePlusExp * onePlusExp) - V1 / onePlusExp;
  };
}

/**
 * Calculate classical probability density for an Eckart potential.
 * P(x) ∝ 1/v(x) = 1/√[2(E - V(x))/m]
 *
 * @param potentialDepth - Potential depth V_0 in Joules
 * @param barrierHeight - Barrier height V_1 in Joules
 * @param wellWidth - Width parameter a in meters
 * @param energy - Energy of the particle in Joules
 * @param mass - Particle mass in kg
 * @param xGrid - Array of x positions in meters
 * @returns Array of normalized classical probability density values (in 1/meters)
 */
export function calculateEckartPotentialClassicalProbability(
  potentialDepth: number,
  barrierHeight: number,
  wellWidth: number,
  energy: number,
  mass: number,
  xGrid: number[],
): number[] {
  const potentialFn = createEckartPotential(
    potentialDepth,
    barrierHeight,
    wellWidth,
  );

  const classicalProbability: number[] = [];
  let integralSum = 0;

  // Calculate unnormalized probability
  for (let i = 0; i < xGrid.length; i++) {
    const kineticEnergy = energy - potentialFn(xGrid[i]);

    if (kineticEnergy <= 0) {
      classicalProbability.push(0);
    } else {
      const epsilon = 1e-10 * Math.abs(potentialDepth);
      const probability =
        1 / Math.sqrt((2 * Math.max(kineticEnergy, epsilon)) / mass);
      classicalProbability.push(probability);

      if (i > 0) {
        const dx = xGrid[i] - xGrid[i - 1];
        integralSum += (probability + classicalProbability[i - 1]) * dx / 2;
      }
    }
  }

  // Normalize
  if (integralSum > 0) {
    for (let i = 0; i < classicalProbability.length; i++) {
      classicalProbability[i] /= integralSum;
    }
  }

  return classicalProbability;
}

/**
 * Calculate the classical turning points for an Eckart potential.
 * Solve E = V_0 / (1 + exp(x/a))² - V_1 / (1 + exp(x/a)) for x numerically
 *
 * @param potentialDepth - Potential depth V_0 in Joules
 * @param barrierHeight - Barrier height V_1 in Joules
 * @param wellWidth - Width parameter a in meters
 * @param energy - Energy of the particle in Joules
 * @returns Object with left and right turning point positions (in meters)
 */
export function calculateEckartPotentialTurningPoints(
  potentialDepth: number,
  barrierHeight: number,
  wellWidth: number,
  energy: number,
): { left: number; right: number } {
  const potentialFn = createEckartPotential(
    potentialDepth,
    barrierHeight,
    wellWidth,
  );
  const a = wellWidth;

  // Search for turning points using bisection
  const searchRange = 20 * a;

  // Find left turning point (negative x)
  let leftLow = -searchRange;
  let leftHigh = 0;

  for (let iter = 0; iter < 50; iter++) {
    const mid = (leftLow + leftHigh) / 2;
    const diff = potentialFn(mid) - energy;

    if (Math.abs(diff) < 1e-14 * Math.abs(energy)) {
      leftLow = mid;
      break;
    }

    if (diff > 0) {
      leftLow = mid;
    } else {
      leftHigh = mid;
    }
  }

  // Find right turning point (positive x)
  let rightLow = 0;
  let rightHigh = searchRange;

  for (let iter = 0; iter < 50; iter++) {
    const mid = (rightLow + rightHigh) / 2;
    const diff = potentialFn(mid) - energy;

    if (Math.abs(diff) < 1e-14 * Math.abs(energy)) {
      rightHigh = mid;
      break;
    }

    if (diff > 0) {
      rightHigh = mid;
    } else {
      rightLow = mid;
    }
  }

  return {
    left: (leftLow + leftHigh) / 2,
    right: (rightLow + rightHigh) / 2,
  };
}

/**
 * Calculate wavefunction zeros for Eckart potential (numerical approach).
 * Finds zeros by detecting sign changes in the wavefunction.
 *
 * @param potentialDepth - Potential depth V_0 in Joules
 * @param barrierHeight - Barrier height V_1 in Joules
 * @param wellWidth - Width parameter a in meters
 * @param mass - Particle mass in kg
 * @param stateIndex - Index of the eigenstate (0 for ground state, etc.)
 * @param searchRange - Range to search for zeros (in meters)
 * @returns Array of x positions (in meters) where wavefunction is zero
 */
export function calculateEckartPotentialWavefunctionZeros(
  potentialDepth: number,
  barrierHeight: number,
  wellWidth: number,
  mass: number,
  stateIndex: number,
  searchRange: number = 20e-9,
): number[] {
  const { HBAR } = QuantumConstants;
  const V0 = potentialDepth;
  const V1 = barrierHeight;
  const a = wellWidth;
  const n = stateIndex;

  const alpha_param = (a * Math.sqrt(2 * mass * V0)) / HBAR;
  const beta_param =
    (V1 * a * Math.sqrt(2 * mass)) / (2 * HBAR * Math.sqrt(V0));

  const s1 = -0.5 + Math.sqrt(0.25 + alpha_param);
  const s2 = -0.5 + Math.sqrt(0.25 + alpha_param - beta_param);

  const alpha_jac = 2 * (s2 - n) - 1;
  const beta_jac = 2 * (s1 - s2 + n);

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

  // Ground state has no zeros
  if (n === 0) {
    return [];
  }

  const zeros: number[] = [];
  const numSamples = 1000;
  const dx = (2 * searchRange) / numSamples;

  let prevX = -searchRange;
  let prevXi = Math.exp(prevX / a);
  let prevXiPlus1 = 1 + prevXi;
  let prevJacobiArg = 1 - (2 * prevXi) / prevXiPlus1;
  let prevJacobi = jacobiPolynomial(n, alpha_jac, beta_jac, prevJacobiArg);
  let prevVal =
    normalization *
    Math.pow(prevXi, s2 - n) *
    Math.pow(prevXiPlus1, -s1 - s2 + n) *
    prevJacobi;

  for (let i = 1; i <= numSamples; i++) {
    const x = -searchRange + i * dx;
    const xi = Math.exp(x / a);
    const xiPlus1 = 1 + xi;
    const jacobiArg = 1 - (2 * xi) / xiPlus1;
    const jacobiPoly = jacobiPolynomial(n, alpha_jac, beta_jac, jacobiArg);
    const val =
      normalization *
      Math.pow(xi, s2 - n) *
      Math.pow(xiPlus1, -s1 - s2 + n) *
      jacobiPoly;

    // Sign change detected
    if (prevVal * val < 0) {
      // Use bisection to refine
      let left = prevX;
      let right = x;
      for (let iter = 0; iter < 20; iter++) {
        const mid = (left + right) / 2;
        const midXi = Math.exp(mid / a);
        const midXiPlus1 = 1 + midXi;
        const midJacobiArg = 1 - (2 * midXi) / midXiPlus1;
        const midJacobi = jacobiPolynomial(
          n,
          alpha_jac,
          beta_jac,
          midJacobiArg,
        );
        const valMid =
          normalization *
          Math.pow(midXi, s2 - n) *
          Math.pow(midXiPlus1, -s1 - s2 + n) *
          midJacobi;

        if (Math.abs(valMid) < 1e-12) {
          zeros.push(mid);
          break;
        }

        if (valMid * prevVal < 0) {
          right = mid;
        } else {
          left = mid;
        }

        if (iter === 19) {
          zeros.push((left + right) / 2);
        }
      }
    }

    prevX = x;
    prevVal = val;
  }

  return zeros;
}

/**
 * Calculate the second derivative of the wavefunction for an Eckart potential.
 * Uses numerical differentiation on the analytical wavefunction.
 *
 * @param potentialDepth - Potential depth V_0 in Joules
 * @param barrierHeight - Barrier height V_1 in Joules
 * @param wellWidth - Width parameter a in meters
 * @param mass - Particle mass in kg
 * @param stateIndex - Index of the eigenstate (0 for ground state, etc.)
 * @param xGrid - Array of x positions in meters where derivatives should be evaluated
 * @returns Array of second derivative values
 */
export function calculateEckartPotentialWavefunctionSecondDerivative(
  potentialDepth: number,
  barrierHeight: number,
  wellWidth: number,
  mass: number,
  stateIndex: number,
  xGrid: number[],
): number[] {
  const { HBAR } = QuantumConstants;
  const V0 = potentialDepth;
  const V1 = barrierHeight;
  const a = wellWidth;
  const n = stateIndex;

  const alpha_param = (a * Math.sqrt(2 * mass * V0)) / HBAR;
  const beta_param =
    (V1 * a * Math.sqrt(2 * mass)) / (2 * HBAR * Math.sqrt(V0));

  const s1 = -0.5 + Math.sqrt(0.25 + alpha_param);
  const s2 = -0.5 + Math.sqrt(0.25 + alpha_param - beta_param);

  const alpha_jac = 2 * (s2 - n) - 1;
  const beta_jac = 2 * (s1 - s2 + n);

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

  const secondDerivative: number[] = [];
  const h = 1e-12; // Small step for numerical differentiation

  for (const x of xGrid) {
    // Evaluate at x-h, x, x+h
    const xMinus = x - h;
    const xPlus = x + h;

    const xiMinus = Math.exp(xMinus / a);
    const xiMinusPlus1 = 1 + xiMinus;
    const jacobiArgMinus = 1 - (2 * xiMinus) / xiMinusPlus1;
    const jacobiMinus = jacobiPolynomial(
      n,
      alpha_jac,
      beta_jac,
      jacobiArgMinus,
    );
    const psiMinus =
      normalization *
      Math.pow(xiMinus, s2 - n) *
      Math.pow(xiMinusPlus1, -s1 - s2 + n) *
      jacobiMinus;

    const xi = Math.exp(x / a);
    const xiPlus1 = 1 + xi;
    const jacobiArg = 1 - (2 * xi) / xiPlus1;
    const jacobi = jacobiPolynomial(n, alpha_jac, beta_jac, jacobiArg);
    const psi =
      normalization *
      Math.pow(xi, s2 - n) *
      Math.pow(xiPlus1, -s1 - s2 + n) *
      jacobi;

    const xiPlus = Math.exp(xPlus / a);
    const xiPlusPlus1 = 1 + xiPlus;
    const jacobiArgPlus = 1 - (2 * xiPlus) / xiPlusPlus1;
    const jacobiPlusVal = jacobiPolynomial(
      n,
      alpha_jac,
      beta_jac,
      jacobiArgPlus,
    );
    const psiPlus =
      normalization *
      Math.pow(xiPlus, s2 - n) *
      Math.pow(xiPlusPlus1, -s1 - s2 + n) *
      jacobiPlusVal;

    // Second derivative using central difference
    const secondDeriv = (psiPlus - 2 * psi + psiMinus) / (h * h);
    secondDerivative.push(secondDeriv);
  }

  return secondDerivative;
}
