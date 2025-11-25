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

/**
 * Create the potential function for a Pöschl-Teller potential.
 * V(x) = -V_0 / cosh²(x/a)
 *
 * @param potentialDepth - Potential depth V_0 in Joules (positive value)
 * @param wellWidth - Width parameter a in meters
 * @returns Potential function V(x) in Joules
 */
export function createPoschlTellerPotential(
  potentialDepth: number,
  wellWidth: number,
): (x: number) => number {
  const V0 = potentialDepth;
  const a = wellWidth;

  return (x: number) => {
    const coshVal = Math.cosh(x / a);
    return -V0 / (coshVal * coshVal);
  };
}

/**
 * Calculate classical probability density for a Pöschl-Teller potential.
 * P(x) ∝ 1/v(x) = 1/√[2(E - V(x))/m]
 *
 * @param potentialDepth - Potential depth V_0 in Joules (positive value)
 * @param wellWidth - Width parameter a in meters
 * @param energy - Energy of the particle in Joules
 * @param mass - Particle mass in kg
 * @param xGrid - Array of x positions in meters
 * @returns Array of normalized classical probability density values (in 1/meters)
 */
export function calculatePoschlTellerClassicalProbability(
  potentialDepth: number,
  wellWidth: number,
  energy: number,
  mass: number,
  xGrid: number[],
): number[] {
  const potentialFn = createPoschlTellerPotential(potentialDepth, wellWidth);

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
 * Calculate the classical turning points for a Pöschl-Teller potential.
 * Solve E = -V_0 / cosh²(x/a) for x
 *
 * @param potentialDepth - Potential depth V_0 in Joules (positive value)
 * @param wellWidth - Width parameter a in meters
 * @param energy - Energy of the particle in Joules
 * @returns Object with left and right turning point positions (in meters)
 */
export function calculatePoschlTellerTurningPoints(
  potentialDepth: number,
  wellWidth: number,
  energy: number,
): { left: number; right: number } {
  const V0 = potentialDepth;
  const a = wellWidth;

  // Solve: E = -V0 / cosh²(x/a)
  // => cosh²(x/a) = -V0 / E
  // => cosh(x/a) = sqrt(-V0 / E)
  // => x/a = ±acosh(sqrt(-V0 / E))
  // => x = ±a * acosh(sqrt(-V0 / E))

  const ratio = Math.sqrt(-V0 / energy);
  const turning = a * Math.acosh(ratio);

  return {
    left: -turning,
    right: turning,
  };
}

/**
 * Calculate wavefunction zeros for Pöschl-Teller potential (numerical approach).
 * Finds zeros by detecting sign changes in the wavefunction.
 *
 * @param potentialDepth - Potential depth V_0 in Joules (positive value)
 * @param wellWidth - Width parameter a in meters
 * @param mass - Particle mass in kg
 * @param stateIndex - Index of the eigenstate (0 for ground state, etc.)
 * @param searchRange - Range to search for zeros (in meters)
 * @returns Array of x positions (in meters) where wavefunction is zero
 */
export function calculatePoschlTellerWavefunctionZeros(
  potentialDepth: number,
  wellWidth: number,
  mass: number,
  stateIndex: number,
  searchRange: number = 20e-9,
): number[] {
  const { HBAR } = QuantumConstants;
  const V0 = potentialDepth;
  const a = wellWidth;
  const n = stateIndex;

  const lambda = (a * Math.sqrt(2 * mass * V0)) / HBAR;
  const alpha = lambda - n - 0.5;
  const normalization =
    Math.sqrt(((1 / a) * (2 * alpha)) / factorial(n)) * Math.sqrt(factorial(n));

  // Ground state has no zeros
  if (n === 0) {
    return [];
  }

  const zeros: number[] = [];
  const numSamples = 1000;
  const dx = (2 * searchRange) / numSamples;

  let prevX = -searchRange;
  let prevTanh = Math.tanh(prevX / a);
  let prevSech = 1.0 / Math.cosh(prevX / a);
  let prevJacobi = jacobiPolynomial(n, alpha, alpha, prevTanh);
  let prevVal = normalization * Math.pow(prevSech, alpha) * prevJacobi;

  for (let i = 1; i <= numSamples; i++) {
    const x = -searchRange + i * dx;
    const tanhVal = Math.tanh(x / a);
    const sechVal = 1.0 / Math.cosh(x / a);
    const jacobiPoly = jacobiPolynomial(n, alpha, alpha, tanhVal);
    const val = normalization * Math.pow(sechVal, alpha) * jacobiPoly;

    // Sign change detected
    if (prevVal * val < 0) {
      // Use bisection to refine
      let left = prevX;
      let right = x;
      for (let iter = 0; iter < 20; iter++) {
        const mid = (left + right) / 2;
        const midTanh = Math.tanh(mid / a);
        const midSech = 1.0 / Math.cosh(mid / a);
        const midJacobi = jacobiPolynomial(n, alpha, alpha, midTanh);
        const valMid = normalization * Math.pow(midSech, alpha) * midJacobi;

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
 * Calculate the second derivative of the wavefunction for a Pöschl-Teller potential.
 * Uses numerical differentiation on the analytical wavefunction.
 *
 * @param potentialDepth - Potential depth V_0 in Joules (positive value)
 * @param wellWidth - Width parameter a in meters
 * @param mass - Particle mass in kg
 * @param stateIndex - Index of the eigenstate (0 for ground state, etc.)
 * @param xGrid - Array of x positions in meters where derivatives should be evaluated
 * @returns Array of second derivative values
 */
export function calculatePoschlTellerWavefunctionSecondDerivative(
  potentialDepth: number,
  wellWidth: number,
  mass: number,
  stateIndex: number,
  xGrid: number[],
): number[] {
  const { HBAR } = QuantumConstants;
  const V0 = potentialDepth;
  const a = wellWidth;
  const n = stateIndex;

  const lambda = (a * Math.sqrt(2 * mass * V0)) / HBAR;
  const alpha = lambda - n - 0.5;
  const normalization =
    Math.sqrt(((1 / a) * (2 * alpha)) / factorial(n)) * Math.sqrt(factorial(n));

  const secondDerivative: number[] = [];
  const h = 1e-12; // Small step for numerical differentiation

  for (const x of xGrid) {
    // Evaluate at x-h, x, x+h
    const xMinus = x - h;
    const xPlus = x + h;

    const tanhMinus = Math.tanh(xMinus / a);
    const sechMinus = 1.0 / Math.cosh(xMinus / a);
    const jacobiMinus = jacobiPolynomial(n, alpha, alpha, tanhMinus);
    const psiMinus = normalization * Math.pow(sechMinus, alpha) * jacobiMinus;

    const tanh = Math.tanh(x / a);
    const sech = 1.0 / Math.cosh(x / a);
    const jacobi = jacobiPolynomial(n, alpha, alpha, tanh);
    const psi = normalization * Math.pow(sech, alpha) * jacobi;

    const tanhPlus = Math.tanh(xPlus / a);
    const sechPlus = 1.0 / Math.cosh(xPlus / a);
    const jacobiPlus = jacobiPolynomial(n, alpha, alpha, tanhPlus);
    const psiPlus = normalization * Math.pow(sechPlus, alpha) * jacobiPlus;

    // Second derivative using central difference
    const secondDeriv = (psiPlus - 2 * psi + psiMinus) / (h * h);
    secondDerivative.push(secondDeriv);
  }

  return secondDerivative;
}
