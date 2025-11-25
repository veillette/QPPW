/**
 * Analytical solution for a quantum harmonic oscillator.
 * V(x) = (1/2) * k * x^2 = (1/2) * m * ω^2 * x^2
 *
 * REFERENCES:
 * - Griffiths, D. J., & Schroeter, D. F. (2018). "Introduction to Quantum Mechanics" (3rd ed.).
 *   Cambridge University Press. Section 2.3, pp. 40-55.
 *   https://doi.org/10.1017/9781316995433
 *   Complete derivation using operator methods and Hermite polynomials.
 *
 * - Shankar, R. (1994). "Principles of Quantum Mechanics" (2nd ed.). Springer.
 *   Section 7.3, pp. 164-182. https://doi.org/10.1007/978-1-4757-0576-8
 *   Algebraic method using ladder operators.
 *
 * - Dirac, P. A. M. (1927). "The Quantum Theory of the Emission and Absorption of Radiation".
 *   Proceedings of the Royal Society A, 114(767), 243-265.
 *   https://doi.org/10.1098/rspa.1927.0039
 *   Introduction of creation and annihilation operators.
 *
 * - Abramowitz, M., & Stegun, I. A. (1964). "Handbook of Mathematical Functions".
 *   National Bureau of Standards. Section 22, pp. 773-802.
 *   https://doi.org/10.1119/1.15378
 *   Properties of Hermite polynomials used in wavefunctions.
 *
 * ENERGY EIGENVALUES:
 *   E_n = ℏω(n + 1/2),  n = 0, 1, 2, ...
 *   where ω = √(k/m)
 *
 * WAVEFUNCTIONS:
 *   ψ_n(x) = (1/√(2^n n!)) · (mω/πℏ)^(1/4) · exp(-mωx²/(2ℏ)) · H_n(√(mω/ℏ) x)
 *   where H_n are the Hermite polynomials
 */

import QuantumConstants from "../QuantumConstants.js";
import {
  BoundStateResult,
  GridConfig,
  PotentialFunction,
} from "../PotentialFunction.js";
import { hermitePolynomial, factorial } from "./math-utilities.js";

/**
 * Create the potential function for a harmonic oscillator.
 * V(x) = (1/2) * k * x^2
 *
 * @param springConstant - Spring constant k in N/m
 * @returns Potential function V(x) in Joules
 */
export function createHarmonicOscillatorPotential(
  springConstant: number,
): PotentialFunction {
  return (x: number) => {
    return 0.5 * springConstant * x * x;
  };
}

/**
 * Calculate classical probability density for a harmonic oscillator.
 * For a classical harmonic oscillator, the probability density is:
 * P(x) = 1 / (π * √(A² - x²))  for |x| < A, where A = √(2E/k) is the amplitude
 *
 * This is one of the cases where the renormalization can be computed analytically.
 * The integral ∫_{-A}^{A} 1/√(A² - x²) dx = π, so the normalization is 1/π.
 *
 * @param springConstant - Spring constant k in N/m
 * @param energy - Energy of the particle in Joules
 * @param mass - Particle mass in kg (unused for harmonic oscillator)
 * @param xGrid - Array of x positions in meters
 * @returns Array of normalized classical probability density values (in 1/meters)
 */
export function calculateHarmonicOscillatorClassicalProbability(
  springConstant: number,
  energy: number,
  _mass: number,
  xGrid: number[],
): number[] {
  const probability: number[] = [];

  // Classical amplitude: A = √(2E/k)
  const amplitude = Math.sqrt((2 * energy) / springConstant);

  // Classical probability density: P(x) = 1 / (π * √(A² - x²))
  // This is already analytically normalized

  for (const x of xGrid) {
    if (Math.abs(x) < amplitude) {
      const arg = amplitude * amplitude - x * x;
      if (arg > 0) {
        const prob = 1 / (Math.PI * Math.sqrt(arg));
        probability.push(prob);
      } else {
        probability.push(0);
      }
    } else {
      probability.push(0);
    }
  }

  return probability;
}

/**
 * Calculate the positions of wavefunction zeros (nodes) for a harmonic oscillator.
 * The zeros are the roots of the Hermite polynomial H_n(ξ) where ξ = √(mω/ℏ) * x
 *
 * Uses numerical root finding to locate all n zeros of the nth Hermite polynomial.
 *
 * @param springConstant - Spring constant k in N/m
 * @param mass - Particle mass in kg
 * @param stateIndex - Index of the eigenstate (0 for ground state, 1 for first excited, etc.)
 * @returns Array of x positions (in meters) where the wavefunction is zero
 */
export function calculateHarmonicOscillatorWavefunctionZeros(
  springConstant: number,
  mass: number,
  stateIndex: number,
): number[] {
  const { HBAR } = QuantumConstants;
  const n = stateIndex; // Quantum number (0, 1, 2, ...)
  const omega = Math.sqrt(springConstant / mass);
  const alpha = Math.sqrt((mass * omega) / HBAR);

  // Ground state has no zeros
  if (n === 0) {
    return [];
  }

  const zeros: number[] = [];

  // Hermite polynomials H_n have n real, distinct zeros
  // They are symmetric about the origin, and for n odd, 0 is a zero
  // Use bisection to find zeros in intervals where H_n changes sign

  // For H_n, all zeros lie in the interval [-√(2n+1), √(2n+1)]
  const searchRange = Math.sqrt(2 * n + 1);
  const numIntervals = 100 * n; // Use many intervals for accuracy
  const dxi = (2 * searchRange) / numIntervals;

  // Search for sign changes
  let prevXi = -searchRange;
  let prevVal = hermitePolynomial(n, prevXi);

  for (let i = 1; i <= numIntervals; i++) {
    const xi = -searchRange + i * dxi;
    const val = hermitePolynomial(n, xi);

    // Sign change detected - use bisection to find zero
    if (prevVal * val < 0) {
      let left = prevXi;
      let right = xi;

      // Bisection method
      for (let iter = 0; iter < 50; iter++) {
        const mid = (left + right) / 2;
        const midVal = hermitePolynomial(n, mid);

        if (Math.abs(midVal) < 1e-12) {
          // Found zero
          zeros.push(mid / alpha); // Convert from ξ to x
          break;
        }

        if (midVal * hermitePolynomial(n, left) < 0) {
          right = mid;
        } else {
          left = mid;
        }

        if (iter === 49) {
          // Use midpoint as approximation
          zeros.push((left + right) / (2 * alpha));
        }
      }
    }

    prevXi = xi;
    prevVal = val;
  }

  return zeros;
}

/**
 * Calculate the classical turning points for a harmonic oscillator.
 * Turning points occur where E = V(x) = (1/2)kx², so x = ±√(2E/k)
 *
 * @param springConstant - Spring constant k in N/m
 * @param energy - Energy of the particle in Joules
 * @returns Object with left and right turning point positions (in meters)
 */
export function calculateHarmonicOscillatorTurningPoints(
  springConstant: number,
  energy: number,
): { left: number; right: number } {
  // E = (1/2)kx² => x = ±√(2E/k)
  const amplitude = Math.sqrt((2 * energy) / springConstant);

  return {
    left: -amplitude,
    right: amplitude,
  };
}

/**
 * Calculate the first and second derivatives of the wavefunction for a harmonic oscillator.
 *
 * For ψ_n(x) = N_n * exp(-αx²/2) * H_n(√α x) where α = mω/ℏ and N_n is normalization:
 * - ψ'_n uses the chain rule and Hermite polynomial recursion: H'_n(ξ) = 2n H_{n-1}(ξ)
 * - ψ''_n uses further differentiation
 *
 * @param springConstant - Spring constant k in N/m
 * @param mass - Particle mass in kg
 * @param stateIndex - Index of the eigenstate (0 for ground state, 1 for first excited, etc.)
 * @param xGrid - Array of x positions in meters where derivatives should be evaluated
 * @returns Object with first and second derivative arrays
 */
export function calculateHarmonicOscillatorWavefunctionDerivatives(
  springConstant: number,
  mass: number,
  stateIndex: number,
  xGrid: number[],
): { firstDerivative: number[]; secondDerivative: number[] } {
  const { HBAR } = QuantumConstants;
  const n = stateIndex; // Quantum number (0, 1, 2, ...)
  const omega = Math.sqrt(springConstant / mass);
  const alpha = Math.sqrt((mass * omega) / HBAR);
  const normalization =
    (1 / Math.sqrt(Math.pow(2, n) * factorial(n))) *
    Math.pow(alpha / Math.PI, 0.25);

  const firstDerivative: number[] = [];
  const secondDerivative: number[] = [];

  for (const x of xGrid) {
    const xi = alpha * x;
    const gaussianFactor = Math.exp((-xi * xi) / 2);
    const hermite_n = hermitePolynomial(n, xi);

    // For first derivative, use: ψ' = N * [-α*x*exp(-αx²/2)*H_n + exp(-αx²/2)*H'_n * α]
    // where H'_n(ξ) = 2n*H_{n-1}(ξ)
    let hermite_derivative = 0;
    if (n > 0) {
      const hermite_n_minus_1 = hermitePolynomial(n - 1, xi);
      hermite_derivative = 2 * n * hermite_n_minus_1;
    }

    const firstDeriv =
      normalization *
      gaussianFactor *
      (-xi * hermite_n + hermite_derivative);
    firstDerivative.push(firstDeriv);

    // For second derivative, differentiate again
    // ψ'' = N * [α²*x²*exp(-αx²/2)*H_n - α*exp(-αx²/2)*H_n
    //            - 2α*x*exp(-αx²/2)*H'_n + α²*exp(-αx²/2)*H''_n]
    // where H''_n(ξ) = 2n*H'_{n-1}(ξ) = 4n(n-1)*H_{n-2}(ξ)

    let hermite_second_derivative = 0;
    if (n > 1) {
      const hermite_n_minus_2 = hermitePolynomial(n - 2, xi);
      hermite_second_derivative = 4 * n * (n - 1) * hermite_n_minus_2;
    }

    const secondDeriv =
      normalization *
      gaussianFactor *
      (xi * xi * hermite_n -
        hermite_n -
        2 * xi * hermite_derivative +
        hermite_second_derivative);
    secondDerivative.push(secondDeriv);
  }

  return { firstDerivative, secondDerivative };
}

/**
 * Analytical solution for a quantum harmonic oscillator.
 * V(x) = (1/2) * k * x^2 = (1/2) * m * ω^2 * x^2
 *
 * @param springConstant - Spring constant k in N/m
 * @param mass - Particle mass in kg
 * @param numStates - Number of energy levels to calculate
 * @param gridConfig - Grid configuration for wavefunction evaluation
 * @returns Bound state results with exact energies and wavefunctions
 */
export function solveHarmonicOscillator(
  springConstant: number,
  mass: number,
  numStates: number,
  gridConfig: GridConfig,
): BoundStateResult {
  const { HBAR } = QuantumConstants;
  const omega = Math.sqrt(springConstant / mass);

  // Calculate energies: E_n = ℏω(n + 1/2) for n = 0, 1, 2, ...
  const energies: number[] = [];
  for (let n = 0; n < numStates; n++) {
    const energy = HBAR * omega * (n + 0.5);
    energies.push(energy);
  }

  // Generate grid
  const numPoints = gridConfig.numPoints;
  const xGrid: number[] = [];
  const dx = (gridConfig.xMax - gridConfig.xMin) / (numPoints - 1);
  for (let i = 0; i < numPoints; i++) {
    xGrid.push(gridConfig.xMin + i * dx);
  }

  // Calculate wavefunctions using Hermite polynomials
  // ψ_n(x) = (1/√(2^n n!)) * (mω/πℏ)^(1/4) * exp(-mωx^2/(2ℏ)) * H_n(√(mω/ℏ) x)
  const wavefunctions: number[][] = [];
  const alpha = Math.sqrt((mass * omega) / HBAR);

  for (let n = 0; n < numStates; n++) {
    const wavefunction: number[] = [];
    const normalization =
      (1 / Math.sqrt(Math.pow(2, n) * factorial(n))) *
      Math.pow(alpha / Math.PI, 0.25);

    for (const x of xGrid) {
      const xi = alpha * x;
      const hermite = hermitePolynomial(n, xi);
      const value = normalization * Math.exp((-xi * xi) / 2) * hermite;
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
 * Calculate coherent state coefficients for harmonic oscillator.
 * Coherent states are eigenstates of the annihilation operator and represent
 * the most classical-like quantum states (Glauber states).
 *
 * For a real displacement x₀, the coherent state is:
 * |α⟩ = e^(-α²/2) Σ (α^n/√n!) |n⟩
 * where α = √(mω/2ℏ) * x₀
 *
 * All coefficients are positive for real α (pure displacement, no momentum).
 *
 * @param displacement - Displacement from equilibrium in meters
 * @param springConstant - Spring constant k in N/m
 * @param mass - Particle mass in kg
 * @param numStates - Number of energy eigenstates to include
 * @returns Object with amplitudes and phases arrays
 */
export function calculateCoherentStateCoefficients(
  displacement: number,
  springConstant: number,
  mass: number,
  numStates: number,
): { amplitudes: number[]; phases: number[] } {
  const { HBAR } = QuantumConstants;
  const omega = Math.sqrt(springConstant / mass);

  // Calculate α = √(mω/2ℏ) * x₀
  const alpha = Math.sqrt((mass * omega) / (2 * HBAR)) * displacement;

  // Calculate coefficients: c_n = e^(-α²/2) * α^n / √n!
  const prefactor = Math.exp((-alpha * alpha) / 2);
  const amplitudes: number[] = [];

  for (let n = 0; n < numStates; n++) {
    const coeff = (prefactor * Math.pow(alpha, n)) / Math.sqrt(factorial(n));
    amplitudes.push(coeff);
  }

  // Normalize the coefficients (should already be normalized, but ensure it)
  const norm = Math.sqrt(amplitudes.reduce((sum, a) => sum + a * a, 0));
  for (let i = 0; i < amplitudes.length; i++) {
    amplitudes[i] /= norm;
  }

  // All phases are zero for a real coherent state (pure displacement)
  const phases = new Array(numStates).fill(0);

  return { amplitudes, phases };
}
