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
import { AnalyticalSolution } from "./AnalyticalSolution.js";

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
 * Calculate the first derivative of the wavefunction for a harmonic oscillator.
 *
 * For ψ_n(x) = N_n * exp(-αx²/2) * H_n(ξ) where ξ = αx and α = mω/ℏ:
 * ψ'_n(x) = N_n * α * exp(-αx²/2) * [H'_n(ξ) - ξ * H_n(ξ)]
 *         = N_n * α * exp(-αx²/2) * [2n H_{n-1}(ξ) - ξ * H_n(ξ)]
 * using the Hermite polynomial recursion relation H'_n(ξ) = 2n H_{n-1}(ξ)
 *
 * @param springConstant - Spring constant k in N/m
 * @param mass - Particle mass in kg
 * @param stateIndex - Index of the eigenstate (0 for ground state, 1 for first excited, etc.)
 * @param xGrid - Array of x positions in meters where derivatives should be evaluated
 * @returns Array of first derivative values
 */
export function calculateHarmonicOscillatorWavefunctionFirstDerivative(
  springConstant: number,
  mass: number,
  stateIndex: number,
  xGrid: number[],
): number[] {
  const { HBAR } = QuantumConstants;
  const n = stateIndex; // Quantum number (0, 1, 2, ...)
  const omega = Math.sqrt(springConstant / mass);
  const alpha = Math.sqrt((mass * omega) / HBAR);
  const normalization =
    (1 / Math.sqrt(Math.pow(2, n) * factorial(n))) *
    Math.pow(alpha / Math.PI, 0.25);

  const firstDerivative: number[] = [];

  for (const x of xGrid) {
    const xi = alpha * x;
    const gaussianFactor = Math.exp((-xi * xi) / 2);
    const hermite_n = hermitePolynomial(n, xi);

    // Calculate H'_n(ξ) = 2n H_{n-1}(ξ)
    let hermite_derivative = 0;
    if (n > 0) {
      const hermite_n_minus_1 = hermitePolynomial(n - 1, xi);
      hermite_derivative = 2 * n * hermite_n_minus_1;
    }

    // ψ' = N * α * exp(-αx²/2) * [H'_n(ξ) - ξ * H_n(ξ)]
    const firstDeriv =
      normalization *
      alpha *
      gaussianFactor *
      (hermite_derivative - xi * hermite_n);
    firstDerivative.push(firstDeriv);
  }

  return firstDerivative;
}

/**
 * Calculate the second derivative of the wavefunction for a harmonic oscillator.
 *
 * For ψ_n(x) = N_n * exp(-αx²/2) * H_n(√α x) where α = mω/ℏ and N_n is normalization:
 * ψ''_n uses the chain rule and Hermite polynomial recursion relations:
 * - H'_n(ξ) = 2n H_{n-1}(ξ)
 * - H''_n(ξ) = 4n(n-1) H_{n-2}(ξ)
 *
 * @param springConstant - Spring constant k in N/m
 * @param mass - Particle mass in kg
 * @param stateIndex - Index of the eigenstate (0 for ground state, 1 for first excited, etc.)
 * @param xGrid - Array of x positions in meters where derivatives should be evaluated
 * @returns Array of second derivative values
 */
export function calculateHarmonicOscillatorWavefunctionSecondDerivative(
  springConstant: number,
  mass: number,
  stateIndex: number,
  xGrid: number[],
): number[] {
  const { HBAR } = QuantumConstants;
  const n = stateIndex; // Quantum number (0, 1, 2, ...)
  const omega = Math.sqrt(springConstant / mass);
  const alpha = Math.sqrt((mass * omega) / HBAR);
  const normalization =
    (1 / Math.sqrt(Math.pow(2, n) * factorial(n))) *
    Math.pow(alpha / Math.PI, 0.25);

  const secondDerivative: number[] = [];

  for (const x of xGrid) {
    const xi = alpha * x;
    const gaussianFactor = Math.exp((-xi * xi) / 2);
    const hermite_n = hermitePolynomial(n, xi);

    // Calculate H'_n(ξ) = 2n H_{n-1}(ξ)
    let hermite_derivative = 0;
    if (n > 0) {
      const hermite_n_minus_1 = hermitePolynomial(n - 1, xi);
      hermite_derivative = 2 * n * hermite_n_minus_1;
    }

    // Calculate H''_n(ξ) = 4n(n-1) H_{n-2}(ξ)
    let hermite_second_derivative = 0;
    if (n > 1) {
      const hermite_n_minus_2 = hermitePolynomial(n - 2, xi);
      hermite_second_derivative = 4 * n * (n - 1) * hermite_n_minus_2;
    }

    // ψ'' = N * [α²x²*exp(-αx²/2)*H_n - α*exp(-αx²/2)*H_n
    //            - 2αx*exp(-αx²/2)*H'_n + α²*exp(-αx²/2)*H''_n]
    const secondDeriv =
      normalization *
      gaussianFactor *
      (xi * xi * hermite_n -
        hermite_n -
        2 * xi * hermite_derivative +
        hermite_second_derivative);
    secondDerivative.push(secondDeriv);
  }

  return secondDerivative;
}

/**
 * Calculate the minimum and maximum values of the wavefunction for a harmonic oscillator.
 *
 * For ψ_n(x) = N_n · exp(-αx²/2) · H_n(√α x), the function is sampled at multiple points
 * in the range [xMin, xMax] to find the extrema.
 *
 * @param springConstant - Spring constant k in N/m
 * @param mass - Particle mass in kg
 * @param stateIndex - Index of the eigenstate (0 for ground state, 1 for first excited, etc.)
 * @param xMin - Left boundary of the region in meters
 * @param xMax - Right boundary of the region in meters
 * @param numPoints - Number of points to sample (default: 1000)
 * @returns Object containing min and max values of the wavefunction
 */
export function calculateHarmonicOscillatorWavefunctionMinMax(
  springConstant: number,
  mass: number,
  stateIndex: number,
  xMin: number,
  xMax: number,
  numPoints: number = 1000,
): { min: number; max: number } {
  const { HBAR } = QuantumConstants;
  const n = stateIndex; // Quantum number (0, 1, 2, ...)
  const omega = Math.sqrt(springConstant / mass);
  const alpha = Math.sqrt((mass * omega) / HBAR);
  const normalization =
    (1 / Math.sqrt(Math.pow(2, n) * factorial(n))) *
    Math.pow(alpha / Math.PI, 0.25);

  let min = Infinity;
  let max = -Infinity;

  const dx = (xMax - xMin) / (numPoints - 1);

  for (let i = 0; i < numPoints; i++) {
    const x = xMin + i * dx;

    const xi = alpha * x;
    const gaussianFactor = Math.exp((-xi * xi) / 2);
    const hermite = hermitePolynomial(n, xi);
    const psi = normalization * gaussianFactor * hermite;

    if (psi < min) min = psi;
    if (psi > max) max = psi;
  }

  return { min, max };
}

/**
 * Calculate the minimum and maximum values of a superposition of wavefunctions
 * for a harmonic oscillator.
 *
 * The superposition is: Ψ(x,t) = Σ cₙ ψₙ(x) exp(-iEₙt/ℏ)
 * We return the min/max of the real part of this complex-valued function.
 *
 * @param springConstant - Spring constant k in N/m
 * @param mass - Particle mass in kg
 * @param coefficients - Complex coefficients for each eigenstate (as [real, imag] pairs)
 * @param energies - Energy eigenvalues in Joules
 * @param time - Time in seconds
 * @param xMin - Left boundary of the region in meters
 * @param xMax - Right boundary of the region in meters
 * @param numPoints - Number of points to sample (default: 1000)
 * @returns Object containing min and max values of the superposition's real part
 */
export function calculateHarmonicOscillatorSuperpositionMinMax(
  springConstant: number,
  mass: number,
  coefficients: Array<[number, number]>,
  energies: number[],
  time: number,
  xMin: number,
  xMax: number,
  numPoints: number = 1000,
): { min: number; max: number } {
  const { HBAR } = QuantumConstants;
  const omega = Math.sqrt(springConstant / mass);
  const alpha = Math.sqrt((mass * omega) / HBAR);

  let min = Infinity;
  let max = -Infinity;

  const dx = (xMax - xMin) / (numPoints - 1);

  for (let i = 0; i < numPoints; i++) {
    const x = xMin + i * dx;

    // Calculate superposition at this position
    let realPart = 0;

    for (let n = 0; n < coefficients.length; n++) {
      const [cReal, cImag] = coefficients[n];
      const energy = energies[n];

      // Calculate wavefunction value
      const normalization =
        (1 / Math.sqrt(Math.pow(2, n) * factorial(n))) *
        Math.pow(alpha / Math.PI, 0.25);
      const xi = alpha * x;
      const gaussianFactor = Math.exp((-xi * xi) / 2);
      const hermite = hermitePolynomial(n, xi);
      const psi = normalization * gaussianFactor * hermite;

      // Time evolution: exp(-iEt/ℏ) = cos(Et/ℏ) - i*sin(Et/ℏ)
      const phase = (-energy * time) / HBAR;
      const cosPhase = Math.cos(phase);
      const sinPhase = Math.sin(phase);

      // Complex multiplication: (cReal + i*cImag) * psi * (cosPhase - i*sinPhase)
      // Real part: cReal * psi * cosPhase + cImag * psi * sinPhase
      realPart += cReal * psi * cosPhase + cImag * psi * sinPhase;
    }

    if (realPart < min) min = realPart;
    if (realPart > max) max = realPart;
  }

  return { min, max };
}

/**
 * Class-based implementation of harmonic oscillator analytical solution.
 * Extends the AnalyticalSolution abstract base class.
 */
export class HarmonicOscillatorSolution extends AnalyticalSolution {
  constructor(
    private springConstant: number,
    private mass: number,
  ) {
    super();
  }

  solve(numStates: number, gridConfig: GridConfig): BoundStateResult {
    return solveHarmonicOscillator(
      this.springConstant,
      this.mass,
      numStates,
      gridConfig,
    );
  }

  createPotential(): PotentialFunction {
    return createHarmonicOscillatorPotential(this.springConstant);
  }

  calculateClassicalProbability(
    energy: number,
    mass: number,
    xGrid: number[],
  ): number[] {
    return calculateHarmonicOscillatorClassicalProbability(
      this.springConstant,
      energy,
      mass,
      xGrid,
    );
  }

  calculateWavefunctionZeros(stateIndex: number, _energy: number): number[] {
    return calculateHarmonicOscillatorWavefunctionZeros(
      this.springConstant,
      this.mass,
      stateIndex,
    );
  }

  calculateTurningPoints(
    energy: number,
  ): Array<{ left: number; right: number }> {
    const points = calculateHarmonicOscillatorTurningPoints(
      this.springConstant,
      energy,
    );
    return [points]; // Return as array with single element for simple single-well potential
  }

  calculateWavefunctionFirstDerivative(
    stateIndex: number,
    xGrid: number[],
  ): number[] {
    return calculateHarmonicOscillatorWavefunctionFirstDerivative(
      this.springConstant,
      this.mass,
      stateIndex,
      xGrid,
    );
  }

  calculateWavefunctionSecondDerivative(
    stateIndex: number,
    xGrid: number[],
  ): number[] {
    return calculateHarmonicOscillatorWavefunctionSecondDerivative(
      this.springConstant,
      this.mass,
      stateIndex,
      xGrid,
    );
  }

  calculateWavefunctionMinMax(
    stateIndex: number,
    xMin: number,
    xMax: number,
    numPoints?: number,
  ): { min: number; max: number } {
    return calculateHarmonicOscillatorWavefunctionMinMax(
      this.springConstant,
      this.mass,
      stateIndex,
      xMin,
      xMax,
      numPoints,
    );
  }

  calculateSuperpositionMinMax(
    coefficients: Array<[number, number]>,
    energies: number[],
    time: number,
    xMin: number,
    xMax: number,
    numPoints?: number,
  ): { min: number; max: number } {
    return calculateHarmonicOscillatorSuperpositionMinMax(
      this.springConstant,
      this.mass,
      coefficients,
      energies,
      time,
      xMin,
      xMax,
      numPoints,
    );
  }
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
