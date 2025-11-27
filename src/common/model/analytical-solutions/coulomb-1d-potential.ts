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
 * The 1D Coulomb problem has been studied in the physics literature, with particular
 * attention to the requirement that only odd-parity states are normalizable due to the
 * boundary conditions at the singular point x=0. The energy eigenvalues and wavefunctions
 * are derived from the Schrödinger equation with appropriate boundary conditions.
 *
 * FREELY AVAILABLE REFERENCES:
 * - Cheng, K.-M., & Lam, C. S. (2009). "The one-dimensional Coulomb Problem."
 *   arXiv:0905.3978. https://arxiv.org/abs/0905.3978
 *   Studies scattering and bound states for the 1D Coulomb potential V(x) = λ/|x|.
 *
 * - Campiglia, M., et al. (2019). "A Distributional Approach for the One-Dimensional Hydrogen Atom."
 *   Frontiers in Physics, 7, 101. https://www.frontiersin.org/articles/10.3389/fphy.2019.00101/full
 *   Addresses the non-integrable singularity at the origin using distribution theory.
 *
 * ENERGY EIGENVALUES:
 *   E_n = -mα²/(2ℏ²(n+1/2)²),  n = 0, 1, 2, ...
 *   (Note the half-integer quantum numbers, different from 3D Hydrogen)
 *
 * WAVEFUNCTIONS (odd parity only):
 *   ψ_n(x) = N_n · sign(x) · ρ · exp(-ρ/2) · L_n^1(ρ)
 *   where ρ = 2|x|/(n+1/2)a₀, a₀ = ℏ²/(mα), and L_n^1 are associated Laguerre polynomials.
 *   The sign(x) factor ensures odd parity: ψ(-x) = -ψ(x)
 *   The ρ factor ensures linear behavior near origin: ψ(x) ≈ Cx as x → 0
 */

import QuantumConstants from "../QuantumConstants.js";
import {
  BoundStateResult,
  GridConfig,
  PotentialFunction,
  FourierTransformResult,
} from "../PotentialFunction.js";
import { associatedLaguerre } from "./math-utilities.js";
import { AnalyticalSolution } from "./AnalyticalSolution.js";
import { computeNumericalFourierTransform } from "./fourier-transform-helper.js";

/**
 * Class-based implementation of 1D Coulomb potential analytical solution.
 * Extends the AnalyticalSolution abstract base class.
 */
export class Coulomb1DPotentialSolution extends AnalyticalSolution {
  constructor(
    private coulombStrength: number,
    private mass: number,
  ) {
    super();
  }

  solve(numStates: number, gridConfig: GridConfig): BoundStateResult {
    return solveCoulomb1DPotential(
      this.coulombStrength,
      this.mass,
      numStates,
      gridConfig,
    );
  }

  createPotential(): PotentialFunction {
    return createCoulomb1DPotential(this.coulombStrength);
  }

  calculateClassicalProbability(
    energy: number,
    mass: number,
    xGrid: number[],
  ): number[] {
    return calculateCoulomb1DClassicalProbability(
      this.coulombStrength,
      energy,
      mass,
      xGrid,
    );
  }

  calculateWavefunctionZeros(stateIndex: number, _energy: number): number[] {
    return calculateCoulomb1DWavefunctionZeros(
      this.coulombStrength,
      this.mass,
      stateIndex,
    );
  }

  calculateTurningPoints(
    energy: number,
  ): Array<{ left: number; right: number }> {
    const points = calculateCoulomb1DTurningPoints(
      this.coulombStrength,
      energy,
    );
    return [points]; // Return as array with single element for simple single-well potential
  }

  calculateWavefunctionFirstDerivative(
    stateIndex: number,
    xGrid: number[],
  ): number[] {
    return calculateCoulomb1DWavefunctionFirstDerivative(
      this.coulombStrength,
      this.mass,
      stateIndex,
      xGrid,
    );
  }

  calculateWavefunctionSecondDerivative(
    stateIndex: number,
    xGrid: number[],
  ): number[] {
    return calculateCoulomb1DWavefunctionSecondDerivative(
      this.coulombStrength,
      this.mass,
      stateIndex,
      xGrid,
    );
  }
  calculateFourierTransform(
    boundStateResult: BoundStateResult,
    mass: number,
    numMomentumPoints?: number,
    pMax?: number,
  ): FourierTransformResult {
    return computeNumericalFourierTransform(
      boundStateResult,
      mass,
      Math.abs(boundStateResult.energies[0]),
      numMomentumPoints,
      pMax,
    );
  }

}

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
    const energy =
      -(mass * alpha * alpha) / (2 * HBAR * HBAR * (n + 0.5) * (n + 0.5));
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
    // For ψ_n(x) = N * ρ * exp(-ρ/2) * L_n^1(ρ) where ρ = 2|x|/a_n
    // with ∫_{-∞}^{∞} |ψ|² dx = 1
    // Computing: ∫_{-∞}^{∞} |sign(x) * N * ρ * exp(-ρ/2) * L_n^1(ρ)|² dx
    //          = 2 ∫_0^∞ |N * ρ * exp(-ρ/2) * L_n^1(ρ)|² * (a_n/2) dρ
    //          = N² * a_n * ∫_0^∞ ρ² * exp(-ρ) * |L_n^1(ρ)|² dρ
    // For L_n^1, the integral ∫_0^∞ ρ² * exp(-ρ) * |L_n^1(ρ)|² dρ = 2(n+1)²
    // Therefore: N² * a_n * 2(n+1)² = 1
    const normalization = Math.sqrt(1.0 / (2 * a_n * (n + 1) * (n + 1)));

    for (const x of xGrid) {
      const absX = Math.abs(x);
      const rho = (2 * absX) / a_n;

      // Wavefunction: ψ(x) = sign(x) * ρ * N * exp(-ρ/2) * L_n^1(ρ)
      // The factor of ρ ensures linear behavior near x=0: ψ(x) ≈ C*x
      // ODD parity: ψ(-x) = -ψ(x), which matches the energy formula E_n = -E_R/(n+1/2)²
      const laguerre = associatedLaguerre(n, 1, rho);
      const radialPart = normalization * rho * Math.exp(-rho / 2) * laguerre;

      // Apply sign(x) for odd parity
      // Since radialPart ~ ρ ~ |x| near origin, and sign(x) * |x| = x,
      // we get ψ(x) ~ x near the origin (linear behavior)
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

/**
 * Create the potential function for a 1D Coulomb potential.
 * V(x) = -α/|x|
 *
 * @param coulombStrength - Coulomb strength parameter α in J·m
 * @returns Potential function V(x) in Joules
 */
export function createCoulomb1DPotential(
  coulombStrength: number,
): (x: number) => number {
  const alpha = coulombStrength;

  return (x: number) => {
    if (x === 0) {
      return -Infinity; // Singularity at x=0
    }
    return -alpha / Math.abs(x);
  };
}

/**
 * Calculate classical probability density for a 1D Coulomb potential.
 * P(x) ∝ 1/v(x) = 1/√[2(E - V(x))/m] = 1/√[2(E + α/|x|)/m]
 *
 * @param coulombStrength - Coulomb strength parameter α in J·m
 * @param energy - Energy of the particle in Joules (negative for bound states)
 * @param mass - Particle mass in kg
 * @param xGrid - Array of x positions in meters
 * @returns Array of normalized classical probability density values (in 1/meters)
 */
export function calculateCoulomb1DClassicalProbability(
  coulombStrength: number,
  energy: number,
  mass: number,
  xGrid: number[],
): number[] {
  const alpha = coulombStrength;
  const classicalProbability: number[] = [];
  let integralSum = 0;

  // Calculate unnormalized probability
  for (let i = 0; i < xGrid.length; i++) {
    const x = xGrid[i];
    const absX = Math.abs(x);

    if (absX < 1e-15) {
      // Near singularity, probability approaches infinity
      classicalProbability.push(0);
      continue;
    }

    const potential = -alpha / absX;
    const kineticEnergy = energy - potential;

    if (kineticEnergy <= 0) {
      classicalProbability.push(0);
    } else {
      const epsilon = 1e-10 * Math.abs(alpha);
      const probability =
        1 / Math.sqrt((2 * Math.max(kineticEnergy, epsilon)) / mass);
      classicalProbability.push(probability);

      if (i > 0) {
        const dx = xGrid[i] - xGrid[i - 1];
        integralSum += ((probability + classicalProbability[i - 1]) * dx) / 2;
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
 * Calculate the classical turning points for a 1D Coulomb potential.
 * Solve E = -α/|x| for x: |x| = -α/E => x = ±α/|E|
 *
 * @param coulombStrength - Coulomb strength parameter α in J·m
 * @param energy - Energy of the particle in Joules (negative for bound states)
 * @returns Object with left and right turning point positions (in meters)
 */
export function calculateCoulomb1DTurningPoints(
  coulombStrength: number,
  energy: number,
): { left: number; right: number } {
  const alpha = coulombStrength;

  // For bound states (E < 0): E = -α/|x| => |x| = -α/E
  const turningRadius = -alpha / energy;

  return {
    left: -turningRadius,
    right: turningRadius,
  };
}

/**
 * Calculate wavefunction zeros for 1D Coulomb potential (numerical approach).
 * The wavefunction is ψ_n(x) = sign(x) * ρ * exp(-ρ/2) * L_n^1(ρ)
 * where ρ = 2|x|/a_n. Zeros occur where L_n^1(ρ) = 0 (excluding ρ=0).
 *
 * @param coulombStrength - Coulomb strength parameter α in J·m
 * @param mass - Particle mass in kg
 * @param stateIndex - Index of the eigenstate (0 for ground state, etc.)
 * @param searchRange - Range to search for zeros (in meters)
 * @returns Array of x positions (in meters) where wavefunction is zero
 */
export function calculateCoulomb1DWavefunctionZeros(
  coulombStrength: number,
  mass: number,
  stateIndex: number,
  searchRange: number = 20e-9,
): number[] {
  const { HBAR } = QuantumConstants;
  const alpha = coulombStrength;
  const n = stateIndex;
  const nEff = n + 0.5;
  const a0 = (HBAR * HBAR) / (mass * alpha);
  const a_n = nEff * a0;

  // Ground state (n=0) has no interior zeros (L_0^1(ρ) = 1 - ρ has one zero at ρ=1)
  // but the ρ factor gives a zero at x=0, which we don't count as an interior zero

  const zeros: number[] = [];
  const numSamples = 1000;
  const dx = (2 * searchRange) / numSamples;

  // Search in positive x only (use symmetry for negative x)
  let prevX = 1e-15; // Start just above zero to avoid singularity
  const prevRho = (2 * prevX) / a_n;
  let prevVal =
    prevRho * Math.exp(-prevRho / 2) * associatedLaguerre(n, 1, prevRho);

  for (let i = 1; i <= numSamples / 2; i++) {
    const x = prevX + i * dx;
    const rho = (2 * x) / a_n;
    const laguerre = associatedLaguerre(n, 1, rho);
    const val = rho * Math.exp(-rho / 2) * laguerre;

    // Sign change detected
    if (prevVal * val < 0 && x > 1e-14) {
      // Use bisection to refine
      let left = prevX;
      let right = x;
      for (let iter = 0; iter < 20; iter++) {
        const mid = (left + right) / 2;
        const midRho = (2 * mid) / a_n;
        const midLaguerre = associatedLaguerre(n, 1, midRho);
        const valMid = midRho * Math.exp(-midRho / 2) * midLaguerre;

        if (Math.abs(valMid) < 1e-12) {
          // Found zero on positive side, add both ±x
          zeros.push(-mid); // Negative side
          zeros.push(mid); // Positive side
          break;
        }

        if (valMid * prevVal < 0) {
          right = mid;
        } else {
          left = mid;
        }

        if (iter === 19) {
          const midFinal = (left + right) / 2;
          zeros.push(-midFinal);
          zeros.push(midFinal);
        }
      }
    }

    prevX = x;
    prevVal = val;
  }

  // Sort zeros
  zeros.sort((a, b) => a - b);
  return zeros;
}

/**
 * Calculate the first derivative of the wavefunction for a 1D Coulomb potential.
 * Uses numerical differentiation on the analytical wavefunction.
 *
 * @param coulombStrength - Coulomb strength parameter α in J·m
 * @param mass - Particle mass in kg
 * @param stateIndex - Index of the eigenstate (0 for ground state, etc.)
 * @param xGrid - Array of x positions in meters where derivatives should be evaluated
 * @returns Array of first derivative values
 */
export function calculateCoulomb1DWavefunctionFirstDerivative(
  coulombStrength: number,
  mass: number,
  stateIndex: number,
  xGrid: number[],
): number[] {
  const { HBAR } = QuantumConstants;
  const alpha = coulombStrength;
  const n = stateIndex;
  const nEff = n + 0.5;
  const a0 = (HBAR * HBAR) / (mass * alpha);
  const a_n = nEff * a0;
  const normalization = Math.sqrt(1.0 / (2 * a_n * (n + 1) * (n + 1)));

  const firstDerivative: number[] = [];
  const h = 1e-12; // Small step for numerical differentiation

  // Helper function to evaluate wavefunction
  const evaluatePsi = (x: number): number => {
    const absX = Math.abs(x);
    const rho = (2 * absX) / a_n;
    const laguerre = associatedLaguerre(n, 1, rho);
    const radialPart = normalization * rho * Math.exp(-rho / 2) * laguerre;
    return Math.sign(x) * radialPart;
  };

  for (const x of xGrid) {
    // Handle near-singularity carefully
    if (Math.abs(x) < 2 * h) {
      // Very close to singularity - use the fact that ψ(x) ~ x near origin
      // So ψ'(0) should be approximately constant (normalization factor)
      firstDerivative.push(normalization);
      continue;
    }

    // Evaluate at x-h and x+h
    const xMinus = x - h;
    const xPlus = x + h;

    const psiMinus = evaluatePsi(xMinus);
    const psiPlus = evaluatePsi(xPlus);

    // First derivative using central difference
    const firstDeriv = (psiPlus - psiMinus) / (2 * h);
    firstDerivative.push(firstDeriv);
  }

  return firstDerivative;
}

/**
 * Calculate the second derivative of the wavefunction for a 1D Coulomb potential.
 * Uses numerical differentiation on the analytical wavefunction.
 *
 * @param coulombStrength - Coulomb strength parameter α in J·m
 * @param mass - Particle mass in kg
 * @param stateIndex - Index of the eigenstate (0 for ground state, etc.)
 * @param xGrid - Array of x positions in meters where derivatives should be evaluated
 * @returns Array of second derivative values
 */
export function calculateCoulomb1DWavefunctionSecondDerivative(
  coulombStrength: number,
  mass: number,
  stateIndex: number,
  xGrid: number[],
): number[] {
  const { HBAR } = QuantumConstants;
  const alpha = coulombStrength;
  const n = stateIndex;
  const nEff = n + 0.5;
  const a0 = (HBAR * HBAR) / (mass * alpha);
  const a_n = nEff * a0;
  const normalization = Math.sqrt(1.0 / (2 * a_n * (n + 1) * (n + 1)));

  const secondDerivative: number[] = [];
  const h = 1e-12; // Small step for numerical differentiation

  // Helper function to evaluate wavefunction
  const evaluatePsi = (x: number): number => {
    const absX = Math.abs(x);
    const rho = (2 * absX) / a_n;
    const laguerre = associatedLaguerre(n, 1, rho);
    const radialPart = normalization * rho * Math.exp(-rho / 2) * laguerre;
    return Math.sign(x) * radialPart;
  };

  for (const x of xGrid) {
    // Handle near-singularity carefully
    if (Math.abs(x) < 2 * h) {
      // Very close to singularity - derivative is undefined/infinite
      // But we can use the fact that ψ(x) ~ x near origin
      // So ψ''(0) should be 0
      secondDerivative.push(0);
      continue;
    }

    // Evaluate at x-h, x, x+h
    const xMinus = x - h;
    const xPlus = x + h;

    const psiMinus = evaluatePsi(xMinus);
    const psi = evaluatePsi(x);
    const psiPlus = evaluatePsi(xPlus);

    // Second derivative using central difference
    const secondDeriv = (psiPlus - 2 * psi + psiMinus) / (h * h);
    secondDerivative.push(secondDeriv);
  }

  return secondDerivative;
}
