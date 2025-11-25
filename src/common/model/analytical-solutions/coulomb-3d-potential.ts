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

/**
 * Create the potential function for a 3D Coulomb potential.
 * V(r) = -α/r (using radial distance r = |x|)
 *
 * @param coulombStrength - Coulomb strength parameter α in J·m
 * @returns Potential function V(r) in Joules
 */
export function createCoulomb3DPotential(
  coulombStrength: number,
): (x: number) => number {
  const alpha = coulombStrength;

  return (x: number) => {
    const r = Math.abs(x);
    if (r === 0) {
      return -Infinity; // Singularity at r=0
    }
    return -alpha / r;
  };
}

/**
 * Calculate classical probability density for a 3D Coulomb potential (radial).
 * P(r) ∝ 1/v(r) = 1/√[2(E - V(r))/m] = 1/√[2(E + α/r)/m]
 *
 * @param coulombStrength - Coulomb strength parameter α in J·m
 * @param energy - Energy of the particle in Joules (negative for bound states)
 * @param mass - Particle mass in kg
 * @param xGrid - Array of x positions in meters (uses r = |x|)
 * @returns Array of normalized classical probability density values (in 1/meters)
 */
export function calculateCoulomb3DClassicalProbability(
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
    const r = Math.abs(xGrid[i]);

    if (r < 1e-15) {
      // Near singularity
      classicalProbability.push(0);
      continue;
    }

    const potential = -alpha / r;
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
 * Calculate the classical turning points for a 3D Coulomb potential.
 * Solve E = -α/r for r: r = -α/E
 *
 * @param coulombStrength - Coulomb strength parameter α in J·m
 * @param energy - Energy of the particle in Joules (negative for bound states)
 * @returns Object with left and right turning point positions (in meters)
 */
export function calculateCoulomb3DTurningPoints(
  coulombStrength: number,
  energy: number,
): { left: number; right: number } {
  const alpha = coulombStrength;

  // For bound states (E < 0): E = -α/r => r = -α/E
  const turningRadius = -alpha / energy;

  return {
    left: -turningRadius, // Symmetric for radial visualization
    right: turningRadius,
  };
}

/**
 * Calculate wavefunction zeros for 3D Coulomb potential (numerical approach).
 * The radial wavefunction is R_n0(r) = N * exp(-ρ/2) * L^1_(n-1)(ρ)
 * where ρ = 2r/(na₀). Zeros occur where L^1_(n-1)(ρ) = 0.
 *
 * @param coulombStrength - Coulomb strength parameter α in J·m
 * @param mass - Particle mass in kg
 * @param stateIndex - Index of the eigenstate (0 for n=1, 1 for n=2, etc.)
 * @param searchRange - Range to search for zeros (in meters)
 * @returns Array of x positions (in meters) where wavefunction is zero
 */
export function calculateCoulomb3DWavefunctionZeros(
  coulombStrength: number,
  mass: number,
  stateIndex: number,
  searchRange: number = 20e-9,
): number[] {
  const { HBAR } = QuantumConstants;
  const alpha = coulombStrength;
  const n = stateIndex + 1; // Convert index to principal quantum number (1, 2, 3, ...)
  const a0 = (HBAR * HBAR) / (mass * alpha);
  const a_n = n * a0;

  // Ground state (n=1) has no zeros (L^1_0(ρ) = 1, constant)
  if (n === 1) {
    return [];
  }

  const zeros: number[] = [];
  const numSamples = 1000;
  const dr = (2 * searchRange) / numSamples;

  // Search in positive r only (use symmetry for negative x)
  let prevR = 1e-15; // Start just above zero
  const prevRho = (2 * prevR) / a_n;
  let prevVal = Math.exp(-prevRho / 2) * associatedLaguerre(n - 1, 1, prevRho);

  for (let i = 1; i <= numSamples / 2; i++) {
    const r = prevR + i * dr;
    const rho = (2 * r) / a_n;
    const laguerre = associatedLaguerre(n - 1, 1, rho);
    const val = Math.exp(-rho / 2) * laguerre;

    // Sign change detected
    if (prevVal * val < 0 && r > 1e-14) {
      // Use bisection to refine
      let left = prevR;
      let right = r;
      for (let iter = 0; iter < 20; iter++) {
        const mid = (left + right) / 2;
        const midRho = (2 * mid) / a_n;
        const midLaguerre = associatedLaguerre(n - 1, 1, midRho);
        const valMid = Math.exp(-midRho / 2) * midLaguerre;

        if (Math.abs(valMid) < 1e-12) {
          // Found zero on positive side, add both ±r for symmetric visualization
          zeros.push(-mid);
          zeros.push(mid);
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

    prevR = r;
    prevVal = val;
  }

  // Sort zeros
  zeros.sort((a, b) => a - b);
  return zeros;
}

/**
 * Calculate the second derivative of the radial wavefunction for a 3D Coulomb potential.
 * Uses numerical differentiation on the analytical wavefunction.
 *
 * @param coulombStrength - Coulomb strength parameter α in J·m
 * @param mass - Particle mass in kg
 * @param stateIndex - Index of the eigenstate (0 for n=1, 1 for n=2, etc.)
 * @param xGrid - Array of x positions in meters where derivatives should be evaluated
 * @returns Array of second derivative values
 */
export function calculateCoulomb3DWavefunctionSecondDerivative(
  coulombStrength: number,
  mass: number,
  stateIndex: number,
  xGrid: number[],
): number[] {
  const { HBAR } = QuantumConstants;
  const alpha = coulombStrength;
  const n = stateIndex + 1; // Convert index to principal quantum number (1, 2, 3, ...)
  const a0 = (HBAR * HBAR) / (mass * alpha);
  const a_n = n * a0;
  const normalization =
    (2.0 / (a_n * Math.sqrt(a_n))) *
    Math.sqrt(factorial(n - 1) / (n * factorial(n)));

  const secondDerivative: number[] = [];
  const h = 1e-12; // Small step for numerical differentiation

  // Helper function to evaluate wavefunction
  const evaluatePsi = (x: number): number => {
    const r = Math.abs(x);
    if (r === 0) {
      return 0;
    }
    const rho = (2 * r) / a_n;
    const laguerre = associatedLaguerre(n - 1, 1, rho);
    return normalization * Math.exp(-rho / 2) * laguerre;
  };

  for (const x of xGrid) {
    // Handle near-singularity carefully
    if (Math.abs(x) < 2 * h) {
      // Very close to singularity
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
