/**
 * Analytical solution for the asymmetric triangle potential.
 * V(x) = V₀ for x < 0
 * V(x) = F·x for 0 < x < a (linear from 0 to F·a)
 * V(x) = V₀ for x > a
 *
 * Where V₀ = F·a (the barrier height)
 * This potential creates a triangular well with the minimum at x=0 where V=0.
 * The solution involves Airy functions.
 */

import QuantumConstants from "../QuantumConstants.js";
import { BoundStateResult, GridConfig } from "../PotentialFunction.js";
import { airyAi } from "./math-utilities.js";

/**
 * Analytical solution for the asymmetric triangle potential with infinite walls.
 * V(x) = ∞ for x < 0
 * V(x) = F·x for 0 < x (linear increasing potential)
 *
 * This is the standard triangular well problem. The eigenvalues are related to
 * the zeros of the Airy function Ai(z).
 *
 * @param slope - Slope parameter F in Joules/meter (field strength)
 * @param wellWidth - Width parameter a in meters (not used for infinite well, but kept for compatibility)
 * @param mass - Particle mass in kg
 * @param numStates - Number of energy levels to calculate
 * @param gridConfig - Grid configuration for wavefunction evaluation
 * @returns Bound state results with exact energies and wavefunctions
 */
export function solveAsymmetricTrianglePotential(
  slope: number,
  _wellWidth: number, // Kept for API compatibility but not used in infinite wall formulation
  mass: number,
  numStates: number,
  gridConfig: GridConfig,
): BoundStateResult {
  const { HBAR } = QuantumConstants;
  const F = slope; // Field strength

  // Define the scaling parameter for Airy functions
  // α = (2mF/ℏ²)^(1/3)
  const alpha = Math.pow((2 * mass * F) / (HBAR * HBAR), 1 / 3);

  // Energy eigenvalues for infinite triangular well:
  // E_n = -z_n · (ℏ²/2m)^(1/3) · F^(2/3)
  // where z_n are the zeros of Ai(z) (negative values)

  const energies: number[] = [];
  const eigenvaluesZ: number[] = []; // Store Airy zeros

  // The first 30 zeros of Ai(z) (accurate values)
  const airyZeros = [
    -2.33810741046, -4.08794944413, -5.52055982810, -6.78670809007,
    -7.94413358712, -9.02265085334, -10.0401743416, -11.0085243037,
    -11.9360155632, -12.8287767529, -13.6914890352, -14.5278299518,
    -15.3407550016, -16.1328355283, -16.9062914467, -17.6629059623,
    -18.4040686682, -19.1316436363, -19.8471121589, -20.5516776771,
    -21.2463345333, -21.9319076370, -22.6091606826, -23.2787913810,
    -23.9414358498, -24.5976826013, -25.2480724088, -25.8930940959,
    -26.5331934669, -27.1687726548
  ];

  /**
   * Asymptotic approximation for Airy zeros when n >= 30
   * z_n ≈ -(3π(4n-1)/8)^(2/3)
   */
  function getApproximateAiryZero(n: number): number {
    return -Math.pow((3 * Math.PI * (4 * n - 1)) / 8, 2 / 3);
  }

  /**
   * Refine Airy zero using Newton-Raphson method
   * We're looking for zeros of Ai(z), so we need Ai(z) = 0
   */
  function refineAiryZero(initialGuess: number, maxIterations: number = 10, tolerance: number = 1e-10): number {
    let z = initialGuess;
    const h = 1e-6; // Small step for numerical derivative

    for (let iter = 0; iter < maxIterations; iter++) {
      const aiZ = airyAi(z);

      // Check convergence
      if (Math.abs(aiZ) < tolerance) {
        break;
      }

      // Numerical derivative: Ai'(z) ≈ (Ai(z+h) - Ai(z-h)) / (2h)
      const aiZPlusH = airyAi(z + h);
      const aiZMinusH = airyAi(z - h);
      const aiPrime = (aiZPlusH - aiZMinusH) / (2 * h);

      // Newton-Raphson step: z_new = z - Ai(z) / Ai'(z)
      if (Math.abs(aiPrime) < 1e-12) {
        // Derivative too small, can't continue
        break;
      }

      const zNew = z - aiZ / aiPrime;

      // Check for convergence
      if (Math.abs(zNew - z) < tolerance) {
        z = zNew;
        break;
      }

      z = zNew;
    }

    return z;
  }

  // Calculate energies using Airy zeros
  for (let n = 0; n < numStates; n++) {
    let z_n: number;

    if (n < airyZeros.length) {
      // Use pre-computed Airy zero
      z_n = airyZeros[n];
    } else {
      // Use asymptotic approximation for states beyond the 30th
      const initialGuess = getApproximateAiryZero(n + 1); // n+1 because zeros are 1-indexed
      z_n = refineAiryZero(initialGuess);
    }

    // Energy eigenvalues: E_n = -z_n · (ℏ²/2m)^(1/3) · F^(2/3)
    // Since z_n < 0, we have -z_n > 0, so E > 0
    const energy = -z_n * Math.pow(HBAR * HBAR / (2 * mass), 1 / 3) * Math.pow(F, 2 / 3);

    energies.push(energy);
    eigenvaluesZ.push(z_n);
  }

  const actualNumStates = energies.length;

  // Generate grid
  const numPoints = gridConfig.numPoints;
  const xGrid: number[] = [];
  const dx = (gridConfig.xMax - gridConfig.xMin) / (numPoints - 1);
  for (let i = 0; i < numPoints; i++) {
    xGrid.push(gridConfig.xMin + i * dx);
  }

  // Calculate wavefunctions using Airy functions
  const wavefunctions: number[][] = [];

  for (let n = 0; n < actualNumStates; n++) {
    const wavefunction: number[] = [];
    const E = energies[n];

    // Classical turning point: x_0 = E/F (where V(x_0) = F·x_0 = E)
    const x0 = E / F;

    // Normalization constant (approximate using numerical integration)
    let normSq = 0;
    for (const x of xGrid) {
      if (x >= 0) {
        const z = alpha * (x - x0);
        const psi = airyAi(z);
        normSq += psi * psi * dx;
      }
    }
    const norm = 1 / Math.sqrt(normSq);

    // Calculate wavefunction on grid
    for (const x of xGrid) {
      if (x < 0) {
        // Region x < 0: infinite wall, ψ = 0
        wavefunction.push(0);
      } else {
        // Region x ≥ 0: Airy function solution
        // ψ(x) = N · Ai(α(x - x_0))
        const z = alpha * (x - x0);
        wavefunction.push(norm * airyAi(z));
      }
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
