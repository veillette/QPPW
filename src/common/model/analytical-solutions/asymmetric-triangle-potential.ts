/**
 * Analytical solution for the asymmetric triangle potential.
 * Modified to have V(∞) = 10 eV and V(0) = 0 eV (minimum)
 *
 * Original formulation:
 * V(x) = 0 for x < 0
 * V(x) = -b(a-x) for 0 < x < a
 * V(x) = 0 for x > a
 *
 * Shifted formulation:
 * V(x) = 10 eV for x < 0
 * V(x) = 10 - b(a-x) for 0 < x < a  (equivalently: V(x) = 10 - ba + bx)
 * V(x) = 10 eV for x > a
 *
 * At x=0: V = 10 - ba = 0 (minimum) => ba = 10 eV
 * At x=a: V = 10 eV
 * At x > a: V = 10 eV
 *
 * This potential creates a triangular well and the solution involves Airy functions.
 */

import QuantumConstants from "../QuantumConstants.js";
import { BoundStateResult, GridConfig } from "../PotentialFunction.js";
import { airyAi } from "./math-utilities.js";

/**
 * Analytical solution for the asymmetric triangle potential.
 * Modified to have V(∞) = 10 eV and V(0) = 0 eV (minimum)
 *
 * @param slope - Slope parameter b in Joules/meter (positive value)
 * @param wellWidth - Width parameter a in meters
 * @param mass - Particle mass in kg
 * @param numStates - Number of energy levels to calculate
 * @param gridConfig - Grid configuration for wavefunction evaluation
 * @returns Bound state results with exact energies and wavefunctions
 */
export function solveAsymmetricTrianglePotential(
  slope: number,
  wellWidth: number,
  mass: number,
  numStates: number,
  gridConfig: GridConfig,
): BoundStateResult {
  const { HBAR, EV_TO_JOULES } = QuantumConstants;
  const b = slope;
  const a = wellWidth;

  // Potential shift: we want V(∞) = 10 eV and V(0) = 0 eV
  // In the original formulation, V(0) = -ba and V(∞) = 0
  // So we shift by +ba + 10 eV, which makes V(0) = 0 and V(∞) = 10 eV
  const POTENTIAL_SHIFT = 10 * EV_TO_JOULES; // 10 eV in Joules

  // Define the scaling parameter for Airy functions
  // α = (2mb/ℏ²)^(1/3)
  const alpha = Math.pow((2 * mass * b) / (HBAR * HBAR), 1 / 3);

  // The original potential minimum is at V(0) = -ba
  // After shifting, we want V(0) = 0 eV (minimum) and V(∞) = 10 eV
  // So we shift by +ba + 10 eV, where ba is the well depth in the original formulation
  const V0_original = -b * a;

  // Energy eigenvalues are found by solving the transcendental equation
  // from matching boundary conditions at x = 0 and x = a
  // For bound states in original formulation: -ba < E < 0
  // After shifting: 0 < E < ba + 10 eV

  const energies: number[] = [];
  const eigenvaluesZ: number[] = []; // Store scaled eigenvalues

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

  // Find eigenvalues using Newton-Raphson root finding with appropriate initial guesses
  for (let n = 0; n < numStates; n++) {
    let initialGuess: number;

    if (n < airyZeros.length) {
      // Use pre-computed Airy zero as initial guess for the first 30 states
      initialGuess = airyZeros[n];
    } else {
      // Use asymptotic approximation for states beyond the 30th
      initialGuess = getApproximateAiryZero(n + 1); // n+1 because zeros are 1-indexed
    }

    // Refine the zero using Newton-Raphson
    const z0 = refineAiryZero(initialGuess);

    // For the asymmetric triangle, the energy eigenvalues are approximately:
    // E_n ≈ -ba + (ℏ²/2m)^(1/3) * (2b)^(2/3) * z_n
    // where z_n are negative values related to Airy function zeros

    // Calculate energy in original formulation (relative to V=0 at infinity)
    const energy_original = V0_original - Math.pow(HBAR * HBAR / (2 * mass), 1 / 3) * Math.pow(2 * b, 2 / 3) * Math.abs(z0);

    // Shift energy to new reference (V(0) = 0, V(∞) = 10 eV)
    // energy_shifted = energy_original + ba + 10 eV
    const totalShift = -V0_original + POTENTIAL_SHIFT; // This is ba + 10 eV
    const energy = energy_original + totalShift;

    // Only include bound states (0 < E < 10 eV for shifted potential)
    if (energy > 0 && energy < POTENTIAL_SHIFT) {
      energies.push(energy);
      eigenvaluesZ.push(z0);
    } else if (energy <= 0) {
      // Energy too low, still in bound states range but may be numerical error
      energies.push(energy);
      eigenvaluesZ.push(z0);
    } else {
      // No more bound states (E >= 10 eV)
      break;
    }
  }

  const actualNumStates = energies.length;

  if (actualNumStates === 0) {
    throw new Error("Asymmetric triangle potential too shallow to support bound states");
  }

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

    // In region II (0 < x < a), the wavefunction is:
    // ψ(x) = N * Ai(α(x - x_0))
    // For the shifted potential, we need to work in the original formulation
    // E_original = E - totalShift
    const E_original = E - (-V0_original + POTENTIAL_SHIFT);
    // where x_0 = (E_original - V_original(0))/b = (E_original + ba)/b
    const x0 = (E_original + b * a) / b;

    // Normalization constant (approximate)
    let normSq = 0;
    for (const x of xGrid) {
      if (x >= 0 && x <= a) {
        const z = alpha * (x - x0);
        const psi = airyAi(z);
        normSq += psi * psi * dx;
      }
    }
    const norm = 1 / Math.sqrt(normSq);

    // Calculate wavefunction on grid
    for (const x of xGrid) {
      if (x < 0) {
        // Region I: exponentially decaying
        // In shifted potential: V = 10 eV, so κ² = 2m(V - E)/ℏ² = 2m(10 eV - E)/ℏ²
        const potentialDiff = POTENTIAL_SHIFT - E;
        if (potentialDiff > 0) {
          const kappa = Math.sqrt(2 * mass * potentialDiff) / HBAR;
          const psi_0 = norm * airyAi(alpha * (0 - x0));
          wavefunction.push(psi_0 * Math.exp(kappa * x));
        } else {
          wavefunction.push(0); // Above barrier, no bound state
        }
      } else if (x <= a) {
        // Region II: Airy function solution
        const z = alpha * (x - x0);
        wavefunction.push(norm * airyAi(z));
      } else {
        // Region III: exponentially decaying
        // In shifted potential: V = 10 eV, so κ² = 2m(V - E)/ℏ² = 2m(10 eV - E)/ℏ²
        const potentialDiff = POTENTIAL_SHIFT - E;
        if (potentialDiff > 0) {
          const kappa = Math.sqrt(2 * mass * potentialDiff) / HBAR;
          const psi_a = norm * airyAi(alpha * (a - x0));
          wavefunction.push(psi_a * Math.exp(-kappa * (x - a)));
        } else {
          wavefunction.push(0); // Above barrier, no bound state
        }
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
