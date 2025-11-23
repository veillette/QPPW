/**
 * RK45 (Runge-Kutta-Fehlberg) method for computing wavefunctions from known energies.
 *
 * Once eigenvalues (energies) have been found using matrix methods (DVR, Spectral, FGH,
 * Matrix Numerov), this solver can compute the corresponding wavefunctions on an
 * arbitrarily fine grid using a standard ODE integration approach.
 *
 * The Schrödinger equation:
 *   -ℏ²/(2m) d²ψ/dx² + V(x)ψ = Eψ
 *
 * Can be rearranged to:
 *   d²ψ/dx² = (2m/ℏ²)(V(x) - E)ψ
 *
 * This is then solved as a system of first-order ODEs using RK45.
 */

import QuantumConstants from "./QuantumConstants.js";
import { GridConfig, PotentialFunction } from "./PotentialFunction.js";
import { normalizeWavefunction } from "./LinearAlgebraUtils.js";
import qppw from "../../QPPWNamespace.js";

/**
 * Result from computing wavefunctions using RK45
 */
export interface WavefunctionResult {
  /** Computed wavefunctions (each row is one eigenstate) */
  wavefunctions: number[][];
  /** Grid x-positions in meters */
  xGrid: number[];
}

/**
 * Compute wavefunctions for given energies using RK45 integration.
 *
 * Uses a shooting method from both boundaries, matching at a classical
 * turning point for numerical stability.
 *
 * @param energies - Energy eigenvalues in Joules (already found by matrix solver)
 * @param potential - Function V(x) that returns potential energy in Joules
 * @param mass - Particle mass in kg
 * @param gridConfig - Grid configuration (can have many more points than used to find energies)
 * @returns Wavefunctions on the specified grid
 */
export function computeWavefunctionsRK45(
  energies: number[],
  potential: PotentialFunction,
  mass: number,
  gridConfig: GridConfig,
): WavefunctionResult {
  const { xMin, xMax, numPoints } = gridConfig;
  const dx = (xMax - xMin) / (numPoints - 1);

  // Generate grid
  const xGrid: number[] = [];
  for (let i = 0; i < numPoints; i++) {
    xGrid.push(xMin + i * dx);
  }

  // Compute each wavefunction
  const wavefunctions: number[][] = [];
  for (const energy of energies) {
    const psi = computeSingleWavefunction(energy, potential, mass, xGrid, dx);
    const normalizedPsi = normalizeWavefunction(psi, dx);
    wavefunctions.push(normalizedPsi);
  }

  return {
    wavefunctions,
    xGrid,
  };
}

/**
 * Compute a single wavefunction for a given energy using RK45.
 *
 * Uses shooting from both boundaries and matching at a turning point.
 */
function computeSingleWavefunction(
  energy: number,
  potential: PotentialFunction,
  mass: number,
  xGrid: number[],
  dx: number,
): number[] {
  const { HBAR } = QuantumConstants;
  const coeff = (2 * mass) / (HBAR * HBAR);
  const N = xGrid.length;

  // Find matching point (classical turning point or middle)
  const matchIndex = findMatchingPoint(energy, potential, xGrid);

  // f(x) = (2m/ℏ²)(V(x) - E)
  const f = (x: number) => coeff * (potential(x) - energy);

  // Shoot from left boundary
  const psiLeft = shootFromLeft(xGrid, dx, f, matchIndex);

  // Shoot from right boundary
  const psiRight = shootFromRight(xGrid, dx, f, matchIndex);

  // Match solutions at matching point
  const psi = matchSolutions(psiLeft, psiRight, matchIndex, N);

  return psi;
}

/**
 * Find a good matching point (preferably near a classical turning point).
 */
function findMatchingPoint(
  energy: number,
  potential: PotentialFunction,
  xGrid: number[],
): number {
  const N = xGrid.length;

  // Look for classical turning point where V(x) = E
  // Start from the middle and search outward
  const middle = Math.floor(N / 2);

  // Search for right turning point (from middle going right)
  let rightTurning = -1;
  for (let i = middle; i < N - 1; i++) {
    if (potential(xGrid[i]) <= energy && potential(xGrid[i + 1]) > energy) {
      rightTurning = i;
      break;
    }
  }

  // Search for left turning point (from middle going left)
  let leftTurning = -1;
  for (let i = middle; i > 0; i--) {
    if (potential(xGrid[i]) <= energy && potential(xGrid[i - 1]) > energy) {
      leftTurning = i;
      break;
    }
  }

  // Choose matching point
  if (rightTurning >= 0 && leftTurning >= 0) {
    // Use the middle of the classical region
    return Math.floor((leftTurning + rightTurning) / 2);
  } else if (rightTurning >= 0) {
    return rightTurning;
  } else if (leftTurning >= 0) {
    return leftTurning;
  } else {
    // No turning point found, use middle
    return middle;
  }
}

/**
 * Shoot from left boundary using RK45.
 *
 * Solves the system:
 *   y1' = y2
 *   y2' = f(x) * y1
 *
 * where y1 = ψ, y2 = ψ'
 */
function shootFromLeft(
  xGrid: number[],
  dx: number,
  f: (x: number) => number,
  matchIndex: number,
): number[] {
  const psi: number[] = new Array(matchIndex + 1).fill(0);

  // Initial conditions: ψ(xMin) = 0, ψ'(xMin) = small value
  let y1 = 0;
  let y2 = 1e-10;
  psi[0] = y1;

  // Integrate forward using RK45
  for (let i = 0; i < matchIndex; i++) {
    const x = xGrid[i];
    const result = rk45Step(x, y1, y2, dx, f);
    y1 = result.y1;
    y2 = result.y2;
    psi[i + 1] = y1;
  }

  return psi;
}

/**
 * Shoot from right boundary using RK45.
 */
function shootFromRight(
  xGrid: number[],
  dx: number,
  f: (x: number) => number,
  matchIndex: number,
): number[] {
  const N = xGrid.length;
  const psi: number[] = new Array(N - matchIndex).fill(0);

  // Initial conditions: ψ(xMax) = 0, ψ'(xMax) = small value
  let y1 = 0;
  let y2 = 1e-10;
  psi[psi.length - 1] = y1;

  // Integrate backward using RK45 (negative step)
  for (let i = N - 1; i > matchIndex; i--) {
    const x = xGrid[i];
    const result = rk45Step(x, y1, y2, -dx, f);
    y1 = result.y1;
    y2 = result.y2;
    psi[i - matchIndex - 1] = y1;
  }

  return psi;
}

/**
 * Match left and right solutions at the matching point.
 */
function matchSolutions(
  psiLeft: number[],
  psiRight: number[],
  matchIndex: number,
  N: number,
): number[] {
  const psi: number[] = new Array(N).fill(0);

  // Get values at matching point
  const leftValue = psiLeft[psiLeft.length - 1];
  const rightValue = psiRight[0];

  // Scale factor to match amplitudes
  const scale = Math.abs(rightValue) > 1e-30 ? leftValue / rightValue : 1;

  // Copy left part
  for (let i = 0; i <= matchIndex; i++) {
    psi[i] = psiLeft[i];
  }

  // Copy scaled right part
  for (let i = matchIndex + 1; i < N; i++) {
    psi[i] = psiRight[i - matchIndex] * scale;
  }

  return psi;
}

/**
 * Single RK45 step for the second-order ODE system:
 *   y1' = y2
 *   y2' = f(x) * y1
 */
function rk45Step(
  x: number,
  y1: number,
  y2: number,
  h: number,
  f: (x: number) => number,
): { y1: number; y2: number } {
  // Butcher tableau coefficients for RK45
  const a2 = 1 / 4,
    a3 = 3 / 8,
    a4 = 12 / 13,
    a5 = 1,
    a6 = 1 / 2;
  const b21 = 1 / 4;
  const b31 = 3 / 32,
    b32 = 9 / 32;
  const b41 = 1932 / 2197,
    b42 = -7200 / 2197,
    b43 = 7296 / 2197;
  const b51 = 439 / 216,
    b52 = -8,
    b53 = 3680 / 513,
    b54 = -845 / 4104;
  const b61 = -8 / 27,
    b62 = 2,
    b63 = -3544 / 2565,
    b64 = 1859 / 4104,
    b65 = -11 / 40;

  // 4th order coefficients
  const c1 = 25 / 216,
    c3 = 1408 / 2565,
    c4 = 2197 / 4104,
    c5 = -1 / 5;

  // k1
  const k1_1 = h * y2;
  const k1_2 = h * f(x) * y1;

  // k2
  const y1_2 = y1 + b21 * k1_1;
  const y2_2 = y2 + b21 * k1_2;
  const k2_1 = h * y2_2;
  const k2_2 = h * f(x + a2 * h) * y1_2;

  // k3
  const y1_3 = y1 + b31 * k1_1 + b32 * k2_1;
  const y2_3 = y2 + b31 * k1_2 + b32 * k2_2;
  const k3_1 = h * y2_3;
  const k3_2 = h * f(x + a3 * h) * y1_3;

  // k4
  const y1_4 = y1 + b41 * k1_1 + b42 * k2_1 + b43 * k3_1;
  const y2_4 = y2 + b41 * k1_2 + b42 * k2_2 + b43 * k3_2;
  const k4_1 = h * y2_4;
  const k4_2 = h * f(x + a4 * h) * y1_4;

  // k5
  const y1_5 = y1 + b51 * k1_1 + b52 * k2_1 + b53 * k3_1 + b54 * k4_1;
  const y2_5 = y2 + b51 * k1_2 + b52 * k2_2 + b53 * k3_2 + b54 * k4_2;
  const k5_1 = h * y2_5;
  const k5_2 = h * f(x + a5 * h) * y1_5;

  // k6 (computed for 5th-order accuracy but not used in 4th-order solution)
  // These intermediate values are needed for the full RK45 embedded method
  const _y1_6 =
    y1 + b61 * k1_1 + b62 * k2_1 + b63 * k3_1 + b64 * k4_1 + b65 * k5_1;
  const _y2_6 =
    y2 + b61 * k1_2 + b62 * k2_2 + b63 * k3_2 + b64 * k4_2 + b65 * k5_2;
  // Suppress unused variable warnings - these are kept for potential future 5th-order error estimation
  void _y1_6;
  void _y2_6;
  void a6;

  // 4th order solution
  const y1_new = y1 + c1 * k1_1 + c3 * k3_1 + c4 * k4_1 + c5 * k5_1;
  const y2_new = y2 + c1 * k1_2 + c3 * k3_2 + c4 * k4_2 + c5 * k5_2;

  return { y1: y1_new, y2: y2_new };
}

/**
 * Convenience function to compute wavefunctions for a single energy.
 *
 * @param energy - Energy eigenvalue in Joules
 * @param potential - Function V(x) that returns potential energy in Joules
 * @param mass - Particle mass in kg
 * @param gridConfig - Grid configuration
 * @returns Single wavefunction on the specified grid
 */
export function computeWavefunctionRK45(
  energy: number,
  potential: PotentialFunction,
  mass: number,
  gridConfig: GridConfig,
): { wavefunction: number[]; xGrid: number[] } {
  const result = computeWavefunctionsRK45(
    [energy],
    potential,
    mass,
    gridConfig,
  );
  return {
    wavefunction: result.wavefunctions[0],
    xGrid: result.xGrid,
  };
}

qppw.register("WavefunctionRK45Solver", {
  computeWavefunctionsRK45,
  computeWavefunctionRK45,
});
