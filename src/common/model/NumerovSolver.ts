/**
 * Numerov method for solving the 1D time-independent Schrödinger equation.
 * Uses the higher-order Numerov formula with O(h^6) error.
 *
 * The TISE is: -ℏ²/(2m) d²ψ/dx² + V(x)ψ = Eψ
 * Rearranged as: d²ψ/dx² = -k²(x)ψ where k²(x) = 2m(E - V(x))/ℏ²
 *
 * Numerov formula: ψ_(j+1) = [(12 - 10f_j)ψ_j - (1+f_(j-1))ψ_(j-1)] / (1+f_(j+1))
 * where f_j = (h²/12) k²(x_j)
 */

import QuantumConstants from "./QuantumConstants.js";
import { BoundStateResult, GridConfig, PotentialFunction } from "./PotentialFunction.js";
import qppw from "../../QPPWNamespace.js";

/**
 * Solve the 1D Schrödinger equation using the Numerov method.
 *
 * @param potential - Function V(x) that returns potential energy in Joules
 * @param mass - Particle mass in kg
 * @param numStates - Number of bound states to find
 * @param gridConfig - Grid configuration
 * @param energyMin - Minimum energy to search (Joules)
 * @param energyMax - Maximum energy to search (Joules)
 * @returns Bound state results
 */
export function solveNumerov(
  potential: PotentialFunction,
  mass: number,
  numStates: number,
  gridConfig: GridConfig,
  energyMin: number,
  energyMax: number,
): BoundStateResult {
  const { xMin, xMax, numPoints } = gridConfig;
  const dx = (xMax - xMin) / (numPoints - 1);

  // Generate grid
  const xGrid: number[] = [];
  for (let i = 0; i < numPoints; i++) {
    xGrid.push(xMin + i * dx);
  }

  // Evaluate potential on grid
  const V = xGrid.map(potential);

  // Find bound state energies using shooting method
  const energies: number[] = [];
  const wavefunctions: number[][] = [];

  // Search for energies where wavefunction satisfies boundary conditions
  const energyStep = (energyMax - energyMin) / 1000;
  let prevSign = 0;

  for (let E = energyMin; E <= energyMax && energies.length < numStates; E += energyStep) {
    const psi = integrateNumerov(E, V, xGrid, dx, mass);
    const endValue = psi[numPoints - 1];

    // Check for sign change (indicates bound state)
    const currentSign = Math.sign(endValue);
    if (prevSign !== 0 && currentSign !== prevSign) {
      // Refine energy using bisection
      const refinedEnergy = refineEnergy(
        E - energyStep,
        E,
        V,
        xGrid,
        dx,
        mass,
      );
      energies.push(refinedEnergy);

      // Calculate normalized wavefunction
      const refinedPsi = integrateNumerov(refinedEnergy, V, xGrid, dx, mass);
      const normalizedPsi = normalizeWavefunction(refinedPsi, dx);
      wavefunctions.push(normalizedPsi);
    }
    prevSign = currentSign;
  }

  return {
    energies,
    wavefunctions,
    xGrid,
    method: "numerov",
  };
}

/**
 * Integrate the Schrödinger equation using Numerov formula.
 *
 * @param E - Energy eigenvalue to test (Joules)
 * @param V - Potential energy array (Joules)
 * @param xGrid - Spatial grid (meters)
 * @param dx - Grid spacing (meters)
 * @param mass - Particle mass (kg)
 * @returns Wavefunction array
 */
export function integrateNumerov(
  E: number,
  V: number[],
  xGrid: number[],
  dx: number,
  mass: number,
): number[] {
  const { HBAR } = QuantumConstants;
  const N = xGrid.length;
  const psi = new Array(N).fill(0);

  // Calculate k²(x) = 2m(E - V(x))/ℏ²
  const k2 = V.map((v) => (2 * mass * (E - v)) / (HBAR * HBAR));

  // Calculate f_j = (h²/12) * k²(x_j)
  const f = k2.map((k) => (dx * dx / 12) * k);

  // Initial conditions (boundary condition: ψ(x_min) = 0)
  psi[0] = 0;
  psi[1] = dx; // Small non-zero value

  // Numerov forward integration
  for (let j = 1; j < N - 1; j++) {
    const numerator = (12 - 10 * f[j]) * psi[j] - (1 + f[j - 1]) * psi[j - 1];
    const denominator = 1 + f[j + 1];
    psi[j + 1] = numerator / denominator;

    // Check for divergence (not a bound state)
    if (Math.abs(psi[j + 1]) > 1e10) {
      // Force large value to indicate divergence
      for (let k = j + 1; k < N; k++) {
        psi[k] = psi[j + 1];
      }
      break;
    }
  }

  return psi;
}

/**
 * Integrate from boundaries inward for symmetric double well potentials.
 * Uses Numerov method with proper boundary conditions to avoid exponential growth.
 *
 * Key insight: For double wells, integrating from the center (in the barrier)
 * outward causes exponential growth. Instead, integrate from the boundaries
 * (where ψ=0) inward to the center, then match using parity conditions.
 *
 * @param E - Energy eigenvalue
 * @param V - Potential array
 * @param xGrid - Spatial grid (must be symmetric around x=0)
 * @param dx - Grid spacing
 * @param mass - Particle mass
 * @param parity - "symmetric" or "antisymmetric"
 * @returns Wavefunction array
 */
export function integrateNumerovFromCenter(
  E: number,
  V: number[],
  xGrid: number[],
  dx: number,
  mass: number,
  parity: "symmetric" | "antisymmetric",
): number[] {
  const { HBAR } = QuantumConstants;
  const N = xGrid.length;
  const psi = new Array(N).fill(0);

  // Find center index (closest to x=0)
  let centerIdx = 0;
  let minDist = Math.abs(xGrid[0]);
  for (let i = 1; i < N; i++) {
    const dist = Math.abs(xGrid[i]);
    if (dist < minDist) {
      minDist = dist;
      centerIdx = i;
    }
  }

  // Calculate k²(x) = 2m(E - V(x))/ℏ²
  const k2 = V.map((v) => (2 * mass * (E - v)) / (HBAR * HBAR));

  // Calculate f_j = (h²/12) * k²(x_j)
  const f = k2.map((k) => (dx * dx / 12) * k);

  // NEW APPROACH: Integrate from RIGHT boundary inward to center
  // Boundary condition: ψ(x_max) = 0, ψ(x_max - dx) = small value  // Use a larger initial value to avoid numerical issues with tiny numbers
  psi[N - 1] = 0;
  psi[N - 2] = 1e-3; // Small but not too tiny

  // Integrate leftward from right boundary to center
  for (let j = N - 2; j > centerIdx; j--) {
    const numerator = (12 - 10 * f[j]) * psi[j] - (1 + f[j + 1]) * psi[j + 1];
    const denominator = 1 + f[j - 1];
    psi[j - 1] = numerator / denominator;

    // Check for divergence - but be more lenient
    if (!isFinite(psi[j - 1]) || Math.abs(psi[j - 1]) > 1e15) {
      // Mark as divergent
      for (let k = 0; k <= j - 1; k++) {
        psi[k] = psi[j - 1];
      }
      return psi;
    }
  }

  // Now use parity to fill the left half
  if (parity === "symmetric") {
    // Symmetric: ψ(-x) = ψ(x)
    for (let i = 0; i < centerIdx; i++) {
      const mirrorIdx = 2 * centerIdx - i;
      if (mirrorIdx < N) {
        psi[i] = psi[mirrorIdx];
      }
    }
  } else {
    // Antisymmetric: ψ(-x) = -ψ(x)
    for (let i = 0; i < centerIdx; i++) {
      const mirrorIdx = 2 * centerIdx - i;
      if (mirrorIdx < N) {
        psi[i] = -psi[mirrorIdx];
      }
    }
  }

  return psi;
}

/**
 * Refine energy eigenvalue using bisection method.
 *
 * @param E1 - Lower energy bound (Joules)
 * @param E2 - Upper energy bound (Joules)
 * @param V - Potential energy array (Joules)
 * @param xGrid - Spatial grid (meters)
 * @param dx - Grid spacing (meters)
 * @param mass - Particle mass (kg)
 * @param tolerance - Energy tolerance (default 1e-10 Joules)
 * @returns Refined energy eigenvalue (Joules)
 */
export function refineEnergy(
  E1: number,
  E2: number,
  V: number[],
  xGrid: number[],
  dx: number,
  mass: number,
  tolerance = 1e-10,
): number {
  const N = xGrid.length;
  let Elow = E1;
  let Ehigh = E2;

  while (Ehigh - Elow > tolerance) {
    const Emid = (Elow + Ehigh) / 2;
    const psi = integrateNumerov(Emid, V, xGrid, dx, mass);
    const endValue = psi[N - 1];

    const psiLow = integrateNumerov(Elow, V, xGrid, dx, mass);
    const endValueLow = psiLow[N - 1];

    if (Math.sign(endValue) === Math.sign(endValueLow)) {
      Elow = Emid;
    } else {
      Ehigh = Emid;
    }
  }

  return (Elow + Ehigh) / 2;
}

/**
 * Normalize a wavefunction using trapezoidal rule.
 *
 * @param psi - Wavefunction array
 * @param dx - Grid spacing (meters)
 * @returns Normalized wavefunction
 */
export function normalizeWavefunction(psi: number[], dx: number): number[] {
  // Calculate ∫|ψ|² dx using trapezoidal rule
  let integral = 0;
  for (let i = 0; i < psi.length - 1; i++) {
    integral += (psi[i] * psi[i] + psi[i + 1] * psi[i + 1]) / 2;
  }
  integral *= dx;

  const normalization = Math.sqrt(integral);
  return psi.map((val) => val / normalization);
}

qppw.register("NumerovSolver", { solveNumerov });
