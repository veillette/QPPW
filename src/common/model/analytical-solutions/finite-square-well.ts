/**
 * Analytical solution for a finite square well.
 * V(x) = -V₀ for |x| < L/2, V(x) = 0 for |x| > L/2
 */

import QuantumConstants from "../QuantumConstants.js";
import { BoundStateResult, GridConfig } from "../PotentialFunction.js";

/**
 * Analytical solution for a finite square well.
 * V(x) = -V₀ for |x| < L/2, V(x) = 0 for |x| > L/2
 *
 * The energy eigenvalues are found by solving transcendental equations:
 * - Even parity: tan(ξ) = η/ξ
 * - Odd parity: -cot(ξ) = η/ξ
 * where ξ = (L/2)√(2mE/ℏ²) and η = (L/2)√(2m(V₀-E)/ℏ²)
 *
 * @param wellWidth - Width of the well (L) in meters
 * @param wellDepth - Depth of the well (V₀) in Joules (positive value)
 * @param mass - Particle mass in kg
 * @param numStates - Number of energy levels to calculate
 * @param gridConfig - Grid configuration for wavefunction evaluation
 * @returns Bound state results with energies and wavefunctions
 */
export function solveFiniteSquareWell(
  wellWidth: number,
  wellDepth: number,
  mass: number,
  numStates: number,
  gridConfig: GridConfig,
): BoundStateResult {
  const { HBAR } = QuantumConstants;
  const L = wellWidth;
  const V0 = wellDepth;
  const halfL = L / 2;

  // Dimensionless parameter: ξ₀ = (L/2)√(2mV₀/ℏ²)
  const xi0 = halfL * Math.sqrt(2 * mass * V0) / HBAR;

  // Maximum number of bound states (approximate)
  const maxStates = Math.floor(xi0 / (Math.PI / 2)) + 1;
  const actualNumStates = Math.min(numStates, maxStates);

  if (actualNumStates <= 0) {
    throw new Error("Finite square well too shallow to support bound states");
  }

  // Find energy eigenvalues by solving transcendental equations
  const energies: number[] = [];
  const parities: ("even" | "odd")[] = [];

  // Search for bound states alternating between even and odd parity
  let evenCount = 0;
  let oddCount = 0;

  for (let n = 0; n < actualNumStates; n++) {
    let xi: number;
    let parity: "even" | "odd";

    // Alternate between even (n=0,2,4,...) and odd (n=1,3,5,...)
    if (n % 2 === 0) {
      // Even parity state: solve tan(ξ) = √((ξ₀/ξ)² - 1)
      xi = findEvenParityState(xi0, evenCount);
      parity = "even";
      evenCount++;
    } else {
      // Odd parity state: solve -cot(ξ) = √((ξ₀/ξ)² - 1)
      xi = findOddParityState(xi0, oddCount);
      parity = "odd";
      oddCount++;
    }

    // Convert ξ to energy: E = (ℏξ/L/2)² / (2m) - V₀
    const E = (HBAR * HBAR * xi * xi) / (2 * mass * halfL * halfL) - V0;
    energies.push(E);
    parities.push(parity);
  }

  // Generate grid
  const numPoints = gridConfig.numPoints;
  const xGrid: number[] = [];
  const dx = (gridConfig.xMax - gridConfig.xMin) / (numPoints - 1);
  for (let i = 0; i < numPoints; i++) {
    xGrid.push(gridConfig.xMin + i * dx);
  }

  // Calculate wavefunctions
  const wavefunctions: number[][] = [];

  for (let n = 0; n < actualNumStates; n++) {
    const wavefunction: number[] = [];
    const E = energies[n];
    const parity = parities[n];

    // Wave numbers
    const k = Math.sqrt(2 * mass * (E + V0)) / HBAR;  // Inside the well
    const kappa = Math.sqrt(-2 * mass * E) / HBAR;     // Outside the well

    // Determine normalization constant
    // For even states: ψ(x) = A cos(kx) inside, B exp(-κ|x|) outside
    // For odd states: ψ(x) = A sin(kx) inside, B sign(x) exp(-κ|x|) outside

    let normalization: number;
    if (parity === "even") {
      // Match boundary conditions at x = L/2
      const cosVal = Math.cos(k * halfL);
      const B = cosVal * Math.exp(kappa * halfL);

      // Normalization integral
      const integral =
        2 * (halfL + Math.sin(2 * k * halfL) / (4 * k)) +
        2 * B * B / (2 * kappa);
      normalization = 1 / Math.sqrt(integral);
    } else {
      // Match boundary conditions at x = L/2
      const sinVal = Math.sin(k * halfL);
      const B = sinVal * Math.exp(kappa * halfL);

      // Normalization integral
      const integral =
        2 * (halfL - Math.sin(2 * k * halfL) / (4 * k)) +
        2 * B * B / (2 * kappa);
      normalization = 1 / Math.sqrt(integral);
    }

    // Evaluate wavefunction on grid
    for (const x of xGrid) {
      let value: number;

      if (Math.abs(x) <= halfL) {
        // Inside the well
        if (parity === "even") {
          value = normalization * Math.cos(k * x);
        } else {
          value = normalization * Math.sin(k * x);
        }
      } else {
        // Outside the well (exponentially decaying)
        if (parity === "even") {
          const cosVal = Math.cos(k * halfL);
          const B = normalization * cosVal * Math.exp(kappa * halfL);
          value = B * Math.exp(-kappa * Math.abs(x));
        } else {
          const sinVal = Math.sin(k * halfL);
          const B = normalization * sinVal * Math.exp(kappa * halfL);
          value = B * Math.sign(x) * Math.exp(-kappa * Math.abs(x));
        }
      }

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
 * Find even parity bound state by solving tan(ξ) = √((ξ₀/ξ)² - 1)
 */
function findEvenParityState(xi0: number, stateIndex: number): number {
  // Search in the interval [(n*π/2), ((n+1)*π/2)] where n is even
  const n = stateIndex * 2;
  const xiMin = n * Math.PI / 2 + 0.001;
  const xiMax = Math.min((n + 1) * Math.PI / 2 - 0.001, xi0);

  if (xiMin >= xiMax) {
    throw new Error(`No even parity state found for index ${stateIndex}`);
  }

  // Bisection method to solve tan(ξ) - η/ξ = 0
  const tolerance = 1e-10;
  let a = xiMin;
  let b = xiMax;

  const f = (xi: number) => {
    const eta = Math.sqrt(xi0 * xi0 - xi * xi);
    return Math.tan(xi) - eta / xi;
  };

  for (let iter = 0; iter < 100; iter++) {
    const c = (a + b) / 2;
    const fc = f(c);

    if (Math.abs(fc) < tolerance || (b - a) / 2 < tolerance) {
      return c;
    }

    if (f(a) * fc < 0) {
      b = c;
    } else {
      a = c;
    }
  }

  return (a + b) / 2;
}

/**
 * Find odd parity bound state by solving -cot(ξ) = √((ξ₀/ξ)² - 1)
 */
function findOddParityState(xi0: number, stateIndex: number): number {
  // Search in the interval [(n*π/2), ((n+1)*π/2)] where n is odd
  const n = stateIndex * 2 + 1;
  const xiMin = n * Math.PI / 2 + 0.001;
  const xiMax = Math.min((n + 1) * Math.PI / 2 - 0.001, xi0);

  if (xiMin >= xiMax) {
    throw new Error(`No odd parity state found for index ${stateIndex}`);
  }

  // Bisection method to solve -cot(ξ) - η/ξ = 0
  const tolerance = 1e-10;
  let a = xiMin;
  let b = xiMax;

  const f = (xi: number) => {
    const eta = Math.sqrt(xi0 * xi0 - xi * xi);
    return -1 / Math.tan(xi) - eta / xi;
  };

  for (let iter = 0; iter < 100; iter++) {
    const c = (a + b) / 2;
    const fc = f(c);

    if (Math.abs(fc) < tolerance || (b - a) / 2 < tolerance) {
      return c;
    }

    if (f(a) * fc < 0) {
      b = c;
    } else {
      a = c;
    }
  }

  return (a + b) / 2;
}
