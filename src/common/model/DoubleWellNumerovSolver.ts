/**
 * Specialized Numerov solver for double square well potentials.
 *
 * This solver uses the key insight that double well energy levels can be estimated
 * from single well energies, then refined using the Numerov shooting method.
 *
 * Physical principle:
 * - When two identical wells are separated, each single-well energy level E_n
 *   splits into a symmetric (E_n^+) and antisymmetric (E_n^-) pair
 * - The splitting depends on barrier width and tunneling probability
 * - Symmetric states: ψ(-x) = ψ(x), lower energy
 * - Antisymmetric states: ψ(-x) = -ψ(x), higher energy
 */

import QuantumConstants from "./QuantumConstants.js";
import { BoundStateResult, GridConfig, PotentialFunction } from "./PotentialFunction.js";
import { integrateNumerov, normalizeWavefunction } from "./NumerovSolver.js";
import qppw from "../../QPPWNamespace.js";

/**
 * Solve the double square well using Numerov method with intelligent energy search.
 *
 * @param wellWidth - Width of each well in meters
 * @param wellDepth - Depth of each well in Joules (positive value)
 * @param wellSeparation - Center-to-center separation between wells in meters
 * @param mass - Particle mass in kg
 * @param numStates - Number of bound states to find
 * @param gridConfig - Grid configuration for wavefunction evaluation
 * @returns Bound state results with energies and wavefunctions
 */
export function solveDoubleWellNumerov(
  wellWidth: number,
  wellDepth: number,
  wellSeparation: number,
  mass: number,
  numStates: number,
  gridConfig: GridConfig,
): BoundStateResult {

  // Create the double square well potential function
  // Wells centered at ±(separation/2 + wellWidth/2)
  const leftWellCenter = -(wellSeparation / 2 + wellWidth / 2);
  const rightWellCenter = wellSeparation / 2 + wellWidth / 2;
  const halfWidth = wellWidth / 2;

  const potential: PotentialFunction = (x: number) => {
    const inLeftWell = Math.abs(x - leftWellCenter) <= halfWidth;
    const inRightWell = Math.abs(x - rightWellCenter) <= halfWidth;

    if (inLeftWell || inRightWell) {
      return wellDepth; // Energy zero at bottom of wells
    } else {
      return 0; // Barrier at zero energy
    }
  };

  // Step 1: Estimate single well energies
  const singleWellEnergies = estimateSingleWellEnergies(
    wellWidth,
    wellDepth,
    mass,
    numStates,
  );

  if (singleWellEnergies.length === 0) {
    console.warn("No bound states found in single well approximation");
    return {
      energies: [],
      wavefunctions: [],
      xGrid: [],
      method: "numerov",
    };
  }

  // Generate spatial grid
  const { xMin, xMax, numPoints } = gridConfig;
  const dx = (xMax - xMin) / (numPoints - 1);
  const xGrid: number[] = [];
  for (let i = 0; i < numPoints; i++) {
    xGrid.push(xMin + i * dx);
  }

  // Evaluate potential on grid
  const V = xGrid.map(potential);

  // Step 2: Find double well energies by searching near single well estimates
  const energies: number[] = [];
  const wavefunctions: number[][] = [];
  const parities: ("symmetric" | "antisymmetric")[] = [];

  // Estimate level splitting based on barrier tunneling
  const barrierWidth = wellSeparation; // Approximate barrier width

  for (let i = 0; i < singleWellEnergies.length && energies.length < numStates; i++) {
    const E_single = singleWellEnergies[i];

    // Skip if energy is above the barrier (not bound)
    if (E_single >= 0) {
      continue;
    }

    // Estimate energy splitting using WKB approximation
    const splitting = estimateEnergySplitting(
      E_single,
      wellDepth,
      barrierWidth,
      mass,
    );

    // Search for symmetric state (lower energy)
    const E_sym_estimate = E_single - splitting / 2;
    const E_antisym_estimate = E_single + splitting / 2;

    // Define search windows around estimates
    const searchWindow = Math.max(splitting * 2, wellDepth * 0.01); // At least 1% of well depth

    // Search for symmetric state
    const E_sym = findEnergyNearEstimate(
      E_sym_estimate,
      searchWindow,
      "symmetric",
      V,
      xGrid,
      dx,
      mass,
      leftWellCenter,
      rightWellCenter,
    );

    if (E_sym !== null && energies.length < numStates) {
      energies.push(E_sym);
      const psi_sym = integrateNumerov(E_sym, V, xGrid, dx, mass);
      const normalized_psi_sym = normalizeWavefunction(psi_sym, dx);
      wavefunctions.push(normalized_psi_sym);
      parities.push("symmetric");
    }

    // Search for antisymmetric state
    const E_antisym = findEnergyNearEstimate(
      E_antisym_estimate,
      searchWindow,
      "antisymmetric",
      V,
      xGrid,
      dx,
      mass,
      leftWellCenter,
      rightWellCenter,
    );

    if (E_antisym !== null && energies.length < numStates) {
      energies.push(E_antisym);
      const psi_antisym = integrateNumerov(E_antisym, V, xGrid, dx, mass);
      const normalized_psi_antisym = normalizeWavefunction(psi_antisym, dx);
      wavefunctions.push(normalized_psi_antisym);
      parities.push("antisymmetric");
    }
  }

  return {
    energies,
    wavefunctions,
    xGrid,
    method: "numerov",
  };
}

/**
 * Estimate single well energy levels using finite square well formulas.
 *
 * For a finite square well: V(x) = -V₀ for |x| < L/2, V(x) = 0 for |x| > L/2
 * Energy eigenvalues satisfy transcendental equations:
 * - Even: tan(ξ) = √((ξ₀/ξ)² - 1)
 * - Odd: -cot(ξ) = √((ξ₀/ξ)² - 1)
 * where ξ = (L/2)√(2m(E+V₀)/ℏ²) and ξ₀ = (L/2)√(2mV₀/ℏ²)
 *
 * @returns Array of estimated energies in Joules (negative values, measured from V=0)
 */
function estimateSingleWellEnergies(
  wellWidth: number,
  wellDepth: number,
  mass: number,
  maxStates: number,
): number[] {
  const { HBAR } = QuantumConstants;
  const L = wellWidth;
  const V0 = wellDepth;
  const halfL = L / 2;

  // Dimensionless parameter: ξ₀ = (L/2)√(2mV₀/ℏ²)
  const xi0 = halfL * Math.sqrt(2 * mass * V0) / HBAR;

  // Maximum number of bound states
  const maxPossibleStates = Math.floor(xi0 / (Math.PI / 2)) + 1;
  const numStates = Math.min(maxStates, maxPossibleStates);

  const energies: number[] = [];

  if (numStates <= 0) {
    return energies;
  }

  // Find energy levels by solving transcendental equations
  for (let n = 0; n < numStates; n++) {
    const isEven = (n % 2 === 0);
    const stateIndex = Math.floor(n / 2);

    let xi: number | null = null;

    if (isEven) {
      // Even parity: solve tan(ξ) = √((ξ₀/ξ)² - 1)
      xi = findEvenParityXi(xi0, stateIndex);
    } else {
      // Odd parity: solve -cot(ξ) = √((ξ₀/ξ)² - 1)
      xi = findOddParityXi(xi0, stateIndex);
    }

    if (xi !== null) {
      // Convert ξ to energy: E = (ℏξ/L/2)² / (2m) - V₀
      const E = (HBAR * HBAR * xi * xi) / (2 * mass * halfL * halfL) - V0;
      energies.push(E);
    }
  }

  return energies;
}

/**
 * Find even parity eigenvalue ξ by solving: tan(ξ) = √((ξ₀/ξ)² - 1)
 */
function findEvenParityXi(xi0: number, stateIndex: number): number | null {
  const n = stateIndex * 2;
  const xiMin = n * Math.PI / 2 + 0.01;
  const xiMax = Math.min((n + 1) * Math.PI / 2 - 0.01, xi0);

  if (xiMin >= xiMax) {
    return null;
  }

  const f = (xi: number) => {
    const eta = Math.sqrt(xi0 * xi0 - xi * xi);
    return Math.tan(xi) - eta / xi;
  };

  return solveBisection(f, xiMin, xiMax, 1e-8, 100);
}

/**
 * Find odd parity eigenvalue ξ by solving: -cot(ξ) = √((ξ₀/ξ)² - 1)
 */
function findOddParityXi(xi0: number, stateIndex: number): number | null {
  const n = stateIndex * 2 + 1;
  const xiMin = n * Math.PI / 2 + 0.01;
  const xiMax = Math.min((n + 1) * Math.PI / 2 - 0.01, xi0);

  if (xiMin >= xiMax) {
    return null;
  }

  const f = (xi: number) => {
    const eta = Math.sqrt(xi0 * xi0 - xi * xi);
    return -1 / Math.tan(xi) - eta / xi;
  };

  return solveBisection(f, xiMin, xiMax, 1e-8, 100);
}

/**
 * Estimate energy splitting between symmetric and antisymmetric states
 * using WKB approximation for tunneling.
 *
 * ΔE ≈ (ℏω) * exp(-∫√(2m(V-E))/ℏ dx)
 * where the integral is over the barrier region.
 *
 * @param energy - Single well energy estimate (Joules, negative)
 * @param wellDepth - Well depth (Joules, positive)
 * @param barrierWidth - Barrier width (meters)
 * @param mass - Particle mass (kg)
 * @returns Estimated energy splitting (Joules, positive)
 */
function estimateEnergySplitting(
  energy: number,
  wellDepth: number,
  barrierWidth: number,
  mass: number,
): number {
  const { HBAR } = QuantumConstants;

  // Energy relative to barrier top (negative in classical forbidden region)
  const E_rel = energy; // E is already negative (below V=0 barrier)
  const V_barrier = 0; // Barrier at zero

  if (E_rel >= V_barrier) {
    // Classical region - large splitting
    return wellDepth * 0.1;
  }

  // WKB tunneling factor
  const kappa = Math.sqrt(2 * mass * (V_barrier - E_rel)) / HBAR;
  const tunneling = Math.exp(-kappa * barrierWidth);

  // Characteristic energy scale
  const E_scale = Math.abs(energy);

  // Splitting proportional to tunneling probability and energy scale
  // Factor of 2 accounts for symmetric/antisymmetric splitting
  const splitting = 2 * E_scale * tunneling;

  // Ensure splitting is reasonable (not too small or too large)
  const minSplitting = wellDepth * 1e-6;
  const maxSplitting = wellDepth * 0.5;

  return Math.max(minSplitting, Math.min(splitting, maxSplitting));
}

/**
 * Find energy eigenvalue near an estimate by scanning and refining.
 *
 * @param energyEstimate - Initial energy estimate (Joules)
 * @param searchWindow - Search range around estimate (Joules)
 * @param parity - "symmetric" or "antisymmetric"
 * @param V - Potential energy array (Joules)
 * @param xGrid - Spatial grid (meters)
 * @param dx - Grid spacing (meters)
 * @param mass - Particle mass (kg)
 * @param leftCenter - Left well center position (meters)
 * @param rightCenter - Right well center position (meters) - reserved for future use
 * @returns Refined energy eigenvalue or null if not found
 */
function findEnergyNearEstimate(
  energyEstimate: number,
  searchWindow: number,
  parity: "symmetric" | "antisymmetric",
  V: number[],
  xGrid: number[],
  dx: number,
  mass: number,
  // eslint-disable-next-line @typescript-eslint/no-unused-vars
  _leftCenter: number,
  // eslint-disable-next-line @typescript-eslint/no-unused-vars
  _rightCenter: number,
): number | null {
  const energyMin = energyEstimate - searchWindow;
  const energyMax = energyEstimate + searchWindow;

  // Ensure we're searching for bound states (E < 0)
  const actualEnergyMax = Math.min(energyMax, -1e-20); // Must be below barrier

  if (energyMin >= actualEnergyMax) {
    return null;
  }

  // Coarse scan to find sign changes
  const numScanPoints = 100;
  const scanStep = (actualEnergyMax - energyMin) / numScanPoints;

  let prevSign = 0;
  let prevEnergy = energyMin;

  for (let i = 0; i <= numScanPoints; i++) {
    const E = energyMin + i * scanStep;
    const psi = integrateNumerov(E, V, xGrid, dx, mass);

    // Check boundary condition: wavefunction should decay at boundaries
    const endValue = psi[xGrid.length - 1];

    // For both symmetric and antisymmetric states, use boundary value
    // The sign change in the boundary value indicates a bound state
    const meritFunction = endValue;

    const currentSign = Math.sign(meritFunction);

    // Detect sign change
    if (prevSign !== 0 && currentSign !== 0 && currentSign !== prevSign) {
      // Refine using bisection
      const refinedEnergy = refineEnergy(
        prevEnergy,
        E,
        V,
        xGrid,
        dx,
        mass,
        parity,
      );

      if (refinedEnergy !== null) {
        return refinedEnergy;
      }
    }

    prevSign = currentSign;
    prevEnergy = E;
  }

  return null;
}

/**
 * Refine energy eigenvalue using bisection method.
 */
function refineEnergy(
  E1: number,
  E2: number,
  V: number[],
  xGrid: number[],
  dx: number,
  mass: number,
  parity: "symmetric" | "antisymmetric",
  tolerance = 1e-10,
): number | null {
  const N = xGrid.length;
  let Elow = E1;
  let Ehigh = E2;

  let iterations = 0;
  const maxIterations = 100;

  while (Ehigh - Elow > tolerance && iterations < maxIterations) {
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

    iterations++;
  }

  // Verify the solution has the correct parity
  const finalE = (Elow + Ehigh) / 2;
  const finalPsi = integrateNumerov(finalE, V, xGrid, dx, mass);

  if (!verifyParity(finalPsi, xGrid, parity)) {
    return null; // Wrong parity, reject
  }

  return finalE;
}

/**
 * Verify that a wavefunction has the expected parity.
 */
function verifyParity(
  psi: number[],
  xGrid: number[],
  expectedParity: "symmetric" | "antisymmetric",
): boolean {
  const N = xGrid.length;
  const centerIdx = Math.floor(N / 2);

  // Check symmetry around center
  let symmetryError = 0;
  let antisymmetryError = 0;
  const numCheckPoints = Math.min(50, Math.floor(N / 4));

  for (let i = 1; i <= numCheckPoints; i++) {
    const leftIdx = centerIdx - i;
    const rightIdx = centerIdx + i;

    if (leftIdx >= 0 && rightIdx < N) {
      const diff = Math.abs(psi[leftIdx] - psi[rightIdx]);
      const sum = Math.abs(psi[leftIdx] + psi[rightIdx]);

      symmetryError += diff;
      antisymmetryError += sum;
    }
  }

  // Normalize errors
  const maxAbs = Math.max(...psi.map(Math.abs));
  symmetryError /= (numCheckPoints * maxAbs + 1e-10);
  antisymmetryError /= (numCheckPoints * maxAbs + 1e-10);

  // Determine parity (threshold: 0.1)
  const isSymmetric = symmetryError < antisymmetryError && symmetryError < 0.1;
  const isAntisymmetric = antisymmetryError < symmetryError && antisymmetryError < 0.1;

  if (expectedParity === "symmetric") {
    return isSymmetric;
  } else {
    return isAntisymmetric;
  }
}

/**
 * Solve equation f(x) = 0 using bisection method.
 */
function solveBisection(
  f: (x: number) => number,
  xMin: number,
  xMax: number,
  tolerance: number,
  maxIterations: number,
): number | null {
  let a = xMin;
  let b = xMax;

  const fa = f(a);
  const fb = f(b);

  // Check if function values have opposite signs
  if (fa * fb > 0) {
    return null;
  }

  for (let iter = 0; iter < maxIterations; iter++) {
    const c = (a + b) / 2;
    const fc = f(c);

    // Check convergence
    if (Math.abs(fc) < tolerance || (b - a) / 2 < tolerance) {
      return c;
    }

    // Update interval
    if (fa * fc < 0) {
      b = c;
    } else {
      a = c;
    }
  }

  return (a + b) / 2;
}

qppw.register("DoubleWellNumerovSolver", { solveDoubleWellNumerov });
