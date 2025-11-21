/**
 * Analytical transcendental equations for symmetric double square well.
 *
 * Energy reference convention (matching the simulation):
 * - Left well:  [-L_outer, -L_inner] with V = 0
 * - Barrier:    [-L_inner, +L_inner] with V = V₀
 * - Right well: [+L_inner, +L_outer] with V = 0
 *
 * Where:
 *   L_inner = wellSeparation / 2
 *   L_outer = L_inner + wellWidth
 *
 * Due to symmetry, we solve on x ≥ 0 only:
 *
 * Region I:   [L_inner, L_outer]  → Inside right well (V = 0)
 * Region II:  [0, L_inner]        → Barrier region (V = V₀)
 * Region III: [L_outer, ∞)        → Outside (V = V₀, exponential decay)
 *
 * Wavefunctions have form:
 *
 * EVEN PARITY (ψ'(0) = 0):
 *   Region II (barrier):  ψ = A cosh(κx)
 *   Region I (well):      ψ = B cos(k(x - L_inner)) + C sin(k(x - L_inner))
 *   Region III (outside): ψ = D exp(-α(x - L_outer))
 *
 * ODD PARITY (ψ(0) = 0):
 *   Region II (barrier):  ψ = A sinh(κx)
 *   Region I (well):      ψ = B cos(k(x - L_inner)) + C sin(k(x - L_inner))
 *   Region III (outside): ψ = D exp(-α(x - L_outer))
 *
 * Where:
 *   k² = 2mE/ℏ²           (oscillatory in well, 0 < E < V₀)
 *   κ² = 2m(V₀ - E)/ℏ²    (exponential in barrier, E < V₀)
 *   α² = 2m(V₀ - E)/ℏ²    (same as κ for bound states)
 */

import QuantumConstants from "../QuantumConstants.js";
import { BoundStateResult, GridConfig } from "../PotentialFunction.js";

/**
 * Analytical solution for symmetric double square well potential.
 *
 * Solves the transcendental equations arising from boundary condition matching
 * to find bound state energies and wavefunctions.
 *
 * Energy convention: Wells at V=0, barrier at V=V₀ (positive).
 * Bound state energies are between 0 and V₀.
 *
 * @param wellWidth - Width of each well (L) in meters
 * @param wellDepth - Height of barrier relative to wells (V₀) in Joules (positive value)
 * @param wellSeparation - Barrier width (edge-to-edge well separation) in meters
 * @param mass - Particle mass in kg
 * @param numStates - Number of energy levels to calculate
 * @param gridConfig - Grid configuration for wavefunction evaluation
 * @returns Bound state results with energies and wavefunctions
 */
export function solveDoubleSquareWellAnalytical(
  wellWidth: number,
  wellDepth: number,
  wellSeparation: number,
  mass: number,
  numStates: number,
  gridConfig: GridConfig,
): BoundStateResult {

  // Geometry
  const Linner = wellSeparation / 2;
  const Louter = Linner + wellWidth;
  const V0 = wellDepth;

  // Search for extra states to ensure we don't miss any roots
  // For a symmetric double well, states MUST alternate in parity (even, odd, even, odd, ...)
  // We search for more states than needed to ensure robust root-finding
  const searchStates = numStates + 4;

  // Find even parity states
  const evenEnergies = findEvenParityDoubleWell(
    Linner, Louter, V0, mass, Math.ceil(searchStates / 2)
  );

  // Find odd parity states
  const oddEnergies = findOddParityDoubleWell(
    Linner, Louter, V0, mass, Math.ceil(searchStates / 2)
  );

  // Combine and sort states by energy
  const combinedStates: Array<{ energy: number; parity: "even" | "odd" }> = [];

  for (const E of evenEnergies) {
    combinedStates.push({ energy: E, parity: "even" });
  }
  for (const E of oddEnergies) {
    combinedStates.push({ energy: E, parity: "odd" });
  }

  // Sort by energy (ascending, lowest energy first)
  combinedStates.sort((a, b) => a.energy - b.energy);

  // For symmetric double well, states MUST alternate in parity
  // Ground state is always even, then odd, then even, etc.
  // Filter to ensure proper alternation
  const selectedStates: Array<{ energy: number; parity: "even" | "odd" }> = [];
  let expectedParity: "even" | "odd" = "even";  // Ground state is even

  for (const state of combinedStates) {
    if (state.parity === expectedParity) {
      selectedStates.push(state);
      // Alternate expected parity
      expectedParity = expectedParity === "even" ? "odd" : "even";

      // Stop when we have enough states
      if (selectedStates.length >= numStates) {
        break;
      }
    }
  }

  // If we didn't find enough alternating states, warn and use what we found
  if (selectedStates.length < numStates) {
    console.warn(
      `Warning: Only found ${selectedStates.length} properly alternating states out of ${numStates} requested. ` +
      `This may indicate numerical issues in root-finding or parameter regime near continuum threshold.`
    );
  }

  // Generate grid
  const numPoints = gridConfig.numPoints;
  const xGrid: number[] = [];
  const dx = (gridConfig.xMax - gridConfig.xMin) / (numPoints - 1);
  for (let i = 0; i < numPoints; i++) {
    xGrid.push(gridConfig.xMin + i * dx);
  }

  // Generate wavefunctions
  const wavefunctions: number[][] = [];
  const energies: number[] = [];

  for (const state of selectedStates) {
    const wf = computeDoubleWellWavefunction(
      state.energy, state.parity, Linner, Louter, V0, mass, xGrid
    );
    wavefunctions.push(wf);
    energies.push(state.energy);
  }

  return {
    energies,
    wavefunctions,
    xGrid,
    method: "analytical",
  };
}

/**
 * Validate that an energy is a true bound state, not a spurious root.
 *
 * Checks that:
 * 1. The wavefunction value at the outer boundary (D) is not too small
 * 2. The derivative matching condition is satisfied (small residual)
 *
 * Spurious roots can occur when D ≈ 0, causing a pole in the transcendental
 * equation that gets mistaken for a zero crossing.
 */
function isValidBoundState(
  E: number,
  Linner: number,
  L: number,
  V0: number,
  mass: number,
  parity: "even" | "odd"
): boolean {
  const { HBAR } = QuantumConstants;

  // Wave numbers
  const k = Math.sqrt(2 * mass * E) / HBAR;
  const kappa = Math.sqrt(2 * mass * (V0 - E)) / HBAR;
  const alpha = kappa;

  // Calculate coefficients
  let B: number, C: number;
  if (parity === "even") {
    B = Math.cosh(kappa * Linner);
    C = (kappa / k) * Math.sinh(kappa * Linner);
  } else {
    B = Math.sinh(kappa * Linner);
    C = (kappa / k) * Math.cosh(kappa * Linner);
  }

  // Calculate D (wavefunction value at outer boundary)
  const D = B * Math.cos(k * L) + C * Math.sin(k * L);

  // Reject if D is too small (near-node at boundary)
  // Use relative threshold compared to B (typical barrier value)
  // A very small D indicates a spurious pole in the transcendental equation
  const relativeDThreshold = 0.001;  // 0.1% of typical barrier amplitude
  if (Math.abs(D) < Math.abs(B) * relativeDThreshold) {
    return false;
  }

  // Verify derivative matching condition
  const numerator = -k * B * Math.sin(k * L) + k * C * Math.cos(k * L);
  const denominator = D;

  const lhs = numerator / denominator;  // ψ'/ψ from inside
  const rhs = -alpha;  // ψ'/ψ from outside (should match)

  // Check that matching condition is satisfied
  // Allow for some numerical error in the root-finding
  const relativeError = Math.abs((lhs - rhs) / rhs);
  const matchingTolerance = 0.05;  // 5% relative error

  if (relativeError > matchingTolerance) {
    return false;
  }

  return true;
}

/**
 * Find even parity bound states: ψ'(0) = 0
 *
 * Matching conditions give us a transcendental equation in E.
 *
 * From boundary matching:
 * - At x = L_inner: A cosh(κL_i) = B, Aκ sinh(κL_i) = -kC
 * - At x = L_outer: B cos(kL) + C sin(kL) = D, -k[B sin(kL) - C cos(kL)] = -αD
 *
 * Where L = L_outer - L_inner = wellWidth
 */
function findEvenParityDoubleWell(
  Linner: number,
  Louter: number,
  V0: number,
  mass: number,
  maxStates: number,
): number[] {

  const { HBAR } = QuantumConstants;
  const energies: number[] = [];
  const L = Louter - Linner; // wellWidth

  /**
   * Transcendental equation for even parity states.
   * Returns f(E) = 0 when E is an eigenvalue.
   */
  const transcendentalEquation = (E: number): number => {
    // Must be bound state (below barrier top)
    if (E >= V0) return Infinity;

    // Must be above well bottom
    if (E <= 0) {
      return Infinity;
    }

    // Wave numbers
    const k = Math.sqrt(2 * mass * E) / HBAR;           // In wells
    const kappa = Math.sqrt(2 * mass * (V0 - E)) / HBAR; // In barrier
    const alpha = kappa; // Same for bound states (outside also at V₀)

    // Hyperbolic functions at barrier/well interface
    const coshKL = Math.cosh(kappa * Linner);
    const sinhKL = Math.sinh(kappa * Linner);

    // From matching conditions:
    // B = cosh(κL_i), C = (κ/k)sinh(κL_i)
    // Substituting into well/outside matching:
    // ψ'_well = -Bk sin(kL) + Ck cos(kL) = -α D
    // ψ_well = B cos(kL) + C sin(kL) = D

    const numerator = -k * coshKL * Math.sin(k * L) + k * (kappa/k) * sinhKL * Math.cos(k * L);
    const denominator = coshKL * Math.cos(k * L) + (kappa/k) * sinhKL * Math.sin(k * L);

    if (Math.abs(denominator) < 1e-15) return Infinity;

    const lhs = numerator / denominator;
    const rhs = alpha;

    // Matching condition: ψ'_well / ψ_well = ψ'_outside / ψ_outside = -α
    return lhs + rhs;
  };

  // Search for roots in energy range
  const Emin = 0;   // Well bottom
  const Emax = V0;  // Barrier top (continuum threshold)

  // Use systematic search with bisection
  // Increased resolution to ensure we don't miss any roots
  const numSearchPoints = 1500;
  const dE = (Emax - Emin) / numSearchPoints;

  for (let i = 0; i < numSearchPoints - 1; i++) {
    const E1 = Emin + i * dE;
    const E2 = Emin + (i + 1) * dE;

    const f1 = transcendentalEquation(E1);
    const f2 = transcendentalEquation(E2);

    // Check for sign change (root exists)
    if (f1 * f2 < 0 && isFinite(f1) && isFinite(f2)) {
      const root = solveBisection(transcendentalEquation, E1, E2, 1e-12, 100);
      if (root !== null) {
        // Validate this is a true bound state, not a spurious root
        if (isValidBoundState(root, Linner, L, V0, mass, "even")) {
          energies.push(root);
          if (energies.length >= maxStates) break;
        }
      }
    }
  }

  return energies;
}

/**
 * Find odd parity bound states: ψ(0) = 0
 *
 * Similar derivation but with sinh instead of cosh in barrier.
 *
 * From boundary matching:
 * - At x = L_inner: A sinh(κL_i) = B, Aκ cosh(κL_i) = -kC
 * - At x = L_outer: B cos(kL) + C sin(kL) = D, -k[B sin(kL) - C cos(kL)] = -αD
 */
function findOddParityDoubleWell(
  Linner: number,
  Louter: number,
  V0: number,
  mass: number,
  maxStates: number,
): number[] {

  const { HBAR } = QuantumConstants;
  const energies: number[] = [];
  const L = Louter - Linner; // wellWidth

  /**
   * Transcendental equation for odd parity states.
   * Returns f(E) = 0 when E is an eigenvalue.
   */
  const transcendentalEquation = (E: number): number => {
    if (E >= V0) return Infinity;

    if (E <= 0) return Infinity;

    const k = Math.sqrt(2 * mass * E) / HBAR;
    const kappa = Math.sqrt(2 * mass * (V0 - E)) / HBAR;
    const alpha = kappa;

    // Hyperbolic functions at barrier/well interface
    const sinhKL = Math.sinh(kappa * Linner);
    const coshKL = Math.cosh(kappa * Linner);

    // From matching conditions:
    // B = sinh(κL_i), C = (κ/k)cosh(κL_i)
    // ψ'_well = -Bk sin(kL) + Ck cos(kL) = -α D
    // ψ_well = B cos(kL) + C sin(kL) = D

    const numerator = -k * sinhKL * Math.sin(k * L) + k * (kappa/k) * coshKL * Math.cos(k * L);
    const denominator = sinhKL * Math.cos(k * L) + (kappa/k) * coshKL * Math.sin(k * L);

    if (Math.abs(denominator) < 1e-15) return Infinity;

    const lhs = numerator / denominator;
    const rhs = alpha;

    // Matching condition: ψ'_well / ψ_well = ψ'_outside / ψ_outside = -α
    return lhs + rhs;
  };

  // Search for roots
  const Emin = 0;
  const Emax = V0;
  // Increased resolution to ensure we don't miss any roots
  const numSearchPoints = 1500;
  const dE = (Emax - Emin) / numSearchPoints;

  for (let i = 0; i < numSearchPoints - 1; i++) {
    const E1 = Emin + i * dE;
    const E2 = Emin + (i + 1) * dE;

    const f1 = transcendentalEquation(E1);
    const f2 = transcendentalEquation(E2);

    if (f1 * f2 < 0 && isFinite(f1) && isFinite(f2)) {
      const root = solveBisection(transcendentalEquation, E1, E2, 1e-12, 100);
      if (root !== null) {
        // Validate this is a true bound state, not a spurious root
        if (isValidBoundState(root, Linner, L, V0, mass, "odd")) {
          energies.push(root);
          if (energies.length >= maxStates) break;
        }
      }
    }
  }

  return energies;
}

/**
 * Compute wavefunction for a given energy and parity on the provided grid.
 *
 * Constructs the full wavefunction by matching boundary conditions and normalizing.
 *
 * Energy convention: Wells at V=0, barrier at V=V₀.
 *
 * @param E - Energy eigenvalue (between 0 and V₀)
 * @param parity - "even" or "odd" parity
 * @param Linner - Inner boundary (half the barrier width)
 * @param Louter - Outer boundary (Linner + wellWidth)
 * @param V0 - Barrier height (positive)
 * @param mass - Particle mass
 * @param xGrid - Position grid
 * @returns Normalized wavefunction values
 */
function computeDoubleWellWavefunction(
  E: number,
  parity: "even" | "odd",
  Linner: number,
  Louter: number,
  V0: number,
  mass: number,
  xGrid: number[],
): number[] {

  const { HBAR } = QuantumConstants;
  const k = Math.sqrt(2 * mass * E) / HBAR;
  const kappa = Math.sqrt(2 * mass * (V0 - E)) / HBAR;
  const alpha = kappa;

  const L = Louter - Linner; // wellWidth

  // Solve for coefficients A, B, C, D from matching conditions
  // Set A = 1 for convenience, then normalize at end

  const A = 1.0;
  let B: number, C: number;

  if (parity === "even") {
    // A cosh(κ L_inner) = B
    B = Math.cosh(kappa * Linner);
    // A κ sinh(κ L_inner) = k C
    C = (kappa / k) * Math.sinh(kappa * Linner);
  } else {
    // A sinh(κ L_inner) = B
    B = Math.sinh(kappa * Linner);
    // A κ cosh(κ L_inner) = k C
    C = (kappa / k) * Math.cosh(kappa * Linner);
  }

  // From well/outside matching:
  // B cos(kL) + C sin(kL) = D
  const D = B * Math.cos(k * L) + C * Math.sin(k * L);

  // Evaluate on grid (not yet normalized)
  const psi: number[] = [];

  for (const x of xGrid) {
    let value: number;

    if (x >= 0) {
      // Right half: evaluate directly
      if (x < Linner) {
        // Barrier region (0 to L_inner)
        if (parity === "even") {
          value = A * Math.cosh(kappa * x);
        } else {
          value = A * Math.sinh(kappa * x);
        }
      } else if (x < Louter) {
        // Inside right well (L_inner to L_outer)
        const xShifted = x - Linner;
        value = B * Math.cos(k * xShifted) + C * Math.sin(k * xShifted);
      } else {
        // Outside right well (x > L_outer)
        value = D * Math.exp(-alpha * (x - Louter));
      }
    } else {
      // Left half: use symmetry/antisymmetry
      const absX = -x;

      if (absX < Linner) {
        // Barrier region
        if (parity === "even") {
          value = A * Math.cosh(kappa * absX);
        } else {
          value = -A * Math.sinh(kappa * absX);  // Antisymmetric
        }
      } else if (absX < Louter) {
        // Inside left well
        const xShifted = absX - Linner;
        const valueRight = B * Math.cos(k * xShifted) + C * Math.sin(k * xShifted);
        value = (parity === "even") ? valueRight : -valueRight;
      } else {
        // Outside left well
        const valueRight = D * Math.exp(-alpha * (absX - Louter));
        value = (parity === "even") ? valueRight : -valueRight;
      }
    }

    psi.push(value);
  }

  // Normalize using trapezoidal rule
  const dx = xGrid[1] - xGrid[0];
  const normSq = psi.reduce((sum, val) => sum + val * val, 0) * dx;
  const norm = Math.sqrt(normSq);

  return psi.map(val => val / norm);
}

/**
 * Solve equation f(x) = 0 using bisection method.
 * Robust root-finding algorithm for continuous functions.
 *
 * @param f - Function to find root of
 * @param xMin - Lower bound of search interval
 * @param xMax - Upper bound of search interval
 * @param tolerance - Convergence tolerance
 * @param maxIterations - Maximum number of iterations
 * @returns Root if found, null otherwise
 */
function solveBisection(
  f: (x: number) => number,
  xMin: number,
  xMax: number,
  tolerance: number,
  maxIterations: number
): number | null {
  let a = xMin;
  let b = xMax;

  // Check if function values have opposite signs
  const fa = f(a);
  const fb = f(b);

  if (fa * fb > 0) {
    // No sign change, bisection won't work
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

  // Return best approximation even if not fully converged
  const c = (a + b) / 2;
  if (Math.abs(f(c)) < tolerance * 10) {
    return c;
  }

  return null;
}
