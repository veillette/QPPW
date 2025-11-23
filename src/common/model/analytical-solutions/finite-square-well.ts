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
 * This implementation uses multiple numerical methods with fallbacks:
 * 1. Bisection method (robust but slower)
 * 2. Newton-Raphson method (fast but requires good initial guess)
 * 3. Secant method (no derivatives needed)
 * 4. Lima's approximation (2020) for difficult cases
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
  const xi0 = (halfL * Math.sqrt(2 * mass * V0)) / HBAR;

  // Maximum number of bound states (approximate)
  const maxStates = Math.floor(xi0 / (Math.PI / 2)) + 1;
  const actualNumStates = Math.min(numStates, maxStates);

  if (actualNumStates <= 0) {
    // Return empty result instead of throwing - well is too shallow
    console.warn("Finite square well too shallow to support bound states");
    return {
      energies: [],
      wavefunctions: [],
      xGrid: [],
      method: "analytical",
    };
  }

  // Find energy eigenvalues by solving transcendental equations
  const energies: number[] = [];
  const parities: ("even" | "odd")[] = [];

  // Search for bound states alternating between even and odd parity
  let evenCount = 0;
  let oddCount = 0;

  for (let n = 0; n < actualNumStates; n++) {
    let xi: number | null = null;
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

    // Skip if state could not be found (fallback returned null)
    if (xi === null) {
      console.warn(`Could not find ${parity} parity state ${n} for xi0=${xi0}`);
      continue;
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
    const k = Math.sqrt(2 * mass * (E + V0)) / HBAR; // Inside the well
    const kappa = Math.sqrt(-2 * mass * E) / HBAR; // Outside the well

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
        (2 * B * B) / (2 * kappa);
      normalization = 1 / Math.sqrt(integral);
    } else {
      // Match boundary conditions at x = L/2
      const sinVal = Math.sin(k * halfL);
      const B = sinVal * Math.exp(kappa * halfL);

      // Normalization integral
      const integral =
        2 * (halfL - Math.sin(2 * k * halfL) / (4 * k)) +
        (2 * B * B) / (2 * kappa);
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
 * Uses multiple numerical methods with fallbacks for robustness.
 *
 * @returns The ξ value for the bound state, or null if not found
 */
function findEvenParityState(xi0: number, stateIndex: number): number | null {
  // Search in the interval [(n*π/2), ((n+1)*π/2)] where n is even
  const n = stateIndex * 2;
  const xiMin = (n * Math.PI) / 2 + 0.001;
  const xiMax = Math.min(((n + 1) * Math.PI) / 2 - 0.001, xi0);

  if (xiMin >= xiMax) {
    // Well is not deep enough for this state
    return null;
  }

  const f = (xi: number) => {
    const eta = Math.sqrt(xi0 * xi0 - xi * xi);
    return Math.tan(xi) - eta / xi;
  };

  const fDerivative = (xi: number) => {
    const eta = Math.sqrt(xi0 * xi0 - xi * xi);
    const tanXi = Math.tan(xi);
    const sec2Xi = 1 + tanXi * tanXi; // sec²(ξ) = 1 + tan²(ξ)
    return sec2Xi + eta / (xi * xi) + xi / eta;
  };

  // Method 1: Try bisection method (most robust)
  const bisectionResult = solveBisection(f, xiMin, xiMax, 1e-10, 100);
  if (bisectionResult !== null) {
    return bisectionResult;
  }

  // Method 2: Try Newton-Raphson with Lima's initial guess
  const limaGuess = getLimaApproximation(xi0, stateIndex, "even");
  if (limaGuess !== null && limaGuess > xiMin && limaGuess < xiMax) {
    const newtonResult = solveNewtonRaphson(
      f,
      fDerivative,
      limaGuess,
      1e-10,
      50,
    );
    if (
      newtonResult !== null &&
      newtonResult >= xiMin &&
      newtonResult <= xiMax
    ) {
      return newtonResult;
    }
  }

  // Method 3: Try secant method with midpoint initial guess
  const midpoint = (xiMin + xiMax) / 2;
  const secantResult = solveSecant(f, xiMin, midpoint, 1e-10, 50);
  if (secantResult !== null && secantResult >= xiMin && secantResult <= xiMax) {
    return secantResult;
  }

  // Method 4: Try Lima's approximation directly (if available)
  if (limaGuess !== null && limaGuess >= xiMin && limaGuess <= xiMax) {
    // Verify it's a reasonable solution
    if (Math.abs(f(limaGuess)) < 1e-6) {
      return limaGuess;
    }
  }

  // All methods failed
  console.warn(`All methods failed to find even parity state ${stateIndex}`);
  return null;
}

/**
 * Find odd parity bound state by solving -cot(ξ) = √((ξ₀/ξ)² - 1)
 * Uses multiple numerical methods with fallbacks for robustness.
 *
 * @returns The ξ value for the bound state, or null if not found
 */
function findOddParityState(xi0: number, stateIndex: number): number | null {
  // Search in the interval [(n*π/2), ((n+1)*π/2)] where n is odd
  const n = stateIndex * 2 + 1;
  const xiMin = (n * Math.PI) / 2 + 0.001;
  const xiMax = Math.min(((n + 1) * Math.PI) / 2 - 0.001, xi0);

  if (xiMin >= xiMax) {
    // Well is not deep enough for this state
    return null;
  }

  const f = (xi: number) => {
    const eta = Math.sqrt(xi0 * xi0 - xi * xi);
    return -1 / Math.tan(xi) - eta / xi;
  };

  const fDerivative = (xi: number) => {
    const eta = Math.sqrt(xi0 * xi0 - xi * xi);
    const sinXi = Math.sin(xi);
    const csc2Xi = 1 / (sinXi * sinXi); // csc²(ξ) = 1/sin²(ξ)
    return csc2Xi + eta / (xi * xi) + xi / eta;
  };

  // Method 1: Try bisection method (most robust)
  const bisectionResult = solveBisection(f, xiMin, xiMax, 1e-10, 100);
  if (bisectionResult !== null) {
    return bisectionResult;
  }

  // Method 2: Try Newton-Raphson with Lima's initial guess
  const limaGuess = getLimaApproximation(xi0, stateIndex, "odd");
  if (limaGuess !== null && limaGuess > xiMin && limaGuess < xiMax) {
    const newtonResult = solveNewtonRaphson(
      f,
      fDerivative,
      limaGuess,
      1e-10,
      50,
    );
    if (
      newtonResult !== null &&
      newtonResult >= xiMin &&
      newtonResult <= xiMax
    ) {
      return newtonResult;
    }
  }

  // Method 3: Try secant method with midpoint initial guess
  const midpoint = (xiMin + xiMax) / 2;
  const secantResult = solveSecant(f, xiMin, midpoint, 1e-10, 50);
  if (secantResult !== null && secantResult >= xiMin && secantResult <= xiMax) {
    return secantResult;
  }

  // Method 4: Try Lima's approximation directly (if available)
  if (limaGuess !== null && limaGuess >= xiMin && limaGuess <= xiMax) {
    // Verify it's a reasonable solution
    if (Math.abs(f(limaGuess)) < 1e-6) {
      return limaGuess;
    }
  }

  // All methods failed
  console.warn(`All methods failed to find odd parity state ${stateIndex}`);
  return null;
}

/**
 * Solve equation f(x) = 0 using bisection method.
 * Improved version with better convergence checking.
 *
 * @returns Root if found, null otherwise
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

  // Check if function values have opposite signs (necessary for bisection)
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
      // fb = fc; // Not needed since we recalculate
    } else {
      a = c;
      // fa = fc; // Not needed since we recalculate
    }
  }

  // Return best approximation even if not fully converged
  const c = (a + b) / 2;
  if (Math.abs(f(c)) < tolerance * 10) {
    return c;
  }

  return null;
}

/**
 * Solve equation f(x) = 0 using Newton-Raphson method.
 * Fast convergence but requires good initial guess and derivative.
 *
 * @returns Root if found, null otherwise
 */
function solveNewtonRaphson(
  f: (x: number) => number,
  fPrime: (x: number) => number,
  x0: number,
  tolerance: number,
  maxIterations: number,
): number | null {
  let x = x0;

  for (let iter = 0; iter < maxIterations; iter++) {
    const fx = f(x);
    const fpx = fPrime(x);

    // Check for zero derivative
    if (Math.abs(fpx) < 1e-15) {
      return null;
    }

    const dx = fx / fpx;
    x = x - dx;

    // Check convergence
    if (Math.abs(dx) < tolerance || Math.abs(fx) < tolerance) {
      return x;
    }

    // Check for divergence
    if (!isFinite(x) || Math.abs(x) > 1e10) {
      return null;
    }
  }

  // Check if we're close enough even without full convergence
  if (Math.abs(f(x)) < tolerance * 10) {
    return x;
  }

  return null;
}

/**
 * Solve equation f(x) = 0 using secant method.
 * Like Newton-Raphson but doesn't require derivative.
 *
 * @returns Root if found, null otherwise
 */
function solveSecant(
  f: (x: number) => number,
  x0: number,
  x1: number,
  tolerance: number,
  maxIterations: number,
): number | null {
  let xPrev = x0;
  let x = x1;

  for (let iter = 0; iter < maxIterations; iter++) {
    const fx = f(x);
    const fxPrev = f(xPrev);

    // Check for same function values (would cause division by zero)
    if (Math.abs(fx - fxPrev) < 1e-15) {
      return null;
    }

    const xNew = x - (fx * (x - xPrev)) / (fx - fxPrev);

    // Check convergence
    if (Math.abs(xNew - x) < tolerance || Math.abs(fx) < tolerance) {
      return xNew;
    }

    // Check for divergence
    if (!isFinite(xNew) || Math.abs(xNew) > 1e10) {
      return null;
    }

    xPrev = x;
    x = xNew;
  }

  // Check if we're close enough even without full convergence
  if (Math.abs(f(x)) < tolerance * 10) {
    return x;
  }

  return null;
}

/**
 * Lima's approximation formula for finite square well eigenvalues.
 * Reference: Lima, F. M. S. (2020). "A simpler graphical solution and an
 * approximate formula for energy eigenvalues in finite square quantum wells."
 * American Journal of Physics, 88(11), 1019.
 *
 * Provides excellent initial guesses for the numerical solvers.
 *
 * @param xi0 - Dimensionless well parameter
 * @param stateIndex - Index of the state (0, 1, 2, ...)
 * @param parity - "even" or "odd"
 * @returns Approximate ξ value, or null if not applicable
 */
function getLimaApproximation(
  xi0: number,
  stateIndex: number,
  parity: "even" | "odd",
): number | null {
  // Determine which interval to search based on parity
  const n = parity === "even" ? stateIndex * 2 : stateIndex * 2 + 1;

  // Check if this state can exist
  if ((n * Math.PI) / 2 >= xi0) {
    return null;
  }

  // Lima's approximation uses an improved formula
  // v_n ≈ (n+1/2)π - arcsin((n+1/2)π/(2*xi0)) when (n+1/2)π < xi0

  const nEffective = n / 2 + 0.5; // Effective quantum number
  const vApprox = nEffective * Math.PI;

  if (vApprox >= xi0) {
    return null;
  }

  // Refined approximation based on Lima (2020)
  // This is a simplified version - the full formula is more complex
  const ratio = vApprox / xi0;

  if (ratio >= 1) {
    return null;
  }

  // Correction term
  const correction = Math.asin(ratio);
  const vRefined = vApprox - correction * 0.5; // Simplified correction

  // Ensure result is in valid range
  const xiMin = (n * Math.PI) / 2;
  const xiMax = ((n + 1) * Math.PI) / 2;

  if (vRefined >= xiMin && vRefined <= Math.min(xiMax, xi0)) {
    return vRefined;
  }

  // Fallback to simple linear interpolation
  const alpha = Math.min(0.5 + 0.3 * (1 - ratio), 0.95);
  const vSimple = xiMin + alpha * (Math.min(xiMax, xi0) - xiMin);

  return vSimple;
}
