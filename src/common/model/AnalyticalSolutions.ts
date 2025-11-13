/**
 * Analytical solutions for well-known quantum potentials.
 * These provide exact solutions without numerical approximation.
 */

import QuantumConstants from "./QuantumConstants.js";
import { BoundStateResult, GridConfig } from "./PotentialFunction.js";
import qppw from "../../QPPWNamespace.js";

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

/**
 * Analytical solution for an infinite square well (particle in a box).
 * V(x) = 0 for 0 < x < L, V(x) = ∞ otherwise
 *
 * @param wellWidth - Width of the well (L) in meters
 * @param mass - Particle mass in kg
 * @param numStates - Number of energy levels to calculate
 * @param gridConfig - Grid configuration for wavefunction evaluation
 * @returns Bound state results with exact energies and wavefunctions
 */
export function solveInfiniteWell(
  wellWidth: number,
  mass: number,
  numStates: number,
  gridConfig: GridConfig,
): BoundStateResult {
  const { HBAR } = QuantumConstants;
  const L = wellWidth;

  // Calculate energies: E_n = (n^2 * π^2 * ℏ^2) / (2 * m * L^2) for n = 1, 2, 3, ...
  const energies: number[] = [];
  for (let n = 1; n <= numStates; n++) {
    const energy =
      (n * n * Math.PI * Math.PI * HBAR * HBAR) / (2 * mass * L * L);
    energies.push(energy);
  }

  // Generate grid
  const numPoints = gridConfig.numPoints;
  const xGrid: number[] = [];
  const dx = (gridConfig.xMax - gridConfig.xMin) / (numPoints - 1);
  for (let i = 0; i < numPoints; i++) {
    xGrid.push(gridConfig.xMin + i * dx);
  }

  // Calculate wavefunctions: ψ_n(x) = sqrt(2/L) * sin(n * π * x / L)
  const wavefunctions: number[][] = [];
  for (let n = 1; n <= numStates; n++) {
    const wavefunction: number[] = [];
    const normalization = Math.sqrt(2 / L);

    for (const x of xGrid) {
      // Well is centered, so shift x coordinate
      const xShifted = x - gridConfig.xMin;
      if (xShifted >= 0 && xShifted <= L) {
        const value = normalization * Math.sin((n * Math.PI * xShifted) / L);
        wavefunction.push(value);
      } else {
        wavefunction.push(0);
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

/**
 * Analytical solution for a quantum harmonic oscillator.
 * V(x) = (1/2) * k * x^2 = (1/2) * m * ω^2 * x^2
 *
 * @param springConstant - Spring constant k in N/m
 * @param mass - Particle mass in kg
 * @param numStates - Number of energy levels to calculate
 * @param gridConfig - Grid configuration for wavefunction evaluation
 * @returns Bound state results with exact energies and wavefunctions
 */
export function solveHarmonicOscillator(
  springConstant: number,
  mass: number,
  numStates: number,
  gridConfig: GridConfig,
): BoundStateResult {
  const { HBAR } = QuantumConstants;
  const omega = Math.sqrt(springConstant / mass);

  // Calculate energies: E_n = ℏω(n + 1/2) for n = 0, 1, 2, ...
  const energies: number[] = [];
  for (let n = 0; n < numStates; n++) {
    const energy = HBAR * omega * (n + 0.5);
    energies.push(energy);
  }

  // Generate grid
  const numPoints = gridConfig.numPoints;
  const xGrid: number[] = [];
  const dx = (gridConfig.xMax - gridConfig.xMin) / (numPoints - 1);
  for (let i = 0; i < numPoints; i++) {
    xGrid.push(gridConfig.xMin + i * dx);
  }

  // Calculate wavefunctions using Hermite polynomials
  // ψ_n(x) = (1/√(2^n n!)) * (mω/πℏ)^(1/4) * exp(-mωx^2/(2ℏ)) * H_n(√(mω/ℏ) x)
  const wavefunctions: number[][] = [];
  const alpha = Math.sqrt((mass * omega) / HBAR);

  for (let n = 0; n < numStates; n++) {
    const wavefunction: number[] = [];
    const normalization =
      (1 / Math.sqrt(Math.pow(2, n) * factorial(n))) *
      Math.pow(alpha / Math.PI, 0.25);

    for (const x of xGrid) {
      const xi = alpha * x;
      const hermite = hermitePolynomial(n, xi);
      const value =
        normalization * Math.exp(-xi * xi / 2) * hermite;
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
 * Calculate the Hermite polynomial H_n(x) using recurrence relation.
 * H_0(x) = 1
 * H_1(x) = 2x
 * H_(n+1)(x) = 2x*H_n(x) - 2n*H_(n-1)(x)
 */
function hermitePolynomial(n: number, x: number): number {
  if (n === 0) return 1;
  if (n === 1) return 2 * x;

  let H_prev = 1;
  let H_curr = 2 * x;

  for (let i = 1; i < n; i++) {
    const H_next = 2 * x * H_curr - 2 * i * H_prev;
    H_prev = H_curr;
    H_curr = H_next;
  }

  return H_curr;
}

/**
 * Calculate factorial n!
 */
function factorial(n: number): number {
  if (n <= 1) return 1;
  let result = 1;
  for (let i = 2; i <= n; i++) {
    result *= i;
  }
  return result;
}

/**
 * Analytical solution for the Morse potential.
 * V(x) = D_e * (1 - exp(-a(x - x_e)))^2
 *
 * The Morse potential describes molecular vibrations more accurately than the harmonic oscillator
 * by including anharmonic effects and bond dissociation.
 *
 * @param dissociationEnergy - Dissociation energy D_e in Joules
 * @param wellWidth - Width parameter a (inverse meters)
 * @param equilibriumPosition - Equilibrium position x_e in meters
 * @param mass - Particle mass in kg
 * @param numStates - Number of energy levels to calculate
 * @param gridConfig - Grid configuration for wavefunction evaluation
 * @returns Bound state results with exact energies and wavefunctions
 */
export function solveMorsePotential(
  dissociationEnergy: number,
  wellWidth: number,
  equilibriumPosition: number,
  mass: number,
  numStates: number,
  gridConfig: GridConfig,
): BoundStateResult {
  const { HBAR } = QuantumConstants;
  const De = dissociationEnergy;
  const a = wellWidth;
  const xe = equilibriumPosition;

  // Calculate the maximum quantum number
  // n_max = floor(sqrt(2*m*D_e)/(a*ℏ) - 1/2)
  const nMax = Math.floor(Math.sqrt(2 * mass * De) / (a * HBAR) - 0.5);
  const actualNumStates = Math.min(numStates, nMax + 1);

  if (actualNumStates <= 0) {
    throw new Error("Morse potential too shallow to support bound states");
  }

  // Calculate the characteristic frequency
  const omega = a * Math.sqrt(2 * De / mass);

  // Calculate energies: E_n = ℏω(n + 1/2) - (ℏω)²(n + 1/2)² / (4*D_e)
  // Relative to the bottom of the well
  const energies: number[] = [];
  for (let n = 0; n < actualNumStates; n++) {
    const term1 = HBAR * omega * (n + 0.5);
    const term2 = (HBAR * HBAR * omega * omega * (n + 0.5) * (n + 0.5)) / (4 * De);
    const energy = term1 - term2 - De; // Energy relative to dissociation limit
    energies.push(energy);
  }

  // Generate grid
  const numPoints = gridConfig.numPoints;
  const xGrid: number[] = [];
  const dx = (gridConfig.xMax - gridConfig.xMin) / (numPoints - 1);
  for (let i = 0; i < numPoints; i++) {
    xGrid.push(gridConfig.xMin + i * dx);
  }

  // Calculate wavefunctions using associated Laguerre polynomials
  // ψ_n(z) = N_n * z^((λ-n-1/2)) * exp(-z/2) * L_n^(2λ-2n-1)(z)
  // where z = 2λ * exp(-a(x-xe)), λ = sqrt(2*m*D_e)/(a*ℏ)
  const wavefunctions: number[][] = [];
  const lambda = Math.sqrt(2 * mass * De) / (a * HBAR);

  for (let n = 0; n < actualNumStates; n++) {
    const wavefunction: number[] = [];

    // Normalization constant
    const alpha = 2 * lambda - 2 * n - 1;
    const normalization = Math.sqrt(
      (factorial(n) * a) / (gamma(2 * lambda - n))
    );

    for (const x of xGrid) {
      const z = 2 * lambda * Math.exp(-a * (x - xe));

      // Calculate wavefunction
      const exponent = lambda - n - 0.5;
      const laguerre = associatedLaguerre(n, alpha, z);
      const value = normalization * Math.pow(z, exponent) * Math.exp(-z / 2) * laguerre;

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
 * Analytical solution for the Pöschl-Teller potential.
 * V(x) = -V_0 / cosh²(ax)
 *
 * This potential is useful for modeling quantum wells and has exact solutions.
 *
 * @param potentialDepth - Potential depth V_0 in Joules (positive value)
 * @param wellWidth - Width parameter a (inverse meters)
 * @param mass - Particle mass in kg
 * @param numStates - Number of energy levels to calculate
 * @param gridConfig - Grid configuration for wavefunction evaluation
 * @returns Bound state results with exact energies and wavefunctions
 */
export function solvePoschlTellerPotential(
  potentialDepth: number,
  wellWidth: number,
  mass: number,
  numStates: number,
  gridConfig: GridConfig,
): BoundStateResult {
  const { HBAR } = QuantumConstants;
  const V0 = potentialDepth;
  const a = wellWidth;

  // Calculate λ = sqrt(2*m*V_0) / (a*ℏ)
  const lambda = Math.sqrt(2 * mass * V0) / (a * HBAR);

  // Maximum number of bound states
  const nMax = Math.floor(lambda - 0.5);
  const actualNumStates = Math.min(numStates, nMax + 1);

  if (actualNumStates <= 0) {
    throw new Error("Pöschl-Teller potential too shallow to support bound states");
  }

  // Calculate energies: E_n = -V_0 + (V_0/λ²)(λ - n - 1/2)²
  const energies: number[] = [];
  for (let n = 0; n < actualNumStates; n++) {
    const term = lambda - n - 0.5;
    const energy = -V0 + (V0 / (lambda * lambda)) * term * term;
    energies.push(energy);
  }

  // Generate grid
  const numPoints = gridConfig.numPoints;
  const xGrid: number[] = [];
  const dx = (gridConfig.xMax - gridConfig.xMin) / (numPoints - 1);
  for (let i = 0; i < numPoints; i++) {
    xGrid.push(gridConfig.xMin + i * dx);
  }

  // Calculate wavefunctions
  // ψ_n(x) = N_n * sech^(λ-n-1/2)(ax) * P_n^(λ-n-1/2, λ-n-1/2)(tanh(ax))
  // where P is the Jacobi polynomial
  const wavefunctions: number[][] = [];

  for (let n = 0; n < actualNumStates; n++) {
    const wavefunction: number[] = [];
    const alpha = lambda - n - 0.5;

    // Normalization (simplified)
    const normalization = Math.sqrt(a * (2 * alpha) / factorial(n)) * Math.sqrt(factorial(n));

    for (const x of xGrid) {
      const tanhVal = Math.tanh(a * x);
      const sechVal = 1.0 / Math.cosh(a * x);

      // Use Legendre polynomials for Jacobi with α=β
      const jacobiPoly = jacobiPolynomial(n, alpha, alpha, tanhVal);
      const value = normalization * Math.pow(sechVal, alpha) * jacobiPoly;

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
 * Calculate the associated Laguerre polynomial L_n^α(x) using recurrence relation.
 */
function associatedLaguerre(n: number, alpha: number, x: number): number {
  if (n === 0) return 1;
  if (n === 1) return 1 + alpha - x;

  let L_prev = 1;
  let L_curr = 1 + alpha - x;

  for (let k = 1; k < n; k++) {
    const L_next =
      ((2 * k + 1 + alpha - x) * L_curr - (k + alpha) * L_prev) / (k + 1);
    L_prev = L_curr;
    L_curr = L_next;
  }

  return L_curr;
}

/**
 * Calculate the Jacobi polynomial P_n^(α,β)(x) using recurrence relation.
 */
function jacobiPolynomial(n: number, alpha: number, beta: number, x: number): number {
  if (n === 0) return 1;
  if (n === 1) return 0.5 * (alpha - beta + (alpha + beta + 2) * x);

  let P_prev = 1;
  let P_curr = 0.5 * (alpha - beta + (alpha + beta + 2) * x);

  for (let k = 1; k < n; k++) {
    const a1 = 2 * (k + 1) * (k + alpha + beta + 1) * (2 * k + alpha + beta);
    const a2 = (2 * k + alpha + beta + 1) * (alpha * alpha - beta * beta);
    const a3 = (2 * k + alpha + beta) * (2 * k + alpha + beta + 1) * (2 * k + alpha + beta + 2);
    const a4 = 2 * (k + alpha) * (k + beta) * (2 * k + alpha + beta + 2);

    const P_next = ((a2 + a3 * x) * P_curr - a4 * P_prev) / a1;
    P_prev = P_curr;
    P_curr = P_next;
  }

  return P_curr;
}

/**
 * Gamma function approximation using Stirling's formula for large n
 * and direct calculation for small n.
 */
function gamma(n: number): number {
  // For integer or half-integer values
  if (n === Math.floor(n)) {
    // Integer
    return factorial(n - 1);
  } else if (n - 0.5 === Math.floor(n - 0.5)) {
    // Half-integer: Γ(n+1/2) = sqrt(π) * (2n)! / (4^n * n!)
    const k = n - 0.5;
    return Math.sqrt(Math.PI) * factorial(2 * k) / (Math.pow(4, k) * factorial(k));
  } else {
    // Use Stirling's approximation: Γ(n) ≈ sqrt(2π/n) * (n/e)^n
    return Math.sqrt(2 * Math.PI / n) * Math.pow(n / Math.E, n);
  }
}

/**
 * Analytical solution for the Rosen-Morse potential.
 * V(x) = -V_0 / cosh²(ax) + V_1 * tanh(ax)
 *
 * This potential is useful for modeling molecular interactions and has exact solutions.
 *
 * @param potentialDepth - Potential depth V_0 in Joules (positive value)
 * @param barrierHeight - Barrier height V_1 in Joules
 * @param wellWidth - Width parameter a (inverse meters)
 * @param mass - Particle mass in kg
 * @param numStates - Number of energy levels to calculate
 * @param gridConfig - Grid configuration for wavefunction evaluation
 * @returns Bound state results with exact energies and wavefunctions
 */
export function solveRosenMorsePotential(
  potentialDepth: number,
  barrierHeight: number,
  wellWidth: number,
  mass: number,
  numStates: number,
  gridConfig: GridConfig,
): BoundStateResult {
  const { HBAR } = QuantumConstants;
  const V0 = potentialDepth;
  const V1 = barrierHeight;
  const a = wellWidth;

  // Calculate dimensionless parameters
  // λ = sqrt(2*m*V_0) / (a*ℏ)
  // μ = V_1 / (2*a*ℏ * sqrt(2*m*V_0))
  const lambda = Math.sqrt(2 * mass * V0) / (a * HBAR);
  const mu = V1 / (2 * a * HBAR * Math.sqrt(2 * mass * V0));

  // For bound states, we need λ > |μ|
  if (lambda <= Math.abs(mu)) {
    throw new Error("Rosen-Morse potential too shallow to support bound states");
  }

  // Calculate the effective parameter for bound states
  const lambdaEff = Math.sqrt(lambda * lambda - mu * mu);

  // Maximum number of bound states
  const nMax = Math.floor(lambdaEff - 0.5);
  const actualNumStates = Math.min(numStates, nMax + 1);

  if (actualNumStates <= 0) {
    throw new Error("Rosen-Morse potential too shallow to support bound states");
  }

  // Calculate energies
  // E_n = -V_0 + (ℏ²a²/2m) * [(λ_eff - n - 1/2)² + μ²]
  const energies: number[] = [];
  const energyFactor = (HBAR * HBAR * a * a) / (2 * mass);

  for (let n = 0; n < actualNumStates; n++) {
    const term = lambdaEff - n - 0.5;
    const energy = -V0 + energyFactor * (term * term + mu * mu);
    energies.push(energy);
  }

  // Generate grid
  const numPoints = gridConfig.numPoints;
  const xGrid: number[] = [];
  const dx = (gridConfig.xMax - gridConfig.xMin) / (numPoints - 1);
  for (let i = 0; i < numPoints; i++) {
    xGrid.push(gridConfig.xMin + i * dx);
  }

  // Calculate wavefunctions
  // ψ_n(x) = N_n * sech^(λ_eff-n-1/2)(ax) * exp(μ*tanh(ax)) * P_n^(α,β)(tanh(ax))
  // where α = λ_eff - n - 1/2 - μ, β = λ_eff - n - 1/2 + μ
  const wavefunctions: number[][] = [];

  for (let n = 0; n < actualNumStates; n++) {
    const wavefunction: number[] = [];
    const s = lambdaEff - n - 0.5;
    const alpha_jac = s - mu;
    const beta_jac = s + mu;

    // Simplified normalization
    const normalization = Math.sqrt(a * (2 * s) / (factorial(n) * Math.exp(logGamma(n + alpha_jac + 1) + logGamma(n + beta_jac + 1) - logGamma(n + alpha_jac + beta_jac + 1))));

    for (const x of xGrid) {
      const tanhVal = Math.tanh(a * x);
      const sechVal = 1.0 / Math.cosh(a * x);

      // Calculate wavefunction
      const jacobiPoly = jacobiPolynomial(n, alpha_jac, beta_jac, tanhVal);
      const value = normalization * Math.pow(sechVal, s) * Math.exp(mu * tanhVal) * jacobiPoly;

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
 * Analytical solution for the Eckart potential.
 * V(x) = V_0 / (1 + exp(ax))² - V_1 / (1 + exp(ax))
 *
 * This potential is useful for modeling molecular barriers and chemical reactions.
 *
 * @param potentialDepth - Potential depth V_0 in Joules
 * @param barrierHeight - Barrier height V_1 in Joules
 * @param wellWidth - Width parameter a (inverse meters)
 * @param mass - Particle mass in kg
 * @param numStates - Number of energy levels to calculate
 * @param gridConfig - Grid configuration for wavefunction evaluation
 * @returns Bound state results with exact energies and wavefunctions
 */
export function solveEckartPotential(
  potentialDepth: number,
  barrierHeight: number,
  wellWidth: number,
  mass: number,
  numStates: number,
  gridConfig: GridConfig,
): BoundStateResult {
  const { HBAR } = QuantumConstants;
  const V0 = potentialDepth;
  const V1 = barrierHeight;
  const a = wellWidth;

  // Calculate dimensionless parameters
  // α = sqrt(2*m*V_0) / (a*ℏ)
  // β = V_1 * sqrt(2*m) / (2*a*ℏ*sqrt(V_0))
  const alpha_param = Math.sqrt(2 * mass * V0) / (a * HBAR);
  const beta_param = V1 * Math.sqrt(2 * mass) / (2 * a * HBAR * Math.sqrt(V0));

  // For bound states, we need certain conditions on α and β
  const s1 = -0.5 + Math.sqrt(0.25 + alpha_param);
  const s2 = -0.5 + Math.sqrt(0.25 + alpha_param - beta_param);

  if (s2 <= 0) {
    throw new Error("Eckart potential too shallow to support bound states");
  }

  // Maximum number of bound states
  const nMax = Math.floor(Math.min(s1, s2));
  const actualNumStates = Math.min(numStates, nMax);

  if (actualNumStates <= 0) {
    throw new Error("Eckart potential too shallow to support bound states");
  }

  // Calculate energies
  // E_n = -(ℏ²a²/2m) * (s_2 - n)²
  const energies: number[] = [];
  const energyFactor = (HBAR * HBAR * a * a) / (2 * mass);

  for (let n = 0; n < actualNumStates; n++) {
    const energy = -energyFactor * Math.pow(s2 - n, 2);
    energies.push(energy);
  }

  // Generate grid
  const numPoints = gridConfig.numPoints;
  const xGrid: number[] = [];
  const dx = (gridConfig.xMax - gridConfig.xMin) / (numPoints - 1);
  for (let i = 0; i < numPoints; i++) {
    xGrid.push(gridConfig.xMin + i * dx);
  }

  // Calculate wavefunctions
  // ψ_n(ξ) = N_n * ξ^(s_2-n) * (1+ξ)^(-s_1-s_2+n) * P_n^(α,β)(1-2ξ/(1+ξ))
  // where ξ = exp(ax)
  const wavefunctions: number[][] = [];

  for (let n = 0; n < actualNumStates; n++) {
    const wavefunction: number[] = [];

    const alpha_jac = 2 * (s2 - n) - 1;
    const beta_jac = 2 * (s1 - s2 + n);

    // Simplified normalization
    const normalization = Math.sqrt(
      a * factorial(n) *
      Math.exp(logGamma(2*s2 - n) + logGamma(alpha_jac + beta_jac + n + 1) -
               logGamma(n + alpha_jac + 1) - logGamma(n + beta_jac + 1))
    );

    for (const x of xGrid) {
      const xi = Math.exp(a * x);
      const xiPlus1 = 1 + xi;

      // Jacobi polynomial argument
      const jacobiArg = 1 - 2 * xi / xiPlus1;
      const jacobiPoly = jacobiPolynomial(n, alpha_jac, beta_jac, jacobiArg);

      // Calculate wavefunction
      const value = normalization *
                   Math.pow(xi, s2 - n) *
                   Math.pow(xiPlus1, -s1 - s2 + n) *
                   jacobiPoly;

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
 * Log-gamma function for better numerical stability.
 */
function logGamma(x: number): number {
  if (x <= 0) {
    throw new Error("logGamma undefined for non-positive arguments");
  }

  // For small integers, use log(factorial)
  if (x === Math.floor(x) && x < 20) {
    return Math.log(factorial(x - 1));
  }

  // Stirling's approximation: log(Γ(x)) ≈ (x-0.5)*log(x) - x + 0.5*log(2π)
  return (x - 0.5) * Math.log(x) - x + 0.5 * Math.log(2 * Math.PI);
}

qppw.register("AnalyticalSolutions", {
  solveFiniteSquareWell,
  solveInfiniteWell,
  solveHarmonicOscillator,
  solveMorsePotential,
  solvePoschlTellerPotential,
  solveRosenMorsePotential,
  solveEckartPotential
});
