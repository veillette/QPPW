/**
 * Analytical solutions for well-known quantum potentials.
 * These provide exact solutions without numerical approximation.
 */

import QuantumConstants from "./QuantumConstants.js";
import { BoundStateResult, GridConfig } from "./PotentialFunction.js";
import qppw from "../../QPPWNamespace.js";

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

/**
 * Analytical solution for the 1D Coulomb potential.
 * V(x) = -α/|x|
 *
 * This potential has a singularity at x=0 and describes a 1D hydrogen-like atom.
 * The energy eigenvalues are given by E_n = -mα²/(2ℏ²(n+1/2)²)
 *
 * @param coulombStrength - Coulomb strength parameter α in J·m
 * @param mass - Particle mass in kg
 * @param numStates - Number of energy levels to calculate
 * @param gridConfig - Grid configuration for wavefunction evaluation
 * @returns Bound state results with exact energies and wavefunctions
 */
export function solveCoulomb1DPotential(
  coulombStrength: number,
  mass: number,
  numStates: number,
  gridConfig: GridConfig,
): BoundStateResult {
  const { HBAR } = QuantumConstants;
  const alpha = coulombStrength;

  // Calculate energies: E_n = -mα²/(2ℏ²(n+1/2)²) for n = 0, 1, 2, ...
  const energies: number[] = [];
  for (let n = 0; n < numStates; n++) {
    const energy = -(mass * alpha * alpha) / (2 * HBAR * HBAR * (n + 0.5) * (n + 0.5));
    energies.push(energy);
  }

  // Generate grid
  const numPoints = gridConfig.numPoints;
  const xGrid: number[] = [];
  const dx = (gridConfig.xMax - gridConfig.xMin) / (numPoints - 1);
  for (let i = 0; i < numPoints; i++) {
    xGrid.push(gridConfig.xMin + i * dx);
  }

  // Calculate effective Bohr radius for 1D: a_0 = ℏ²/(mα)
  const a0 = (HBAR * HBAR) / (mass * alpha);

  // Calculate wavefunctions
  // For 1D Coulomb, we use a hydrogen-like form with modified quantum numbers
  // ψ_n(x) ∝ exp(-|x|/n*a_0) * L_n(2|x|/(n*a_0))
  // where the effective n is (n + 1/2) for the 1D case
  const wavefunctions: number[][] = [];

  for (let n = 0; n < numStates; n++) {
    const wavefunction: number[] = [];

    // Effective principal quantum number for 1D
    const nEff = n + 0.5;
    const a_n = nEff * a0;

    // Normalization constant (simplified)
    const normalization = Math.sqrt(1.0 / (a_n * gamma(2 * n + 2)));

    for (const x of xGrid) {
      const absX = Math.abs(x);
      const rho = 2 * absX / a_n;

      // Wavefunction: ψ(x) = N * exp(-ρ/2) * L_n^1(ρ)
      const laguerre = associatedLaguerre(n, 1, rho);
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
 * Analytical solution for the 3D Coulomb potential (hydrogen atom, radial equation with L=0).
 * V(r) = -α/r
 *
 * This solves the radial Schrödinger equation for the hydrogen atom with L=0 (s-waves).
 * The energy eigenvalues are given by E_n = -mα²/(2ℏ²n²) for n = 1, 2, 3, ...
 *
 * @param coulombStrength - Coulomb strength parameter α in J·m
 * @param mass - Particle mass in kg
 * @param numStates - Number of energy levels to calculate
 * @param gridConfig - Grid configuration for radial wavefunction evaluation (r > 0)
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

  // Generate grid (must be r > 0 for radial equation)
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
    // N_n0 = 2/(n²a_0^3) * sqrt(1/n)
    const normalization = 2.0 / (a_n * Math.sqrt(a_n)) * Math.sqrt(factorial(n - 1) / (2 * n * factorial(n)));

    for (const r of xGrid) {
      if (r <= 0) {
        // Radial wavefunction must be zero at r=0 for L>0, but for L=0 it's finite
        // However, we'll handle r <= 0 by setting to 0
        wavefunction.push(0);
        continue;
      }

      const rho = 2 * r / a_n;

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

qppw.register("AnalyticalSolutions", {
  solveInfiniteWell,
  solveHarmonicOscillator,
  solveMorsePotential,
  solvePoschlTellerPotential,
  solveRosenMorsePotential,
  solveEckartPotential,
  solveCoulomb1DPotential,
  solveCoulomb3DPotential
});
