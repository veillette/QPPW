/**
 * Analytical solution for the Rosen-Morse potential.
 * V(x) = -V_0 / cosh²(ax) + V_1 * tanh(ax)
 *
 * This potential is useful for modeling molecular interactions and has exact solutions.
 */

import QuantumConstants from "../QuantumConstants.js";
import { BoundStateResult, GridConfig } from "../PotentialFunction.js";
import { jacobiPolynomial, factorial, logGamma } from "./math-utilities.js";

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
