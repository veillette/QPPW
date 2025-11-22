/**
 * Analytical solution for the Rosen-Morse potential.
 * V(x) = -V_0 / cosh²(x/a) + V_1 * tanh(x/a)
 *
 * This potential is useful for modeling molecular interactions and has exact solutions.
 */

import QuantumConstants from "../QuantumConstants.js";
import { BoundStateResult, GridConfig } from "../PotentialFunction.js";
import { jacobiPolynomial, factorial, logGamma } from "./math-utilities.js";

/**
 * Analytical solution for the Rosen-Morse potential.
 * V(x) = -V_0 / cosh²(x/a) + V_1 * tanh(x/a)
 *
 * This potential is useful for modeling molecular interactions and has exact solutions.
 *
 * @param potentialDepth - Potential depth V_0 in Joules (positive value)
 * @param barrierHeight - Barrier height V_1 in Joules
 * @param wellWidth - Width parameter a in meters
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

  // Calculate dimensionless parameters (with x/a substitution)
  // λ = a * sqrt(2*m*V_0) / ℏ
  // μ = V_1 / (2 * sqrt(V_0 * ℏ²/(2*m*a²)))
  // The correct formula for μ should make it dimensionless and scale as V1/sqrt(V0)
  const lambda = (a * Math.sqrt(2 * mass * V0)) / HBAR;

  // Corrected formula: mu should be dimensionless and proportional to V1/sqrt(V0)
  // mu = (a * sqrt(2m) / (2*hbar)) * V1 / sqrt(V0)
  const mu = (a * Math.sqrt(2 * mass) * V1) / (2 * HBAR * Math.sqrt(V0));

  console.log(`Rosen-Morse parameters: lambda=${lambda.toFixed(3)}, mu=${mu.toFixed(3)}, |mu|=${Math.abs(mu).toFixed(3)}`);
  console.log(`V0=${(V0 * 6.242e18).toFixed(3)} eV, V1=${(V1 * 6.242e18).toFixed(3)} eV`);

  // For bound states, we need λ > |μ|
  if (lambda <= Math.abs(mu)) {
    throw new Error(`Rosen-Morse potential: λ=${lambda.toFixed(3)} ≤ |μ|=${Math.abs(mu).toFixed(3)}. Increase well depth or decrease barrier height.`);
  }

  // Calculate the effective parameter for bound states
  const lambdaEff = Math.sqrt(lambda * lambda - mu * mu);

  // Maximum number of bound states
  const nMax = Math.floor(lambdaEff - 0.5);
  const actualNumStates = Math.min(numStates, nMax + 1);

  if (actualNumStates <= 0) {
    throw new Error("Rosen-Morse potential: No bound states available. Increase well depth or decrease barrier height.");
  }

  // Calculate energies
  // Standard Rosen-Morse energy formula from original 1932 Rosen & Morse paper (Phys. Rev. 42, 210)
  // For hyperbolic potential V(x) = -V₀/cosh²(x/a) + V₁·tanh(x/a)
  // Energy eigenvalues: E_n = -(ℏ²/2ma²) * (λ_eff - n - 1/2)²
  // where λ_eff = √(λ² - μ²), λ = a√(2mV₀)/ℏ, μ = aV₁√(2m)/(2ℏ√V₀)
  // Ground state (n=0) has lowest energy, all states have E < 0 for the symmetric case
  const energies: number[] = [];
  const energyFactor = (HBAR * HBAR) / (2 * mass * a * a);

  for (let n = 0; n < actualNumStates; n++) {
    const s_n = lambdaEff - n - 0.5;
    // Standard formula: E_n = -(ℏ²/2ma²) * s_n²
    const energy = -energyFactor * s_n * s_n;

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
  // ψ_n(x) = N_n * sech^(λ_eff-n-1/2)(x/a) * exp(μ*tanh(x/a)) * P_n^(α,β)(tanh(x/a))
  // where α = λ_eff - n - 1/2 - μ, β = λ_eff - n - 1/2 + μ
  const wavefunctions: number[][] = [];

  for (let n = 0; n < actualNumStates; n++) {
    const wavefunction: number[] = [];
    const s = lambdaEff - n - 0.5;
    const alpha_jac = s - mu;
    const beta_jac = s + mu;

    // Normalization (with 1/a factor from variable change)
    const normalization = Math.sqrt((1 / a) * (2 * s) / (factorial(n) * Math.exp(logGamma(n + alpha_jac + 1) + logGamma(n + beta_jac + 1) - logGamma(n + alpha_jac + beta_jac + 1))));

    for (const x of xGrid) {
      const tanhVal = Math.tanh(x / a);
      const sechVal = 1.0 / Math.cosh(x / a);

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
