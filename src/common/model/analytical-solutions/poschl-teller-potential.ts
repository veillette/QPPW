/**
 * Analytical solution for the Pöschl-Teller potential.
 * V(x) = -V_0 / cosh²(ax)
 *
 * This potential is useful for modeling quantum wells and has exact solutions.
 */

import QuantumConstants from "../QuantumConstants.js";
import { BoundStateResult, GridConfig } from "../PotentialFunction.js";
import { jacobiPolynomial, factorial } from "./math-utilities.js";

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
