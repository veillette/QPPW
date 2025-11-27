/**
 * Helper functions for computing Fourier transforms of wavefunctions.
 *
 * This module provides default implementations for analytical solutions
 * that don't have closed-form Fourier transforms.
 */

import {
  BoundStateResult,
  FourierTransformResult,
} from "../PotentialFunction.js";
import { numericalFourierTransform } from "../LinearAlgebraUtils.js";

/**
 * Compute numerical Fourier transform for a bound state result.
 *
 * This is a generic implementation that uses FFT to transform position-space
 * wavefunctions to momentum space. It can be used by any analytical solution
 * that doesn't have a simple closed-form Fourier transform.
 *
 * @param boundStateResult - Position-space wavefunction results
 * @param mass - Particle mass in kg
 * @param energyScale - Typical energy scale for the potential (used to estimate pMax)
 * @param numMomentumPoints - Optional number of momentum points
 * @param pMax - Optional maximum momentum
 * @returns Fourier transform result
 */
export function computeNumericalFourierTransform(
  boundStateResult: BoundStateResult,
  mass: number,
  energyScale: number,
  numMomentumPoints?: number,
  pMax?: number,
): FourierTransformResult {
  const numStates = boundStateResult.energies.length;

  // Determine number of momentum points
  const nMomentum = numMomentumPoints || boundStateResult.xGrid.length;

  // Determine pMax if not provided
  // Use energy scale to estimate typical momentum: p ~ sqrt(2 * m * E)
  const defaultPMax = Math.sqrt(2 * mass * Math.abs(energyScale)) * 5;
  const actualPMax = pMax || defaultPMax;

  // Use numerical Fourier transform for each state
  const pGrid: number[] = [];
  const momentumWavefunctions: number[][] = [];

  for (let stateIndex = 0; stateIndex < numStates; stateIndex++) {
    const psi = boundStateResult.wavefunctions[stateIndex];
    const xGrid = boundStateResult.xGrid;

    const { pGrid: pGridState, phiP } = numericalFourierTransform(
      psi,
      xGrid,
      mass,
      nMomentum,
      actualPMax,
    );

    if (stateIndex === 0) {
      // Use pGrid from first state
      pGrid.push(...pGridState);
    }

    momentumWavefunctions.push(phiP);
  }

  return {
    pGrid,
    momentumWavefunctions,
    method: "numerical",
  };
}
