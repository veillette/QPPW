/**
 * Abstract base class for potentials with analytical solutions.
 *
 * This class extends BasePotential and wraps an AnalyticalSolution instance
 * to provide access to exact analytical solutions. It provides:
 * - Exact energy eigenvalues
 * - Exact wavefunctions
 * - Classical probability density
 * - Wavefunction zeros (nodes)
 * - Classical turning points
 * - Wavefunction derivatives
 *
 * @abstract
 */

import { BasePotential } from "./BasePotential.js";
import {
  BoundStateResult,
  GridConfig,
  PotentialFunction,
} from "../PotentialFunction.js";
import { AnalyticalSolution } from "../analytical-solutions/AnalyticalSolution.js";

export abstract class AnalyticalPotential extends BasePotential {
  protected solution: AnalyticalSolution;

  /**
   * Create a new analytical potential.
   * @param solution - The analytical solution instance for this potential
   * @param mass - Particle mass in kg
   */
  constructor(solution: AnalyticalSolution, mass: number) {
    super(mass);
    this.solution = solution;
  }

  /**
   * Solve the Schrödinger equation for this potential analytically.
   * Delegates to the wrapped AnalyticalSolution instance.
   *
   * @param numStates - Number of energy levels to calculate
   * @param gridConfig - Grid configuration for wavefunction evaluation
   * @returns Bound state results with exact energies and wavefunctions
   */
  public solve(numStates: number, gridConfig: GridConfig): BoundStateResult {
    return this.solution.solve(numStates, gridConfig);
  }

  /**
   * Create the potential function V(x).
   * Delegates to the wrapped AnalyticalSolution instance.
   */
  public createPotential(): PotentialFunction {
    return this.solution.createPotential();
  }

  /**
   * Calculate classical probability density for this potential.
   * Delegates to the wrapped AnalyticalSolution instance.
   *
   * @param energy - Energy of the particle in Joules
   * @param xGrid - Array of x positions in meters
   * @returns Array of normalized classical probability density values (in 1/meters)
   */
  public calculateClassicalProbability(
    energy: number,
    xGrid: number[],
  ): number[] {
    return this.solution.calculateClassicalProbability(energy, this.mass, xGrid);
  }

  /**
   * Calculate the positions of wavefunction zeros (nodes).
   * Delegates to the wrapped AnalyticalSolution instance.
   *
   * @param stateIndex - Index of the eigenstate (0 for ground state, 1 for first excited, etc.)
   * @param energy - Energy of the eigenstate in Joules (may be needed for some potentials)
   * @returns Array of x positions (in meters) where the wavefunction is zero
   */
  public calculateWavefunctionZeros(
    stateIndex: number,
    energy: number,
  ): number[] {
    return this.solution.calculateWavefunctionZeros(stateIndex, energy);
  }

  /**
   * Calculate the classical turning points.
   * Delegates to the wrapped AnalyticalSolution instance.
   *
   * @param energy - Energy of the particle in Joules
   * @returns Array of turning point pairs
   */
  public calculateTurningPoints(
    energy: number,
  ): Array<{ left: number; right: number }> {
    return this.solution.calculateTurningPoints(energy);
  }

  /**
   * Calculate the second derivative of the wavefunction.
   * Delegates to the wrapped AnalyticalSolution instance.
   *
   * @param stateIndex - Index of the eigenstate (0 for ground state, 1 for first excited, etc.)
   * @param xGrid - Array of x positions in meters where derivatives should be evaluated
   * @returns Array of second derivative values
   */
  public calculateWavefunctionSecondDerivative(
    stateIndex: number,
    xGrid: number[],
  ): number[] {
    return this.solution.calculateWavefunctionSecondDerivative(
      stateIndex,
      xGrid,
    );
  }

  /**
   * Always returns true for analytical potentials.
   */
  public hasAnalyticalSolution(): boolean {
    return true;
  }

  /**
   * Get the wrapped analytical solution instance.
   * This provides direct access to the solution for advanced use cases.
   */
  public getAnalyticalSolution(): AnalyticalSolution {
    return this.solution;
  }
}
