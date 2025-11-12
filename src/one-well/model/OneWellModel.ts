/**
 * OneWellModel represents the physics model for a single quantum potential well.
 * It handles the quantum mechanical calculations for a particle in a single well.
 */

import { NumberProperty, Property } from "scenerystack/axon";
import Schrodinger1DSolver from "../../common/model/Schrodinger1DSolver.js";
import { PotentialType, BoundStateResult } from "../../common/model/PotentialFunction.js";
import QuantumConstants from "../../common/model/QuantumConstants.js";
import QPPWPreferences from "../../QPPWPreferences.js";

export class OneWellModel {
  // Well parameters
  public readonly wellWidthProperty: NumberProperty;
  public readonly wellDepthProperty: NumberProperty;

  // Energy level
  public readonly energyLevelProperty: NumberProperty;

  // Simulation state
  public readonly isPlayingProperty: Property<boolean>;
  public readonly timeProperty: NumberProperty;

  // Solver for quantum calculations
  private readonly solver: Schrodinger1DSolver;

  // Cached bound state results
  private boundStateResult: BoundStateResult | null = null;

  public constructor() {
    // Initialize well parameters with default values
    this.wellWidthProperty = new NumberProperty(10); // in nanometers
    this.wellDepthProperty = new NumberProperty(5); // in eV

    // Initialize energy level (ground state by default)
    this.energyLevelProperty = new NumberProperty(0);

    // Initialize simulation state
    this.isPlayingProperty = new Property<boolean>(false);
    this.timeProperty = new NumberProperty(0);

    // Initialize solver with user's preferred method
    this.solver = new Schrodinger1DSolver();

    // Update solver method when preference changes
    QPPWPreferences.numericalMethodProperty.link((method) => {
      this.solver.setNumericalMethod(method);
      this.boundStateResult = null; // Invalidate cache
    });

    // Recalculate bound states when well parameters change
    this.wellWidthProperty.link(() => {
      this.boundStateResult = null;
    });
    this.wellDepthProperty.link(() => {
      this.boundStateResult = null;
    });
  }

  /**
   * Resets the model to its initial state.
   */
  public reset(): void {
    this.wellWidthProperty.reset();
    this.wellDepthProperty.reset();
    this.energyLevelProperty.reset();
    this.isPlayingProperty.reset();
    this.timeProperty.reset();
  }

  /**
   * Steps the model forward in time.
   * @param dt - The time step in seconds
   */
  public step(dt: number): void {
    if (this.isPlayingProperty.value) {
      this.timeProperty.value += dt;
      // Add quantum mechanical time evolution here
    }
  }

  /**
   * Calculates the energy eigenvalues for the current well parameters.
   * Uses the Schrödinger solver with analytical solution for infinite well.
   * @param n - The quantum number (1, 2, 3, ...)
   * @returns The energy in eV
   */
  public getEnergyLevel(n: number): number {
    // Ensure bound states are calculated
    if (!this.boundStateResult) {
      this.calculateBoundStates();
    }

    // Return energy for quantum number n (1-indexed)
    if (this.boundStateResult && n > 0 && n <= this.boundStateResult.energies.length) {
      const energyJoules = this.boundStateResult.energies[n - 1];
      return Schrodinger1DSolver.joulesToEV(energyJoules);
    }

    // Fallback to analytical formula if solver fails
    const L = this.wellWidthProperty.value * QuantumConstants.NM_TO_M;
    const energy =
      (n * n * Math.PI * Math.PI * QuantumConstants.HBAR * QuantumConstants.HBAR) /
      (2 * QuantumConstants.ELECTRON_MASS * L * L);
    return Schrodinger1DSolver.joulesToEV(energy);
  }

  /**
   * Get all bound state energies and wavefunctions for the current well.
   * @returns BoundStateResult with energies (in eV) and wavefunctions
   */
  public getBoundStates(): BoundStateResult | null {
    if (!this.boundStateResult) {
      this.calculateBoundStates();
    }
    return this.boundStateResult;
  }

  /**
   * Calculate bound states using the Schrödinger solver.
   * Results are cached until well parameters change.
   */
  private calculateBoundStates(): void {
    const wellWidth = this.wellWidthProperty.value * QuantumConstants.NM_TO_M;
    const numStates = 10; // Calculate first 10 states

    const gridConfig = {
      xMin: 0,
      xMax: wellWidth,
      numPoints: 200,
    };

    try {
      // Use analytical solution for infinite well
      this.boundStateResult = this.solver.solveAnalyticalIfPossible(
        {
          type: PotentialType.INFINITE_WELL,
          wellWidth: wellWidth,
        },
        QuantumConstants.ELECTRON_MASS,
        numStates,
        gridConfig,
      );
    } catch (error) {
      console.error("Error calculating bound states:", error);
      this.boundStateResult = null;
    }
  }

  /**
   * Get the wavefunction for a specific quantum number.
   * @param n - The quantum number (1, 2, 3, ...)
   * @returns Array of wavefunction values at grid points, or null if unavailable
   */
  public getWavefunction(n: number): number[] | null {
    if (!this.boundStateResult) {
      this.calculateBoundStates();
    }

    if (this.boundStateResult && n > 0 && n <= this.boundStateResult.wavefunctions.length) {
      return this.boundStateResult.wavefunctions[n - 1];
    }

    return null;
  }

  /**
   * Get the spatial grid points for wavefunctions.
   * @returns Array of x positions in nanometers
   */
  public getXGrid(): number[] | null {
    if (!this.boundStateResult) {
      this.calculateBoundStates();
    }

    if (this.boundStateResult) {
      // Convert from meters to nanometers
      return this.boundStateResult.xGrid.map((x) => x * QuantumConstants.M_TO_NM);
    }

    return null;
  }
}
