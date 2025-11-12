/**
 * OneWellModel represents the physics model for a single quantum potential well.
 * It handles the quantum mechanical calculations for a particle in a single well.
 */

import { NumberProperty, Property } from "scenerystack/axon";

export class OneWellModel {
  // Well parameters
  public readonly wellWidthProperty: NumberProperty;
  public readonly wellDepthProperty: NumberProperty;

  // Energy level
  public readonly energyLevelProperty: NumberProperty;

  // Simulation state
  public readonly isPlayingProperty: Property<boolean>;
  public readonly timeProperty: NumberProperty;

  public constructor() {
    // Initialize well parameters with default values
    this.wellWidthProperty = new NumberProperty(10); // in nanometers
    this.wellDepthProperty = new NumberProperty(5); // in eV

    // Initialize energy level (ground state by default)
    this.energyLevelProperty = new NumberProperty(0);

    // Initialize simulation state
    this.isPlayingProperty = new Property<boolean>(false);
    this.timeProperty = new NumberProperty(0);
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
   * For an infinite square well: E_n = (n^2 * h^2) / (8 * m * L^2)
   * @param n - The quantum number (1, 2, 3, ...)
   * @returns The energy in eV
   */
  public getEnergyLevel(n: number): number {
    const hbar = 1.054571817e-34; // Reduced Planck constant (J·s)
    const m = 9.10938356e-31; // Electron mass (kg)
    const L = this.wellWidthProperty.value * 1e-9; // Convert nm to m
    const eV = 1.602176634e-19; // Electron volt in joules

    const energy = (n * n * hbar * hbar * Math.PI * Math.PI) / (2 * m * L * L);
    return energy / eV; // Convert to eV
  }
}
