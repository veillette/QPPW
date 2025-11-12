/**
 * ManyWellsModel represents the physics model for multiple quantum potential wells.
 * It demonstrates energy bands and band gaps in a periodic potential (solid-state physics).
 */

import { NumberProperty, Property } from "scenerystack/axon";

export class ManyWellsModel {
  // Well parameters
  public readonly wellWidthProperty: NumberProperty;
  public readonly wellDepthProperty: NumberProperty;
  public readonly numberOfWellsProperty: NumberProperty;

  // Barrier parameters
  public readonly barrierWidthProperty: NumberProperty;
  public readonly barrierHeightProperty: NumberProperty;

  // Lattice parameter
  public readonly latticeConstantProperty: NumberProperty;

  // Energy bands
  public readonly selectedBandProperty: NumberProperty;

  // Simulation state
  public readonly isPlayingProperty: Property<boolean>;
  public readonly timeProperty: NumberProperty;

  public constructor() {
    // Initialize well parameters with default values
    this.wellWidthProperty = new NumberProperty(10); // in nanometers
    this.wellDepthProperty = new NumberProperty(5); // in eV
    this.numberOfWellsProperty = new NumberProperty(5); // number of wells

    // Initialize barrier parameters
    this.barrierWidthProperty = new NumberProperty(2); // in nanometers
    this.barrierHeightProperty = new NumberProperty(3); // in eV

    // Initialize lattice constant (well width + barrier width)
    this.latticeConstantProperty = new NumberProperty(12); // in nanometers

    // Initialize selected energy band
    this.selectedBandProperty = new NumberProperty(0);

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
    this.numberOfWellsProperty.reset();
    this.barrierWidthProperty.reset();
    this.barrierHeightProperty.reset();
    this.latticeConstantProperty.reset();
    this.selectedBandProperty.reset();
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
      // Add Bloch wave dynamics here
    }
  }

  /**
   * Updates the lattice constant based on well and barrier widths.
   */
  public updateLatticeConstant(): void {
    this.latticeConstantProperty.value =
      this.wellWidthProperty.value + this.barrierWidthProperty.value;
  }

  /**
   * Calculates the energy band structure for the periodic potential.
   * This is a simplified model - full band structure requires numerical methods.
   * @param bandIndex - The band index (0 for lowest band)
   * @param k - The wave vector in the first Brillouin zone
   * @returns The energy in eV
   */
  public getEnergyBand(bandIndex: number, k: number): number {
    const hbar = 1.054571817e-34; // Reduced Planck constant (J·s)
    const m = 9.10938356e-31; // Electron mass (kg)
    const a = this.latticeConstantProperty.value * 1e-9; // Lattice constant in meters
    const eV = 1.602176634e-19; // Electron volt in joules

    // Simplified dispersion relation for a periodic potential
    // E(k) = E_0 + (hbar^2 * k^2) / (2m) + V_0 * cos(k * a)
    const E0 = (bandIndex + 1) * 2; // Base energy for this band
    const kineticEnergy =
      (hbar * hbar * k * k) / (2 * m * eV);
    const periodicTerm =
      (this.wellDepthProperty.value / 4) * Math.cos(k * a);

    return E0 + kineticEnergy + periodicTerm;
  }

  /**
   * Gets the number of energy bands to display.
   */
  public getNumberOfBands(): number {
    // Show bands up to a reasonable energy
    return Math.min(5, Math.floor(this.numberOfWellsProperty.value / 2) + 1);
  }
}
