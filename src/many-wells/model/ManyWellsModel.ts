/**
 * ManyWellsModel represents the physics model for multiple quantum potential wells.
 * It demonstrates energy bands and band gaps in a periodic potential (solid-state physics).
 */

import { NumberProperty } from "scenerystack/axon";
import { Range } from "scenerystack/dot";
import { BaseModel } from "../../common/model/BaseModel.js";
import { NumericalMethod } from "../../common/model/Schrodinger1DSolver.js";

export class ManyWellsModel extends BaseModel {
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

  public constructor() {
    super();
    // Initialize well parameters with default values
    this.wellWidthProperty = new NumberProperty(1.0, { range: new Range(0.1, 6.0) }); // in nanometers (max 6 nm)
    this.wellDepthProperty = new NumberProperty(5.0, { range: new Range(0.1, 15.0) }); // in eV
    this.numberOfWellsProperty = new NumberProperty(5, { range: new Range(2, 20) }); // number of wells

    // Initialize barrier parameters
    this.barrierWidthProperty = new NumberProperty(2); // in nanometers
    this.barrierHeightProperty = new NumberProperty(3); // in eV

    // Initialize lattice constant (well width + barrier width)
    this.latticeConstantProperty = new NumberProperty(12); // in nanometers

    // Initialize selected energy band
    this.selectedBandProperty = new NumberProperty(0);
  }

  /**
   * Called when the solver method changes.
   * @param _method - The new numerical method (unused)
   */
   
  protected onSolverMethodChanged(_method: NumericalMethod): void {
    // No cached results to invalidate for many wells model yet
  }

  /**
   * Resets the model to its initial state.
   * This is the public API method that delegates to resetAll().
   */
  public reset(): void {
    this.resetAll();
  }

  /**
   * Resets all properties to their initial state.
   * Override from BaseModel.
   */
  public override resetAll(): void {
    super.resetAll();
    this.wellWidthProperty.reset();
    this.wellDepthProperty.reset();
    this.numberOfWellsProperty.reset();
    this.barrierWidthProperty.reset();
    this.barrierHeightProperty.reset();
    this.latticeConstantProperty.reset();
    this.selectedBandProperty.reset();
  }

  /**
   * Steps the model forward in time.
   * @param dt - The time step in seconds (can be negative for backward stepping)
   * @param forced - If true, steps even when paused (for manual stepping buttons)
   */
  public override step(dt: number, forced = false): void {
    super.step(dt, forced);
    // Add Bloch wave dynamics here
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
    const electronMass = 9.10938356e-31; // Electron mass (kg)
    const latticeConstantMeters = this.latticeConstantProperty.value * 1e-9; // Lattice constant in meters
    const eV = 1.602176634e-19; // Electron volt in joules

    // Simplified dispersion relation for a periodic potential
    // E(k) = E_0 + (hbar^2 * k^2) / (2m) + V_0 * cos(k * a)
    const E0 = (bandIndex + 1) * 2; // Base energy for this band
    const kineticEnergy =
      (hbar * hbar * k * k) / (2 * electronMass * eV);
    const periodicTerm =
      (this.wellDepthProperty.value / 4) * Math.cos(k * latticeConstantMeters);

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
