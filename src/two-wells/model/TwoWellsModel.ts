/**
 * TwoWellsModel represents the physics model for a double quantum potential well.
 * It handles quantum tunneling and double-well dynamics.
 */

import { NumberProperty } from "scenerystack/axon";
import { Range } from "scenerystack/dot";
import { BaseModel } from "../../common/model/BaseModel.js";
import { NumericalMethod } from "../../common/model/Schrodinger1DSolver.js";

export class TwoWellsModel extends BaseModel {
  // Well parameters
  public readonly wellWidthProperty: NumberProperty;
  public readonly wellDepthProperty: NumberProperty;
  public readonly wellSeparationProperty: NumberProperty;

  // Barrier parameters
  public readonly barrierHeightProperty: NumberProperty;
  public readonly barrierWidthProperty: NumberProperty;

  // Energy level
  public readonly energyLevelProperty: NumberProperty;

  // Tunneling visualization
  public readonly tunnelingProbabilityProperty: NumberProperty;

  public constructor() {
    super();
    // Initialize well parameters with default values
    this.wellWidthProperty = new NumberProperty(1.0, { range: new Range(0.1, 6.0) }); // in nanometers (max 6 nm)
    this.wellDepthProperty = new NumberProperty(5.0, { range: new Range(0.1, 15.0) }); // in eV
    this.wellSeparationProperty = new NumberProperty(5.0, { range: new Range(0.1, 10.0) }); // in nanometers

    // Initialize barrier parameters
    this.barrierHeightProperty = new NumberProperty(3); // in eV
    this.barrierWidthProperty = new NumberProperty(2); // in nanometers

    // Initialize energy level
    this.energyLevelProperty = new NumberProperty(0);

    // Initialize tunneling probability
    this.tunnelingProbabilityProperty = new NumberProperty(0);
  }

  /**
   * Called when the solver method changes.
   * @param _method - The new numerical method (unused)
   */
  // eslint-disable-next-line @typescript-eslint/no-unused-vars
  protected onSolverMethodChanged(_method: NumericalMethod): void {
    // No cached results to invalidate for two wells model yet
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
    this.wellSeparationProperty.reset();
    this.barrierHeightProperty.reset();
    this.barrierWidthProperty.reset();
    this.energyLevelProperty.reset();
    this.tunnelingProbabilityProperty.reset();
  }

  /**
   * Steps the model forward in time.
   * @param dt - The time step in seconds (can be negative for backward stepping)
   * @param forced - If true, steps even when paused (for manual stepping buttons)
   */
  public override step(dt: number, forced = false): void {
    super.step(dt, forced);
    if (this.isPlayingProperty.value || forced) {
      this.updateTunnelingProbability();
      // Add quantum tunneling dynamics here
    }
  }

  /**
   * Updates the tunneling probability based on current parameters.
   * Uses WKB approximation for quantum tunneling.
   */
  private updateTunnelingProbability(): void {
    const energy = this.energyLevelProperty.value;
    const barrierHeight = this.barrierHeightProperty.value;
    const barrierWidth = this.barrierWidthProperty.value * 1e-9; // Convert to meters

    if (energy >= barrierHeight) {
      // Classical regime - particle goes over the barrier
      this.tunnelingProbabilityProperty.value = 1.0;
    } else {
      // Quantum tunneling regime
      const hbar = 1.054571817e-34; // Reduced Planck constant (J·s)
      const electronMass = 9.10938356e-31; // Electron mass (kg)
      const eV = 1.602176634e-19; // Electron volt in joules

      const V0 = barrierHeight * eV;
      const E = energy * eV;
      const kappa = Math.sqrt((2 * electronMass * (V0 - E)) / (hbar * hbar));
      const transmissionCoeff = Math.exp(-2 * kappa * barrierWidth);

      this.tunnelingProbabilityProperty.value = transmissionCoeff;
    }
  }
}
