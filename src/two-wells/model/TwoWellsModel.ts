/**
 * TwoWellsModel represents the physics model for a double quantum potential well.
 * It handles quantum tunneling and double-well dynamics.
 */

import { NumberProperty, Property } from "scenerystack/axon";

export class TwoWellsModel {
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

  // Simulation state
  public readonly isPlayingProperty: Property<boolean>;
  public readonly timeProperty: NumberProperty;

  public constructor() {
    // Initialize well parameters with default values
    this.wellWidthProperty = new NumberProperty(10); // in nanometers
    this.wellDepthProperty = new NumberProperty(5); // in eV
    this.wellSeparationProperty = new NumberProperty(5); // in nanometers

    // Initialize barrier parameters
    this.barrierHeightProperty = new NumberProperty(3); // in eV
    this.barrierWidthProperty = new NumberProperty(2); // in nanometers

    // Initialize energy level
    this.energyLevelProperty = new NumberProperty(0);

    // Initialize tunneling probability
    this.tunnelingProbabilityProperty = new NumberProperty(0);

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
    this.wellSeparationProperty.reset();
    this.barrierHeightProperty.reset();
    this.barrierWidthProperty.reset();
    this.energyLevelProperty.reset();
    this.tunnelingProbabilityProperty.reset();
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
      const m = 9.10938356e-31; // Electron mass (kg)
      const eV = 1.602176634e-19; // Electron volt in joules

      const V0 = barrierHeight * eV;
      const E = energy * eV;
      const kappa = Math.sqrt((2 * m * (V0 - E)) / (hbar * hbar));
      const transmissionCoeff = Math.exp(-2 * kappa * barrierWidth);

      this.tunnelingProbabilityProperty.value = transmissionCoeff;
    }
  }
}
