/**
 * OneWellModel represents the physics model for a single quantum potential well.
 * It handles the quantum mechanical calculations for a particle in a single well.
 */

import { NumberProperty, Property } from "scenerystack/axon";
import { Range } from "scenerystack/dot";
import Schrodinger1DSolver, { WellParameters } from "../../common/model/Schrodinger1DSolver.js";
import { PotentialType, BoundStateResult } from "../../common/model/PotentialFunction.js";
import QuantumConstants from "../../common/model/QuantumConstants.js";
import QPPWPreferences from "../../QPPWPreferences.js";

export type DisplayMode = "probabilityDensity" | "waveFunction";
export type SimulationSpeed = "normal" | "fast";

export class OneWellModel {
  // Potential type selection
  public readonly potentialTypeProperty: Property<PotentialType>;

  // Well parameters (used for different potential types)
  public readonly wellWidthProperty: NumberProperty;
  public readonly wellDepthProperty: NumberProperty;
  public readonly wellOffsetProperty: NumberProperty; // For asymmetric wells

  // Particle properties
  public readonly particleMassProperty: NumberProperty; // In units of electron mass

  // Energy level selection
  public readonly selectedEnergyLevelIndexProperty: NumberProperty; // 0-indexed
  public readonly energyLevelProperty: NumberProperty; // Deprecated, kept for compatibility

  // Display settings
  public readonly displayModeProperty: Property<DisplayMode>;
  public readonly showRealPartProperty: Property<boolean>;
  public readonly showImaginaryPartProperty: Property<boolean>;
  public readonly showMagnitudeProperty: Property<boolean>;
  public readonly showPhaseProperty: Property<boolean>;

  // Chart visibility
  public readonly showTotalEnergyProperty: Property<boolean>;
  public readonly showPotentialEnergyProperty: Property<boolean>;

  // Simulation state
  public readonly isPlayingProperty: Property<boolean>;
  public readonly timeProperty: NumberProperty; // In femtoseconds
  public readonly simulationSpeedProperty: Property<SimulationSpeed>;

  // Solver for quantum calculations
  private readonly solver: Schrodinger1DSolver;

  // Cached bound state results
  private boundStateResult: BoundStateResult | null = null;

  public constructor() {
    // Initialize potential type (square/infinite well by default)
    this.potentialTypeProperty = new Property<PotentialType>(PotentialType.INFINITE_WELL);

    // Initialize well parameters with default values
    this.wellWidthProperty = new NumberProperty(1.0, { range: new Range(0.1, 10.0) }); // in nanometers
    this.wellDepthProperty = new NumberProperty(5.0, { range: new Range(0.1, 20.0) }); // in eV
    this.wellOffsetProperty = new NumberProperty(0.5, { range: new Range(0.0, 1.0) }); // normalized position

    // Initialize particle mass (1.0 = electron mass)
    this.particleMassProperty = new NumberProperty(1.0, { range: new Range(0.1, 10.0) });

    // Initialize energy level selection (ground state by default)
    this.selectedEnergyLevelIndexProperty = new NumberProperty(0, { range: new Range(0, 9) });
    this.energyLevelProperty = new NumberProperty(0); // Deprecated

    // Initialize display settings
    this.displayModeProperty = new Property<DisplayMode>("probabilityDensity");
    this.showRealPartProperty = new Property<boolean>(true);
    this.showImaginaryPartProperty = new Property<boolean>(false);
    this.showMagnitudeProperty = new Property<boolean>(false);
    this.showPhaseProperty = new Property<boolean>(false);

    // Initialize chart visibility
    this.showTotalEnergyProperty = new Property<boolean>(true);
    this.showPotentialEnergyProperty = new Property<boolean>(true);

    // Initialize simulation state
    this.isPlayingProperty = new Property<boolean>(false);
    this.timeProperty = new NumberProperty(0); // in femtoseconds
    this.simulationSpeedProperty = new Property<SimulationSpeed>("normal");

    // Initialize solver with user's preferred method
    this.solver = new Schrodinger1DSolver();

    // Update solver method when preference changes
    QPPWPreferences.numericalMethodProperty.link((method) => {
      this.solver.setNumericalMethod(method);
      this.boundStateResult = null; // Invalidate cache
    });

    // Recalculate bound states when parameters change
    const invalidateCache = () => {
      this.boundStateResult = null;
    };

    this.potentialTypeProperty.link(invalidateCache);
    this.wellWidthProperty.link(invalidateCache);
    this.wellDepthProperty.link(invalidateCache);
    this.wellOffsetProperty.link(invalidateCache);
    this.particleMassProperty.link(invalidateCache);
  }

  /**
   * Resets the model to its initial state.
   */
  public reset(): void {
    this.potentialTypeProperty.reset();
    this.wellWidthProperty.reset();
    this.wellDepthProperty.reset();
    this.wellOffsetProperty.reset();
    this.particleMassProperty.reset();
    this.selectedEnergyLevelIndexProperty.reset();
    this.energyLevelProperty.reset();
    this.displayModeProperty.reset();
    this.showRealPartProperty.reset();
    this.showImaginaryPartProperty.reset();
    this.showMagnitudeProperty.reset();
    this.showPhaseProperty.reset();
    this.showTotalEnergyProperty.reset();
    this.showPotentialEnergyProperty.reset();
    this.isPlayingProperty.reset();
    this.timeProperty.reset();
    this.simulationSpeedProperty.reset();
  }

  /**
   * Steps the model forward in time.
   * @param dt - The time step in seconds
   */
  public step(dt: number): void {
    if (this.isPlayingProperty.value) {
      // Convert dt to femtoseconds and apply speed multiplier
      const speedMultiplier = this.simulationSpeedProperty.value === "fast" ? 10 : 1;
      const dtFemtoseconds = (dt * 1e15) * speedMultiplier; // seconds to femtoseconds
      this.timeProperty.value += dtFemtoseconds;
      // Quantum mechanical time evolution is handled in the view layer
    }
  }

  /**
   * Restarts the simulation (resets time to zero).
   */
  public restart(): void {
    this.timeProperty.value = 0;
    this.isPlayingProperty.value = false;
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
    const wellDepth = this.wellDepthProperty.value * QuantumConstants.EV_TO_JOULES;
    const mass = this.particleMassProperty.value * QuantumConstants.ELECTRON_MASS;
    const numStates = 10; // Calculate first 10 states

    // Grid configuration with buffer on sides for finite potentials
    const bufferFactor = 2.0; // Extra space on each side
    const gridConfig = {
      xMin: -wellWidth * bufferFactor,
      xMax: wellWidth * (1 + bufferFactor),
      numPoints: 400,
    };

    try {
      // Build potential parameters based on type
      const potentialParams: WellParameters = {
        type: this.potentialTypeProperty.value,
        wellWidth: wellWidth,
      };

      // Add type-specific parameters
      switch (this.potentialTypeProperty.value) {
        case PotentialType.INFINITE_WELL:
          // No additional parameters needed
          break;
        case PotentialType.FINITE_WELL:
          potentialParams.wellDepth = wellDepth;
          break;
        case PotentialType.HARMONIC_OSCILLATOR:
          // Convert well depth to spring constant: k = mω² = m(4V₀/mL²) = 4V₀/L²
          potentialParams.springConstant = (4 * wellDepth) / (wellWidth * wellWidth);
          break;
        case PotentialType.ASYMMETRIC_TRIANGLE:
          // Slope is the field strength
          potentialParams.slope = wellDepth / wellWidth;
          break;
        case PotentialType.COULOMB_1D:
        case PotentialType.COULOMB_3D: {
          // For Coulomb potentials, use coulombStrength parameter α = k*e²
          // where k = 1/(4πε₀) ≈ 8.9875517923e9 N·m²/C²
          // α ≈ 2.307e-28 J·m for electron charge
          const k = 8.9875517923e9; // Coulomb's constant in N·m²/C²
          potentialParams.coulombStrength = k * QuantumConstants.ELEMENTARY_CHARGE * QuantumConstants.ELEMENTARY_CHARGE;
          break;
        }
        default:
          // For other potential types, use numerical solution
          break;
      }

      // Attempt analytical solution first
      this.boundStateResult = this.solver.solveAnalyticalIfPossible(
        potentialParams,
        mass,
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
