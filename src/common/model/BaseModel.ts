/**
 * BaseModel is an abstract base class for all physics models in QPPW.
 * It provides common functionality for time evolution, simulation control, and solver configuration.
 */

import {
  EnumerationProperty,
  NumberProperty,
  Property,
} from "scenerystack/axon";
import { Range } from "scenerystack/dot";
import { TimeSpeed } from "scenerystack";
import Schrodinger1DSolver, { NumericalMethod } from "./Schrodinger1DSolver.js";
import QPPWPreferences from "../../QPPWPreferences.js";
import { PotentialType, BoundStateResult } from "./PotentialFunction.js";
import QuantumConstants from "./QuantumConstants.js";
import { SuperpositionType, SuperpositionConfig } from "./SuperpositionType.js";

export abstract class BaseModel {
  // ==================== CONSTANTS ====================

  /**
   * Default time step for manual stepping (in seconds).
   * Corresponds to approximately 1 frame at 60 FPS.
   */
  public static readonly MANUAL_STEP_SIZE = 0.016;

  /**
   * Minimum value for well width in nanometers.
   * Ensures the well is large enough for meaningful quantum behavior.
   */
  private static readonly WELL_WIDTH_MIN = 0.1;

  /**
   * Maximum value for well width in nanometers.
   * Limits the computational domain to reasonable sizes.
   */
  private static readonly WELL_WIDTH_MAX = 6.0;

  /**
   * Minimum value for well depth in electron volts.
   * Ensures at least one bound state can exist.
   */
  private static readonly WELL_DEPTH_MIN = 0.1;

  /**
   * Maximum value for well depth in electron volts.
   * Limits energy to reasonable values for atomic-scale systems.
   */
  private static readonly WELL_DEPTH_MAX = 15.0;

  /**
   * Minimum value for well offset (normalized position).
   * Used for asymmetric potential configurations.
   */
  private static readonly WELL_OFFSET_MIN = 0.0;

  /**
   * Maximum value for well offset (normalized position).
   * Used for asymmetric potential configurations.
   */
  private static readonly WELL_OFFSET_MAX = 1.0;

  /**
   * Minimum particle mass in units of electron mass.
   * Allows lighter particles (e.g., 0.5 m_e).
   */
  private static readonly PARTICLE_MASS_MIN = 0.5;

  /**
   * Maximum particle mass in units of electron mass.
   * Limits to slightly heavier than electron mass.
   */
  private static readonly PARTICLE_MASS_MAX = 1.1;

  /**
   * Minimum energy level index (0 = ground state).
   * Energy levels are 0-indexed internally.
   */
  private static readonly ENERGY_LEVEL_INDEX_MIN = 0;

  /**
   * Maximum energy level index.
   * Caps the number of accessible quantum states.
   */
  private static readonly ENERGY_LEVEL_INDEX_MAX = 99;

  /**
   * Minimum kinetic energy threshold as fraction of maximum kinetic energy.
   * Used in classical probability calculations to prevent singularities at turning points.
   * Value of 0.01 (1%) prevents infinities while preserving probability distribution shape.
   */
  private static readonly MIN_KINETIC_ENERGY_FRACTION = 0.01;

  /**
   * Speed multiplier for slow time evolution.
   * Slows down the simulation by a factor of 10 for better observation.
   */
  private static readonly SLOW_SPEED_MULTIPLIER = 0.1;

  /**
   * Divisor for trapezoidal integration (averaging adjacent values).
   * Used in numerical integration: (f(x_i) + f(x_{i+1})) / 2.
   */
  private static readonly TRAPEZOIDAL_DIVISOR = 2;

  // ==================== PROPERTIES ====================

  // Simulation state properties
  public readonly isPlayingProperty: Property<boolean>;
  public readonly timeProperty: NumberProperty; // In femtoseconds
  public readonly timeSpeedProperty: EnumerationProperty<TimeSpeed>;

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

  // Display settings
  public readonly showRealPartProperty: Property<boolean>;
  public readonly showImaginaryPartProperty: Property<boolean>;
  public readonly showMagnitudeProperty: Property<boolean>;
  public readonly showPhaseProperty: Property<boolean>;
  public readonly showClassicalProbabilityProperty: Property<boolean>;
  public readonly showZerosProperty: Property<boolean>;

  // Chart visibility
  public readonly showTotalEnergyProperty: Property<boolean>;
  public readonly showPotentialEnergyProperty: Property<boolean>;

  // Superposition state
  public readonly superpositionTypeProperty: Property<SuperpositionType>;
  public readonly superpositionConfigProperty: Property<SuperpositionConfig>;

  // Cached bound state results
  protected boundStateResult: BoundStateResult | null = null;

  // Solver for quantum calculations
  protected readonly solver: Schrodinger1DSolver;

  // Guard flag to prevent reentry in step method
  private isStepping: boolean = false;

  protected constructor() {
    // Initialize simulation state
    this.isPlayingProperty = new Property<boolean>(false);
    this.timeProperty = new NumberProperty(0); // in femtoseconds
    this.timeSpeedProperty = new EnumerationProperty(TimeSpeed.NORMAL);

    // Initialize potential type
    this.potentialTypeProperty = new Property<PotentialType>(
      PotentialType.INFINITE_WELL,
    );

    // Initialize well parameters with default values
    this.wellWidthProperty = new NumberProperty(1.0, {
      range: new Range(BaseModel.WELL_WIDTH_MIN, BaseModel.WELL_WIDTH_MAX),
    }); // in nanometers
    this.wellDepthProperty = new NumberProperty(5.0, {
      range: new Range(BaseModel.WELL_DEPTH_MIN, BaseModel.WELL_DEPTH_MAX),
    }); // in eV
    this.wellOffsetProperty = new NumberProperty(0.5, {
      range: new Range(BaseModel.WELL_OFFSET_MIN, BaseModel.WELL_OFFSET_MAX),
    }); // normalized position

    // Initialize particle mass (1.0 = electron mass)
    this.particleMassProperty = new NumberProperty(1.0, {
      range: new Range(
        BaseModel.PARTICLE_MASS_MIN,
        BaseModel.PARTICLE_MASS_MAX,
      ),
    }); // 0.5 to 1.1 times electron mass

    // Initialize energy level selection (ground state by default)
    this.selectedEnergyLevelIndexProperty = new NumberProperty(0, {
      range: new Range(
        BaseModel.ENERGY_LEVEL_INDEX_MIN,
        BaseModel.ENERGY_LEVEL_INDEX_MAX,
      ),
    });

    // Initialize display settings
    this.showRealPartProperty = new Property<boolean>(true);
    this.showImaginaryPartProperty = new Property<boolean>(false);
    this.showMagnitudeProperty = new Property<boolean>(false);
    this.showPhaseProperty = new Property<boolean>(false);
    this.showClassicalProbabilityProperty = new Property<boolean>(false);
    this.showZerosProperty = new Property<boolean>(false);

    // Initialize chart visibility
    this.showTotalEnergyProperty = new Property<boolean>(true);
    this.showPotentialEnergyProperty = new Property<boolean>(true);

    // Initialize superposition state
    this.superpositionTypeProperty = new Property<SuperpositionType>(
      SuperpositionType.SINGLE,
    );
    this.superpositionConfigProperty = new Property<SuperpositionConfig>({
      type: SuperpositionType.SINGLE,
      amplitudes: [1.0],
      phases: [0],
    });

    // Initialize solver with user's preferred method
    this.solver = new Schrodinger1DSolver();

    // Update solver method when preference changes
    QPPWPreferences.numericalMethodProperty.link((method: NumericalMethod) => {
      this.solver.setNumericalMethod(method);
      this.onSolverMethodChanged(method);
    });

    // Invalidate caches when grid points preference changes
    // Use lazyLink to avoid triggering during initialization
    QPPWPreferences.gridPointsProperty.lazyLink(() => {
      this.onSolverMethodChanged(this.solver.getNumericalMethod());
    });

    // Note: setupCacheInvalidation() must be called by subclasses
    // after all their properties are initialized
  }

  /**
   * Setup cache invalidation listeners for common properties.
   * Subclasses can override this to add additional invalidation listeners.
   */
  protected setupCacheInvalidation(): void {
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
   * Called when the solver method changes.
   * Subclasses should override this to invalidate caches or recalculate results.
   * @param method - The new numerical method
   */
  protected abstract onSolverMethodChanged(method: NumericalMethod): void;

  /**
   * Resets the model to its initial state.
   * This is the public API method that delegates to resetAll().
   */
  public reset(): void {
    this.resetAll();
  }

  /**
   * Resets all properties to their initial state.
   * Subclasses should override this method and call super.resetAll().
   */
  public resetAll(): void {
    this.isPlayingProperty.reset();
    this.timeProperty.reset();
    this.timeSpeedProperty.reset();
    this.potentialTypeProperty.reset();
    this.wellWidthProperty.reset();
    this.wellDepthProperty.reset();
    this.wellOffsetProperty.reset();
    this.particleMassProperty.reset();
    this.selectedEnergyLevelIndexProperty.reset();
    this.showRealPartProperty.reset();
    this.showImaginaryPartProperty.reset();
    this.showMagnitudeProperty.reset();
    this.showPhaseProperty.reset();
    this.showClassicalProbabilityProperty.reset();
    this.showZerosProperty.reset();
    this.showTotalEnergyProperty.reset();
    this.showPotentialEnergyProperty.reset();
    this.superpositionTypeProperty.reset();
    this.superpositionConfigProperty.reset();
  }

  /**
   * Steps the model forward in time.
   * @param dt - The time step in seconds (can be negative for backward stepping)
   * @param forced - If true, steps even when paused (for manual stepping buttons)
   */
  public step(dt: number, forced = false): void {
    // Prevent reentry
    if (this.isStepping) {
      return;
    }

    this.isStepping = true;
    try {
      if (this.isPlayingProperty.value || forced) {
        // Convert dt to femtoseconds and apply speed multiplier (only when playing normally)
        const speedMultiplier = forced
          ? 1
          : this.timeSpeedProperty.value === TimeSpeed.SLOW
            ? BaseModel.SLOW_SPEED_MULTIPLIER
            : 1;
        const dtFemtoseconds = dt * speedMultiplier; // seconds to femtoseconds
        this.timeProperty.value += dtFemtoseconds;
        // Quantum mechanical time evolution is handled in the view layer
      }
    } finally {
      this.isStepping = false;
    }
  }

  /**
   * Calculate classical probability density from potential energy and energy level.
   * This is a common calculation used by all model classes to compute the classical
   * probability density for display alongside quantum probability.
   *
   * The classical probability density is P(x) ∝ 1/v(x) where v(x) is the classical velocity.
   * For a particle with total energy E in potential V(x):
   *   v(x) = √[2(E - V(x))/m]
   *   P(x) = 1/v(x) = √[m/(2(E - V(x)))] = 1/√[2(E - V(x))/m]
   *
   * To prevent singularities at turning points (where E ≈ V(x) and v → 0), we use a
   * minimum kinetic energy threshold (5% of total energy). This prevents display issues
   * when classical probability is plotted with quantum probability on the same scale.
   *
   * @param potential - Array of potential energy values at each grid point (in Joules)
   * @param energy - Total energy of the particle (in Joules)
   * @param mass - Particle mass (in kg)
   * @param xGrid - Array of x positions (in meters)
   * @returns Normalized classical probability density array (in 1/meters)
   */
  protected calculateClassicalProbabilityDensity(
    potential: number[],
    energy: number,
    mass: number,
    xGrid: number[],
  ): number[] {
    const classicalProbability: number[] = [];
    let integralSum = 0;

    // First pass: find maximum kinetic energy to set appropriate threshold
    let maxKE = 0;
    for (let i = 0; i < xGrid.length; i++) {
      const ke = energy - potential[i];
      if (ke > maxKE) {
        maxKE = ke;
      }
    }

    // Use minimum kinetic energy threshold to prevent singularities
    // This preserves the shape while avoiding infinities at turning points
    const minKE = BaseModel.MIN_KINETIC_ENERGY_FRACTION * maxKE;

    // Second pass: calculate probability with threshold
    for (let i = 0; i < xGrid.length; i++) {
      const kineticEnergy = energy - potential[i];

      let probability = 0;
      if (kineticEnergy > 0) {
        // Use minimum kinetic energy to prevent singularities at turning points
        // This prevents display issues when classical probability is plotted with quantum probability
        const safeKE = Math.max(kineticEnergy, minKE);
        probability = 1 / Math.sqrt((2 * safeKE) / mass);
      }
      classicalProbability.push(probability);

      // Trapezoidal integration for normalization
      if (i > 0) {
        const dx = xGrid[i] - xGrid[i - 1];
        integralSum +=
          ((probability + classicalProbability[i - 1]) * dx) /
          BaseModel.TRAPEZOIDAL_DIVISOR;
      }
    }

    // Normalize so that ∫P(x)dx = 1
    if (integralSum > 0) {
      for (let i = 0; i < classicalProbability.length; i++) {
        classicalProbability[i] /= integralSum;
      }
    }

    return classicalProbability;
  }

  /**
   * Abstract method for calculating bound states.
   * Subclasses must implement this method to perform the actual bound state calculations.
   */
  protected abstract calculateBoundStates(): void;

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
   * Calculates the energy eigenvalues for the current well parameters.
   * @param n - The quantum number (1, 2, 3, ...)
   * @returns The energy in eV
   */
  public getEnergyLevel(n: number): number {
    // Ensure bound states are calculated
    if (!this.boundStateResult) {
      this.calculateBoundStates();
    }

    // Return energy for quantum number n (1-indexed)
    if (
      this.boundStateResult &&
      n > 0 &&
      n <= this.boundStateResult.energies.length
    ) {
      const energyJoules = this.boundStateResult.energies[n - 1];
      return Schrodinger1DSolver.joulesToEV(energyJoules);
    }

    // Fallback to analytical formula if solver fails
    const L = this.wellWidthProperty.value * QuantumConstants.NM_TO_M;
    const energy =
      (n *
        n *
        Math.PI *
        Math.PI *
        QuantumConstants.HBAR *
        QuantumConstants.HBAR) /
      (2 * QuantumConstants.ELECTRON_MASS * L * L);
    return Schrodinger1DSolver.joulesToEV(energy);
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

    if (
      this.boundStateResult &&
      n > 0 &&
      n <= this.boundStateResult.wavefunctions.length
    ) {
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
      return this.boundStateResult.xGrid.map(
        (x) => x * QuantumConstants.M_TO_NM,
      );
    }

    return null;
  }

  /**
   * Calculate the classical probability density for a given energy level.
   * This method provides a default implementation that can be overridden by subclasses.
   * @param energyIndex - Index of the energy level (0-indexed)
   * @returns Array of classical probability density values, or null if unavailable
   */
  public getClassicalProbabilityDensity(_energyIndex: number): number[] | null {
    // Subclasses must override this method to provide potential-specific implementations
    return null;
  }
}
