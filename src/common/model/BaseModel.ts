/**
 * BaseModel is an abstract base class for all physics models in QPPW.
 * It provides common functionality for time evolution, simulation control, and solver configuration.
 */

import {
  EnumerationProperty,
  NumberProperty,
  Property,
} from "scenerystack/axon";
import { TimeSpeed } from "scenerystack";
import Schrodinger1DSolver, { NumericalMethod } from "./Schrodinger1DSolver.js";
import QPPWPreferences from "../../QPPWPreferences.js";

export abstract class BaseModel {
  // Default time step for manual stepping (in seconds, ~1 frame at 60 FPS)
  public static readonly MANUAL_STEP_SIZE = 0.016;

  // Simulation state properties
  public readonly isPlayingProperty: Property<boolean>;
  public readonly timeProperty: NumberProperty; // In femtoseconds
  public readonly timeSpeedProperty: EnumerationProperty<TimeSpeed>;

  // Display options
  public readonly showClassicalProbabilityProperty: Property<boolean>;

  // Solver for quantum calculations
  protected readonly solver: Schrodinger1DSolver;

  // Guard flag to prevent reentry in step method
  private isStepping: boolean = false;

  protected constructor() {
    // Initialize simulation state
    this.isPlayingProperty = new Property<boolean>(false);
    this.timeProperty = new NumberProperty(0); // in femtoseconds
    this.timeSpeedProperty = new EnumerationProperty(TimeSpeed.NORMAL);

    // Initialize display options
    this.showClassicalProbabilityProperty = new Property<boolean>(false);

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
  }

  /**
   * Called when the solver method changes.
   * Subclasses should override this to invalidate caches or recalculate results.
   * @param method - The new numerical method
   */
  protected abstract onSolverMethodChanged(method: NumericalMethod): void;

  /**
   * Resets the model to its initial state.
   * Subclasses should override this method and call super.resetAll().
   */
  public resetAll(): void {
    this.isPlayingProperty.reset();
    this.timeProperty.reset();
    this.timeSpeedProperty.reset();
    this.showClassicalProbabilityProperty.reset();
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
            ? 1 / 10
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

    // Use 1% of maximum kinetic energy as threshold to prevent singularities
    // This preserves the shape while avoiding infinities at turning points
    const minKE = 0.01 * maxKE;

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
        integralSum += (probability + classicalProbability[i - 1]) * dx / 2;
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
}
