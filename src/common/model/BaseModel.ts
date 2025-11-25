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
}
