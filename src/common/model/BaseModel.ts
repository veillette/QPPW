/**
 * BaseModel is an abstract base class for all physics models in QPPW.
 * It provides common functionality for time evolution, simulation control, and solver configuration.
 */

import { EnumerationProperty, NumberProperty, Property } from "scenerystack/axon";
import { TimeSpeed } from "scenerystack";
import Schrodinger1DSolver, { NumericalMethod } from "./Schrodinger1DSolver.js";
import QPPWPreferences from "../../QPPWPreferences.js";

export type SimulationSpeed = "normal" | "fast";

export abstract class BaseModel {
  // Simulation state properties
  public readonly isPlayingProperty: Property<boolean>;
  public readonly timeProperty: NumberProperty; // In femtoseconds
  public readonly timeSpeedProperty: EnumerationProperty<TimeSpeed>;

  // Solver for quantum calculations
  protected readonly solver: Schrodinger1DSolver;

  protected constructor() {
    // Initialize simulation state
    this.isPlayingProperty = new Property<boolean>(false);
    this.timeProperty = new NumberProperty(0); // in femtoseconds
    this.timeSpeedProperty = new EnumerationProperty(TimeSpeed.NORMAL);

    // Initialize solver with user's preferred method
    this.solver = new Schrodinger1DSolver();

    // Update solver method when preference changes
    QPPWPreferences.numericalMethodProperty.link((method) => {
      this.solver.setNumericalMethod(method);
      this.onSolverMethodChanged(method);
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
  }

  /**
   * Steps the model forward in time.
   * @param dt - The time step in seconds (can be negative for backward stepping)
   * @param forced - If true, steps even when paused (for manual stepping buttons)
   */
  public step(dt: number, forced = false): void {
    if (this.isPlayingProperty.value || forced) {
      // Convert dt to femtoseconds and apply speed multiplier (only when playing normally)
      const speedMultiplier = forced ? 1 : (this.timeSpeedProperty.value === TimeSpeed.SLOW ? 1/10 : 1);
      const dtFemtoseconds = (dt) * speedMultiplier; // seconds to femtoseconds
      this.timeProperty.value += dtFemtoseconds;
      // Quantum mechanical time evolution is handled in the view layer
    }
  }
}
