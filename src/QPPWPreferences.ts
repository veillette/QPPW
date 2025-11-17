/**
 * Preferences for the Quantum Bound States (QPPW) Simulation
 */

import { BooleanProperty, NumberProperty, Property } from "scenerystack/axon";
import { Range } from "scenerystack";
import { Tandem, StringIO } from "scenerystack/tandem";
import qppw from "./QPPWNamespace.js";
import { NumericalMethod } from "./common/model/Schrodinger1DSolver.js";

const QPPWPreferences = {
  // Simulation Preferences

  /**
   * Whether to automatically pause the simulation when the browser tab is hidden
   */
  autoPauseWhenTabHiddenProperty: new BooleanProperty(true, {
    tandem: Tandem.PREFERENCES.createTandem("autoPauseWhenTabHiddenProperty"),
    phetioFeatured: true,
  }),

  /**
   * Numerical method for solving the Schrödinger equation
   * Options: 'numerov', 'matrix_numerov', 'dvr', 'fgh', 'spectral'
   */
  numericalMethodProperty: new Property<NumericalMethod>(NumericalMethod.FGH, {
    tandem: Tandem.PREFERENCES.createTandem("numericalMethodProperty"),
    phetioFeatured: true,
    phetioValueType: StringIO,
    validValues: [NumericalMethod.NUMEROV, NumericalMethod.MATRIX_NUMEROV, NumericalMethod.DVR, NumericalMethod.FGH, NumericalMethod.SPECTRAL],
  }),

  /**
   * Number of grid points for numerical solvers
   * Range: 64-2000 points. Higher values give more accurate results but slower computation.
   * Default: 128 points (fast computation for interactive exploration)
   */
  gridPointsProperty: new NumberProperty(128, {
    tandem: Tandem.PREFERENCES.createTandem("gridPointsProperty"),
    phetioFeatured: true,
    range: new Range(64, 2000),
  }),

  // Visual Preferences

  /**
   * Whether the user has reduced motion enabled in their OS settings
   */
  reducedMotionProperty: new BooleanProperty(
    window.matchMedia("(prefers-reduced-motion: reduce)").matches,
    {
      tandem: Tandem.PREFERENCES.createTandem("reducedMotionProperty"),
      phetioFeatured: true,
    },
  ),

  /**
   * Whether high contrast mode is enabled for improved accessibility
   */
  highContrastModeProperty: new BooleanProperty(false, {
    tandem: Tandem.PREFERENCES.createTandem("highContrastModeProperty"),
    phetioFeatured: true,
  }),

  // Audio/Voicing Preferences

  /**
   * Whether to announce parameter changes (e.g., well depth, width, etc.) via voicing
   */
  announceParameterChangesProperty: new BooleanProperty(false, {
    tandem: Tandem.PREFERENCES.createTandem("announceParameterChangesProperty"),
  }),

  /**
   * Whether to announce state changes (e.g., play, pause, reset, time speed) via voicing
   */
  announceStateChangesProperty: new BooleanProperty(false, {
    tandem: Tandem.PREFERENCES.createTandem("announceStateChangesProperty"),
  }),

  /**
   * Whether to announce drag interactions with objects via voicing
   */
  announceDragInteractionsProperty: new BooleanProperty(false, {
    tandem: Tandem.PREFERENCES.createTandem("announceDragInteractionsProperty"),
  }),
};

qppw.register("QPPWPreferences", QPPWPreferences);

export default QPPWPreferences;
