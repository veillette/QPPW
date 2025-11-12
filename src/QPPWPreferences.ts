/**
 * Preferences for the Quantum Bound States (QPPW) Simulation
 */

import { BooleanProperty } from "scenerystack/axon";
import { Tandem } from "scenerystack/tandem";
import qppw from "./QPPWNamespace.js";

const QPPWPreferences = {
  // Simulation Preferences

  /**
   * Whether to automatically pause the simulation when the browser tab is hidden
   */
  autoPauseWhenTabHiddenProperty: new BooleanProperty(true, {
    tandem: Tandem.PREFERENCES.createTandem("autoPauseWhenTabHiddenProperty"),
    phetioFeatured: true,
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
