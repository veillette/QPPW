/**
 * IntroViewState holds display-related properties specific to the Intro screen.
 */

import { Property } from "scenerystack/axon";
import { BaseViewState } from "../../common/view/BaseViewState.js";

// Intro screen only supports these two modes (no "phaseColor")
export type IntroDisplayMode = "probabilityDensity" | "waveFunction";

export class IntroViewState extends BaseViewState {
  // Intro-specific display mode (simplified, no "phaseColor" mode)
  public readonly displayModeProperty: Property<IntroDisplayMode>;

  public constructor() {
    super();

    // Initialize intro-specific display settings
    this.displayModeProperty = new Property<IntroDisplayMode>("probabilityDensity");
  }

  /**
   * Resets all view state properties to their initial values.
   */
  public override reset(): void {
    super.reset();
    this.displayModeProperty.reset();
  }
}
