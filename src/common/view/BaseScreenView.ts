/**
 * BaseScreenView is an abstract base class for all screen views in the QPPW simulation.
 * It provides common functionality and requires subclasses to implement screen-specific methods.
 */

import { ScreenView, ScreenViewOptions } from "scenerystack/sim";
import { ResetAllButton } from "scenerystack/scenery-phet";
import { Node } from "scenerystack/scenery";
import QPPWColors from "../../QPPWColors.js";

export abstract class BaseScreenView extends ScreenView {
  protected readonly resetAllButton: ResetAllButton;

  protected constructor(resetCallback: () => void, options?: ScreenViewOptions) {
    super(options);

    // Create the reset all button in the bottom-right corner
    this.resetAllButton = new ResetAllButton({
      listener: resetCallback,
      right: this.layoutBounds.maxX - 10,
      bottom: this.layoutBounds.maxY - 10,
      baseColor: QPPWColors.panelFillProperty,
    });
    this.addChild(this.resetAllButton);
  }

  /**
   * Creates the content for the info dialog.
   * Subclasses must implement this method to provide screen-specific information.
   */
  public abstract createInfoDialogContent(): Node;

  /**
   * Creates the screen summary content for accessibility.
   * Subclasses must implement this method to provide screen-specific summary.
   */
  public abstract createScreenSummaryContent(): Node;

  /**
   * Resets the screen view to its initial state.
   * Subclasses should override this method to add screen-specific reset logic.
   */
  public reset(): void {
    // Base implementation - subclasses should call super.reset() and add their own logic
  }

  /**
   * Steps the screen view forward in time.
   * Subclasses should override this method to add screen-specific step logic.
   * @param dt - The time step in seconds
   */
  public step(dt: number): void {
    // Base implementation - subclasses should override
  }
}
