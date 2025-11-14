/**
 * OneWellScreenView is the main view for the One Well screen.
 * It displays a single quantum potential well and its wave functions.
 */

import { BaseScreenView } from "../../common/view/BaseScreenView.js";
import { OneWellModel } from "../model/OneWellModel.js";
import { ScreenViewOptions } from "scenerystack/sim";

export class OneWellScreenView extends BaseScreenView {
  public constructor(model: OneWellModel, options?: ScreenViewOptions) {
    super(
      model,
      () => {
        model.resetAll();
        this.reset();
      },
      options,
    );

    // Create the standard quantum well layout
    this.createStandardLayout(model);
  }

  /**
   * Resets the screen view to its initial state.
   */
  public override reset(): void {
    super.reset();
    // Add screen-specific reset logic here
  }

  /**
   * Steps the screen view forward in time.
   * @param dt - The time step in seconds
   */
  public override step(dt: number): void {
    super.step(dt);
    if (this.model) {
      (this.model as OneWellModel).step(dt);
    }
  }
}
