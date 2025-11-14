/**
 * OneWellScreenView is the main view for the One Well screen.
 * It displays a single quantum potential well and its wave functions.
 */

import { BaseScreenView } from "../../common/view/BaseScreenView.js";
import { OneWellModel } from "../model/OneWellModel.js";
import { ScreenViewOptions } from "scenerystack/sim";
import { Node, Text } from "scenerystack/scenery";
import { PhetFont } from "scenerystack/scenery-phet";
import QPPWColors from "../../QPPWColors.js";
import stringManager from "../../i18n/StringManager.js";

export class OneWellScreenView extends BaseScreenView {
  public constructor(model: OneWellModel, options?: ScreenViewOptions) {
    super(model, options);

    // Create the standard quantum well layout
    this.createStandardLayout(model);
  }

  /**
   * Creates the content for the info dialog.
   */
  public createInfoDialogContent(): Node {
    const text = new Text(stringManager.oneWellDescriptionStringProperty, {
      font: new PhetFont(14),
      fill: QPPWColors.textFillProperty,
    });
    return new Node({ children: [text] });
  }

  /**
   * Creates the screen summary content for accessibility.
   */
  public createScreenSummaryContent(): Node {
    const text = new Text(stringManager.oneWellSummaryStringProperty, {
      font: new PhetFont(14),
      fill: QPPWColors.textFillProperty,
    });
    return new Node({ children: [text] });
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
