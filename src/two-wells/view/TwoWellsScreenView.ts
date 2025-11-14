/**
 * TwoWellsScreenView is the main view for the Two Wells screen.
 * It displays a double quantum potential well and demonstrates quantum tunneling.
 */

import { BaseScreenView } from "../../common/view/BaseScreenView.js";
import { TwoWellsModel } from "../model/TwoWellsModel.js";
import { TwoWellsControlPanelNode } from "./TwoWellsControlPanelNode.js";
import { ScreenViewOptions } from "scenerystack/sim";
import { Node, Text } from "scenerystack/scenery";
import { PhetFont } from "scenerystack/scenery-phet";
import QPPWColors from "../../QPPWColors.js";
import stringManager from "../../i18n/StringManager.js";

export class TwoWellsScreenView extends BaseScreenView {
  public constructor(model: TwoWellsModel, options?: ScreenViewOptions) {
    super(model, options);

    // Create the listBoxParent for ComboBox popup BEFORE creating control panel
    this.listBoxParent = new Node();

    // Create the custom control panel for Two Wells (without mass slider, limited potential types)
    const customControlPanel = new TwoWellsControlPanelNode(model, this.listBoxParent);

    // Create the standard quantum well layout with custom control panel
    this.createStandardLayout(model, customControlPanel);
  }

  /**
   * Creates the content for the info dialog.
   */
  public createInfoDialogContent(): Node {
    return new Text(stringManager.twoWellsDescriptionStringProperty, {
      font: new PhetFont(14),
      fill: QPPWColors.textFillProperty,
    });
  }

  /**
   * Creates the screen summary content for accessibility.
   */
  public createScreenSummaryContent(): Node {
    return new Text(stringManager.twoWellsSummaryStringProperty, {
      font: new PhetFont(14),
      fill: QPPWColors.textFillProperty,
    });
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
    (this.model as TwoWellsModel).step(dt);
  }
}
