/**
 * TwoWellsScreenView is the main view for the Two Wells screen.
 * It displays a double quantum potential well and demonstrates quantum tunneling.
 */

import { BaseScreenView } from "../../common/view/BaseScreenView.js";
import { TwoWellsModel } from "../model/TwoWellsModel.js";
import { PotentialType } from "../../common/model/PotentialFunction.js";
import { ScreenViewOptions } from "scenerystack/sim";
import { Node, Text, VBox, RichText } from "scenerystack/scenery";
import { PhetFont } from "scenerystack/scenery-phet";
import QPPWColors from "../../QPPWColors.js";
import stringManager from "../../i18n/StringManager.js";

export class TwoWellsScreenView extends BaseScreenView {
  public constructor(model: TwoWellsModel, options?: ScreenViewOptions) {
    super(model, options);

    // Create the standard quantum well layout with custom control panel options
    // - Hide particle mass slider
    // - Allow Coulomb 1D and Double Square Well potential types
    this.createStandardLayout(model, {
      showParticleMass: false,
      allowedPotentialTypes: [PotentialType.COULOMB_1D, PotentialType.DOUBLE_SQUARE_WELL],
    });
  }

  /**
   * Creates the content for the info dialog.
   */
  public createInfoDialogContent(): Node {
    const titleText = new Text(stringManager.twoWellsStringProperty, {
      font: new PhetFont({ size: 18, weight: "bold" }),
      fill: QPPWColors.textFillProperty,
    });

    const descriptionText = new RichText(stringManager.twoWellsDescriptionStringProperty, {
      font: new PhetFont(14),
      fill: QPPWColors.textFillProperty,
      maxWidth: 500,
    });

    const keyConceptsTitle = new Text(stringManager.keyConceptsTitleStringProperty, {
      font: new PhetFont({ size: 14, weight: "bold" }),
      fill: QPPWColors.textFillProperty,
    });

    const keyConceptsList = new RichText(
      stringManager.twoWellsKeyConceptsStringProperty,
      {
        font: new PhetFont(13),
        fill: QPPWColors.textFillProperty,
        maxWidth: 500,
      }
    );

    const interactionTitle = new Text(stringManager.interactionsTitleStringProperty, {
      font: new PhetFont({ size: 14, weight: "bold" }),
      fill: QPPWColors.textFillProperty,
    });

    const interactionsList = new RichText(
      stringManager.twoWellsInteractionsStringProperty,
      {
        font: new PhetFont(13),
        fill: QPPWColors.textFillProperty,
        maxWidth: 500,
      }
    );

    return new VBox({
      spacing: 12,
      align: "left",
      children: [
        titleText,
        descriptionText,
        keyConceptsTitle,
        keyConceptsList,
        interactionTitle,
        interactionsList,
      ],
    });
  }

  /**
   * Creates the screen summary content for accessibility.
   */
  public createScreenSummaryContent(): Node {
    const summaryText = new RichText(
      stringManager.twoWellsEducationalContentStringProperty,
      {
        font: new PhetFont(13),
        fill: QPPWColors.textFillProperty,
        maxWidth: 600,
      }
    );

    return new VBox({
      spacing: 10,
      align: "left",
      children: [summaryText],
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
