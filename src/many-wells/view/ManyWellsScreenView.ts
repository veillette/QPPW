/**
 * ManyWellsScreenView is the main view for the Many Wells screen.
 * It displays multiple quantum potential wells (1-10) with two potential types:
 * - Multi-square well (generalization of double square well)
 * - Multi-Coulomb 1D (multiple Coulomb centers)
 */

import { BaseScreenView } from "../../common/view/BaseScreenView.js";
import { ManyWellsModel } from "../model/ManyWellsModel.js";
import { PotentialType } from "../../common/model/PotentialFunction.js";
import { ScreenViewOptions } from "scenerystack/sim";
import { Node, Text, VBox, RichText } from "scenerystack/scenery";
import { PhetFont } from "scenerystack/scenery-phet";
import QPPWColors from "../../QPPWColors.js";
import stringManager from "../../i18n/StringManager.js";

export class ManyWellsScreenView extends BaseScreenView {
  public constructor(model: ManyWellsModel, options?: ScreenViewOptions) {
    super(model, options);

    // Create the standard quantum well layout with custom control panel options
    // - Hide particle mass slider (use electron mass for simplicity)
    // - Allow Multi-Square Well and Multi-Coulomb 1D potential types
    // - Show number of wells slider (1-10)
    this.createStandardLayout(model, {
      showParticleMass: false,
      allowedPotentialTypes: [
        PotentialType.MULTI_SQUARE_WELL,
        PotentialType.MULTI_COULOMB_1D,
      ],
    });
  }

  /**
   * Creates the content for the info dialog.
   */
  public createInfoDialogContent(): Node {
    const titleText = new Text(stringManager.manyWellsStringProperty, {
      font: new PhetFont({ size: 18, weight: "bold" }),
      fill: QPPWColors.textFillProperty,
    });

    const descriptionText = new RichText(
      stringManager.manyWellsDescriptionStringProperty,
      {
        font: new PhetFont(14),
        fill: QPPWColors.textFillProperty,
        maxWidth: 500,
      },
    );

    const keyConceptsTitle = new Text(
      stringManager.keyConceptsTitleStringProperty,
      {
        font: new PhetFont({ size: 14, weight: "bold" }),
        fill: QPPWColors.textFillProperty,
      },
    );

    const keyConceptsList = new RichText(
      stringManager.manyWellsKeyConceptsStringProperty,
      {
        font: new PhetFont(13),
        fill: QPPWColors.textFillProperty,
        maxWidth: 500,
      },
    );

    const interactionTitle = new Text(
      stringManager.interactionsTitleStringProperty,
      {
        font: new PhetFont({ size: 14, weight: "bold" }),
        fill: QPPWColors.textFillProperty,
      },
    );

    const interactionsList = new RichText(
      stringManager.manyWellsInteractionsStringProperty,
      {
        font: new PhetFont(13),
        fill: QPPWColors.textFillProperty,
        maxWidth: 500,
      },
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
      stringManager.manyWellsEducationalContentStringProperty,
      {
        font: new PhetFont(13),
        fill: QPPWColors.textFillProperty,
        maxWidth: 600,
      },
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
    (this.model as ManyWellsModel).step(dt);
  }
}
