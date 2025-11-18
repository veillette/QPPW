/**
 * OneWellScreenView is the main view for the One Well screen.
 * It displays a single quantum potential well and its wave functions.
 */

import { BaseScreenView } from "../../common/view/BaseScreenView.js";
import { OneWellModel } from "../model/OneWellModel.js";
import { ScreenViewOptions } from "scenerystack/sim";
import { Node, Text, VBox, RichText } from "scenerystack/scenery";
import { PhetFont } from "scenerystack/scenery-phet";
import QPPWColors from "../../QPPWColors.js";
import stringManager from "../../i18n/StringManager.js";
import { PotentialType } from "../../common/model/PotentialFunction.js";

export class OneWellScreenView extends BaseScreenView {
  public constructor(model: OneWellModel, options?: ScreenViewOptions) {
    super(model, options);

    // Create the standard quantum well layout, excluding Double Square Well
    // which is meant for the Two Wells screen
    const allowedPotentials = [
      PotentialType.INFINITE_WELL,
      PotentialType.FINITE_WELL,
      PotentialType.HARMONIC_OSCILLATOR,
      PotentialType.MORSE,
      PotentialType.POSCHL_TELLER,
      PotentialType.ROSEN_MORSE,
      PotentialType.ECKART,
      PotentialType.ASYMMETRIC_TRIANGLE,
      PotentialType.TRIANGULAR,
      PotentialType.COULOMB_1D,
      PotentialType.COULOMB_3D,
    ];

    this.createStandardLayout(model, { allowedPotentialTypes: allowedPotentials });
  }

  /**
   * Creates the content for the info dialog.
   */
  public createInfoDialogContent(): Node {
    const titleText = new Text(stringManager.oneWellStringProperty, {
      font: new PhetFont({ size: 18, weight: "bold" }),
      fill: QPPWColors.textFillProperty,
    });

    const descriptionText = new RichText(stringManager.oneWellDescriptionStringProperty, {
      font: new PhetFont(14),
      fill: QPPWColors.textFillProperty,
      maxWidth: 500,
    });

    const keyConceptsTitle = new Text(stringManager.keyConceptsTitleStringProperty, {
      font: new PhetFont({ size: 14, weight: "bold" }),
      fill: QPPWColors.textFillProperty,
    });

    const keyConceptsList = new RichText(
      stringManager.oneWellKeyConceptsStringProperty,
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
      stringManager.oneWellInteractionsStringProperty,
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
      stringManager.oneWellEducationalContentStringProperty,
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
    if (this.model) {
      (this.model as OneWellModel).step(dt);
    }
  }
}
