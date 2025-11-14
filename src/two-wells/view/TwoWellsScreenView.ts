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

    const keyConceptsTitle = new Text("Key Concepts:", {
      font: new PhetFont({ size: 14, weight: "bold" }),
      fill: QPPWColors.textFillProperty,
    });

    const keyConceptsList = new RichText(
      "• Quantum tunneling through potential barriers<br>" +
      "• Symmetric and antisymmetric states in double wells<br>" +
      "• Energy level splitting due to coupling<br>" +
      "• Barrier penetration and transmission probability<br>" +
      "• Bonding and antibonding states",
      {
        font: new PhetFont(13),
        fill: QPPWColors.textFillProperty,
        maxWidth: 500,
      }
    );

    const interactionTitle = new Text("Interactions:", {
      font: new PhetFont({ size: 14, weight: "bold" }),
      fill: QPPWColors.textFillProperty,
    });

    const interactionsList = new RichText(
      "• Adjust barrier height to control tunneling rate<br>" +
      "• Change barrier width to see its effect on coupling<br>" +
      "• Modify well separation to observe state splitting<br>" +
      "• Compare energy levels of coupled vs. isolated wells<br>" +
      "• Observe wave function localization and delocalization",
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
      stringManager.twoWellsSummaryStringProperty.value +
      "<br><br>" +
      "<strong>Screen Overview:</strong><br>" +
      "This screen demonstrates quantum tunneling in a double potential well system. " +
      "Two quantum wells are separated by a tunable potential barrier, allowing particles to tunnel between them. " +
      "The visualization shows how barrier properties affect energy level splitting and wave function coupling.<br><br>" +
      "<strong>Available Controls:</strong><br>" +
      "• Barrier Height slider: Controls the energy barrier between wells<br>" +
      "• Barrier Width slider: Adjusts the thickness of the potential barrier<br>" +
      "• Well Separation slider: Changes the distance between the two wells<br>" +
      "• Quantum State selector: Choose different energy eigenstates<br>" +
      "• Play/Pause/Step controls: Control time evolution<br>" +
      "• Reset button: Returns all parameters to initial values<br><br>" +
      "<strong>Learning Objectives:</strong><br>" +
      "Understand quantum tunneling as a uniquely quantum mechanical phenomenon. " +
      "Observe how barrier properties affect tunneling probability and coupling strength. " +
      "Learn about bonding and antibonding states in coupled quantum systems. " +
      "Explore the relationship between well separation and energy level splitting.",
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
