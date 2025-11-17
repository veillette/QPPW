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

    const keyConceptsTitle = new Text("Key Concepts:", {
      font: new PhetFont({ size: 14, weight: "bold" }),
      fill: QPPWColors.textFillProperty,
    });

    const keyConceptsList = new RichText(
      "• Quantum confinement in a single potential well<br>" +
      "• Discrete energy levels and eigenstates<br>" +
      "• Wave function normalization and probability density<br>" +
      "• Energy quantization in bound states",
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
      "• Adjust well width to see how energy levels change<br>" +
      "• Select different quantum states (n=1,2,3...)<br>" +
      "• Observe the wave function and probability density<br>" +
      "• Compare energy levels on the energy chart",
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
      stringManager.oneWellSummaryStringProperty.value +
      "<br><br>" +
      "<strong>Screen Overview:</strong><br>" +
      "This screen simulates a single quantum potential well. " +
      "The top chart displays energy levels and the potential well shape. " +
      "The bottom chart shows the wave function and probability density for the selected quantum state.<br><br>" +
      "<strong>Available Controls:</strong><br>" +
      "• Well Width slider: Adjusts the width of the potential well<br>" +
      "• Quantum State selector: Choose between different energy eigenstates (n=1, 2, 3, etc.)<br>" +
      "• Visualization options: Toggle between different display modes<br>" +
      "• Play/Pause/Step controls: Control the time evolution of the wave function<br>" +
      "• Reset button: Returns all parameters to their initial values<br><br>" +
      "<strong>Learning Objectives:</strong><br>" +
      "Explore how particle confinement leads to quantized energy levels. " +
      "Observe how changing the well width affects the energy spectrum. " +
      "Understand the relationship between energy levels and quantum numbers.",
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
