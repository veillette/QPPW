/**
 * OneWellScreenView is the main view for the One Well screen.
 * It displays a single quantum potential well and its wave functions.
 */

import {
  BaseScreenView,
  ScreenStringProperties,
} from "../../common/view/BaseScreenView.js";
import { OneWellModel } from "../model/OneWellModel.js";
import { ScreenViewOptions } from "scenerystack/sim";
import { TReadOnlyProperty } from "scenerystack/axon";
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

    this.createStandardLayout(model, {
      allowedPotentialTypes: allowedPotentials,
    });
  }

  /**
   * Get screen-specific string properties for creating dialog content.
   */
  protected getScreenStringProperties(): ScreenStringProperties {
    return {
      titleStringProperty: stringManager.oneWellStringProperty,
      descriptionStringProperty: stringManager.oneWellDescriptionStringProperty,
      keyConceptsStringProperty: stringManager.oneWellKeyConceptsStringProperty,
      interactionsStringProperty:
        stringManager.oneWellInteractionsStringProperty,
      educationalContentStringProperty:
        stringManager.oneWellEducationalContentStringProperty,
    };
  }

  /**
   * Get the common "Key Concepts" title string property.
   */
  protected getKeyConceptsTitleStringProperty(): TReadOnlyProperty<string> {
    return stringManager.keyConceptsTitleStringProperty;
  }

  /**
   * Get the common "Interactions" title string property.
   */
  protected getInteractionsTitleStringProperty(): TReadOnlyProperty<string> {
    return stringManager.interactionsTitleStringProperty;
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
