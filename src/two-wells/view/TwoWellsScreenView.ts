/**
 * TwoWellsScreenView is the main view for the Two Wells screen.
 * It displays a double quantum potential well and demonstrates quantum tunneling.
 */

import {
  BaseScreenView,
  ScreenStringProperties,
} from "../../common/view/BaseScreenView.js";
import { TwoWellsModel } from "../model/TwoWellsModel.js";
import { PotentialType } from "../../common/model/PotentialFunction.js";
import { ScreenViewOptions } from "scenerystack/sim";
import { TReadOnlyProperty } from "scenerystack/axon";
import stringManager from "../../i18n/StringManager.js";

export class TwoWellsScreenView extends BaseScreenView {
  public constructor(model: TwoWellsModel, options?: ScreenViewOptions) {
    super(model, options);

    // Create the standard quantum well layout with custom control panel options
    // - Hide particle mass slider
    // - Allow Coulomb 1D and Double Square Well potential types
    this.createStandardLayout(model, {
      showParticleMass: false,
      allowedPotentialTypes: [
        PotentialType.COULOMB_1D,
        PotentialType.DOUBLE_SQUARE_WELL,
      ],
    });
  }

  /**
   * Get screen-specific string properties for creating dialog content.
   */
  protected getScreenStringProperties(): ScreenStringProperties {
    return {
      titleStringProperty: stringManager.twoWellsStringProperty,
      descriptionStringProperty:
        stringManager.twoWellsDescriptionStringProperty,
      keyConceptsStringProperty:
        stringManager.twoWellsKeyConceptsStringProperty,
      interactionsStringProperty:
        stringManager.twoWellsInteractionsStringProperty,
      educationalContentStringProperty:
        stringManager.twoWellsEducationalContentStringProperty,
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
    (this.model as TwoWellsModel).step(dt);
  }
}
