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
import { TwoWellsViewState } from "./TwoWellsViewState.js";

export class TwoWellsScreenView extends BaseScreenView {
  private readonly viewState: TwoWellsViewState;

  public constructor(model: TwoWellsModel, options?: ScreenViewOptions) {
    super(model, options);

    // Create the view state for display properties
    this.viewState = new TwoWellsViewState();

    // Create the standard quantum well layout with custom control panel options
    // - Hide particle mass slider
    // - Allow Coulomb 1D and Double Square Well potential types
    this.createStandardLayout(model, this.viewState, {
      showParticleMass: false,
      allowedPotentialTypes: [
        PotentialType.COULOMB_1D,
        PotentialType.DOUBLE_SQUARE_WELL,
      ],
    });

    // Set up PDOM (Parallel DOM) structure for accessibility
    this.setupAccessibility(model);
  }

  /**
   * Sets up the three-section PDOM structure for accessibility.
   */
  private setupAccessibility(_model: TwoWellsModel): void {
    // Create screen summary (Section 1)
    this.screenSummaryNode = this.createScreenSummaryNode({
      screenName: "Two Wells",
      screenDescription:
        "Two Wells screen for exploring quantum tunneling and energy level splitting in double potential wells.",
    });

    // Create play area node (Section 2) and add charts to it
    this.playAreaNode = this.createPlayAreaNode();
    if (this.chartsContainer) {
      this.playAreaNode.addChild(this.chartsContainer);
    }

    // Create control area node (Section 3) and add controls to it
    this.controlAreaNode = this.createControlAreaNode();
    if (this.controlPanel) {
      this.controlAreaNode.addChild(this.controlPanel);
    }
    if (this.simulationControlBar) {
      this.controlAreaNode.addChild(this.simulationControlBar);
    }

    // Set PDOM navigation order
    this.setupPDOMStructure();
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
    this.viewState.reset();
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
