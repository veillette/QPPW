/**
 * IntroScreenView is the view for the Intro screen.
 * It displays quantum wells without play/pause controls or superposition options.
 */

import { ScreenViewOptions } from "scenerystack/sim";
import { Node } from "scenerystack/scenery";
import { TReadOnlyProperty } from "scenerystack/axon";
import {
  BaseScreenView,
  ScreenStringProperties,
} from "../../common/view/BaseScreenView.js";
import { IntroModel } from "../model/IntroModel.js";
import { EnergyChartNode } from "../../common/view/EnergyChartNode.js";
import { WaveFunctionChartNode } from "../../common/view/WaveFunctionChartNode.js";
import { IntroControlPanelNode } from "./IntroControlPanelNode.js";
import stringManager from "../../i18n/StringManager.js";

export class IntroScreenView extends BaseScreenView {
  private introControlPanel: IntroControlPanelNode;

  public constructor(model: IntroModel, options?: ScreenViewOptions) {
    super(model, options);

    // Calculate layout dimensions
    const margin = 10;

    // Fixed chart dimensions
    const chartsWidth = 600;
    const energyChartHeight = 300;
    const waveFunctionChartHeight = 180;

    // Create the energy chart (top plot)
    // Cast to OneWellModel since IntroModel has all necessary properties
    this.energyChart = new EnergyChartNode(
      model as unknown as import("../../one-well/model/OneWellModel.js").OneWellModel,
      {
        width: chartsWidth,
        height: energyChartHeight,
      },
    );

    // Create the wave function chart (bottom plot)
    this.waveFunctionChart = new WaveFunctionChartNode(
      model as unknown as import("../../one-well/model/OneWellModel.js").OneWellModel,
      {
        width: chartsWidth,
        height: waveFunctionChartHeight,
      },
    );

    // Position charts
    this.energyChart.left = margin;
    this.energyChart.top = 10;

    this.waveFunctionChart.left = margin;
    this.waveFunctionChart.top = margin + energyChartHeight + 30;

    // Create listbox parent node for ComboBox popups
    this.listBoxParent = new Node();

    // Create control panel (simplified for intro)
    this.introControlPanel = new IntroControlPanelNode(
      model,
      this.listBoxParent,
      this.waveFunctionChart,
    );
    this.introControlPanel.left = chartsWidth + margin * 2;
    this.introControlPanel.top = margin;

    // Add all components to the view
    this.addChild(this.energyChart);
    this.addChild(this.waveFunctionChart);
    this.addChild(this.introControlPanel);
    this.addChild(this.listBoxParent); // ListBox parent must be added last for proper z-ordering
  }

  /**
   * Get screen-specific string properties for creating dialog content.
   */
  protected getScreenStringProperties(): ScreenStringProperties {
    return {
      titleStringProperty: stringManager.introStringProperty,
      descriptionStringProperty: stringManager.introDescriptionStringProperty,
      keyConceptsStringProperty: stringManager.introKeyConceptsStringProperty,
      interactionsStringProperty: stringManager.introInteractionsStringProperty,
      educationalContentStringProperty:
        stringManager.introEducationalContentStringProperty,
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
   * Steps the screen view forward in time.
   * Note: Intro screen doesn't have play/pause controls, but model still steps
   * @param dt - The time step in seconds
   */
  public override step(dt: number): void {
    super.step(dt);
    if (this.model) {
      this.model.step(dt);
    }
  }
}
