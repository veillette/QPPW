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
import { WavenumberChartNode } from "../../common/view/WavenumberChartNode.js";
import { IntroControlPanelNode } from "./IntroControlPanelNode.js";
import stringManager from "../../i18n/StringManager.js";
import { IntroViewState } from "./IntroViewState.js";

export class IntroScreenView extends BaseScreenView {
  private readonly viewState: IntroViewState;
  private introControlPanel: IntroControlPanelNode;
  private probabilityChart: WaveFunctionChartNode;
  private wavenumberChart: WavenumberChartNode;

  public constructor(model: IntroModel, options?: ScreenViewOptions) {
    super(model, options);

    // Create the view state for display properties
    this.viewState = new IntroViewState();

    // Calculate layout dimensions
    const margin = 10;
    const chartSpacing = 15;

    // Fixed chart dimensions
    const chartsWidth = 600;
    const energyChartHeight = 200;
    const probabilityChartHeight = 150;
    const waveFunctionChartHeight = 150;
    const wavenumberChartHeight = 140;

    // Create the energy chart (top plot)
    this.energyChart = new EnergyChartNode(model, this.viewState, {
      width: chartsWidth,
      height: energyChartHeight,
    });

    // Create the probability density chart (middle plot) - always shows probability density
    this.probabilityChart = new WaveFunctionChartNode(model, this.viewState, {
      width: chartsWidth,
      height: probabilityChartHeight,
      fixedDisplayMode: "probabilityDensity",
    });

    // Create the wave function chart - always shows wavefunction
    this.waveFunctionChart = new WaveFunctionChartNode(model, this.viewState, {
      width: chartsWidth,
      height: waveFunctionChartHeight,
      fixedDisplayMode: "waveFunction",
    });

    // Create the wavenumber chart (bottom plot) - shows |φ(k)|²
    this.wavenumberChart = new WavenumberChartNode(model, {
      width: chartsWidth,
      height: wavenumberChartHeight,
    });

    // Position charts stacked vertically
    this.energyChart.left = margin;
    this.energyChart.top = 10;

    this.probabilityChart.left = margin;
    this.probabilityChart.top =
      this.energyChart.top + energyChartHeight + chartSpacing;

    this.waveFunctionChart.left = margin;
    this.waveFunctionChart.top =
      this.probabilityChart.top + probabilityChartHeight + chartSpacing;

    this.wavenumberChart.left = margin;
    this.wavenumberChart.top =
      this.waveFunctionChart.top + waveFunctionChartHeight + chartSpacing;

    // Create listbox parent node for ComboBox popups
    this.listBoxParent = new Node();

    // Create control panel (simplified for intro)
    this.introControlPanel = new IntroControlPanelNode(
      model,
      this.listBoxParent,
      this.probabilityChart,
      this.waveFunctionChart,
    );
    this.introControlPanel.left = chartsWidth + margin * 2;
    this.introControlPanel.top = margin;

    // Add all components to the view
    this.addChild(this.energyChart);
    this.addChild(this.probabilityChart);
    this.addChild(this.waveFunctionChart);
    this.addChild(this.wavenumberChart);
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
   * Resets the screen view to its initial state.
   */
  public override reset(): void {
    super.reset();
    this.viewState.reset();
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
