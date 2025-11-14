/**
 * SuperpositionDialog allows users to configure superposition state amplitudes.
 */

import { Node, Text, VBox, HBox, RichText } from "scenerystack/scenery";
import { Dialog, HSlider, RectangularPushButton } from "scenerystack/sun";
import { Dimension2 } from "scenerystack/dot";
import { Property, NumberProperty } from "scenerystack/axon";
import { PhetFont } from "scenerystack/scenery-phet";
import QPPWColors from "../../QPPWColors.js";
import stringManager from "../../i18n/StringManager.js";
import { SuperpositionConfig } from "../model/SuperpositionType.js";
import { BoundStateResult } from "../model/PotentialFunction.js";

export class SuperpositionDialog {
  private readonly dialog: Dialog;
  private readonly configProperty: Property<SuperpositionConfig>;
  private readonly amplitudeProperties: NumberProperty[];

  public constructor(
    configProperty: Property<SuperpositionConfig>,
    boundStateResult: BoundStateResult | null,
  ) {
    const content = SuperpositionDialog.createContent(configProperty, boundStateResult);

    this.dialog = new Dialog(content, {
      xSpacing: 20,
      ySpacing: 20,
      cornerRadius: 8,
      fill: QPPWColors.panelFillProperty,
      stroke: QPPWColors.panelStrokeProperty,
      closeButtonColor: QPPWColors.textFillProperty,
    });

    this.configProperty = configProperty;
    this.amplitudeProperties = [];
  }

  /**
   * Shows the dialog.
   */
  public show(): void {
    this.dialog.show();
  }

  /**
   * Hides the dialog.
   */
  public hide(): void {
    this.dialog.hide();
  }

  /**
   * Creates the content for the dialog.
   */
  private static createContent(
    configProperty: Property<SuperpositionConfig>,
    boundStateResult: BoundStateResult | null,
  ): Node {
    const titleText = new Text(stringManager.superpositionDialogTitleStringProperty, {
      font: new PhetFont({ size: 18, weight: "bold" }),
      fill: QPPWColors.textFillProperty,
    });

    const descriptionText = new RichText(
      "Adjust the amplitude of each eigenstate in the superposition.<br>The sum of squared amplitudes should equal 1.",
      {
        font: new PhetFont(12),
        fill: QPPWColors.textFillProperty,
        maxWidth: 400,
      },
    );

    // Determine how many states to show based on the current config
    const config = configProperty.value;
    const numStates = Math.min(config.amplitudes.length, boundStateResult?.energies.length || 5);

    // Create sliders for each amplitude
    const sliderNodes: Node[] = [];
    const amplitudeProperties: NumberProperty[] = [];

    for (let i = 0; i < numStates; i++) {
      const amplitude = config.amplitudes[i] || 0;
      const amplitudeProperty = new NumberProperty(amplitude, {
        range: { min: 0, max: 1 },
      });
      amplitudeProperties.push(amplitudeProperty);

      const stateLabel = new RichText(`ψ<sub>${i}</sub>:`, {
        font: new PhetFont(14),
        fill: QPPWColors.textFillProperty,
      });

      const valueText = new Text("", {
        font: new PhetFont(12),
        fill: QPPWColors.textFillProperty,
        minWidth: 40,
      });

      amplitudeProperty.link((value) => {
        valueText.string = value.toFixed(2);
      });

      const slider = new HSlider(amplitudeProperty, amplitudeProperty.range!, {
        trackSize: new Dimension2(200, 4),
        thumbSize: new Dimension2(15, 30),
      });

      const sliderRow = new HBox({
        spacing: 10,
        children: [stateLabel, slider, valueText],
      });

      sliderNodes.push(sliderRow);
    }

    // Update config when amplitudes change
    const updateConfig = () => {
      const amplitudes = amplitudeProperties.map((prop) => prop.value);
      const phases = config.phases.slice(0, numStates);

      // Extend phases array if needed
      while (phases.length < amplitudes.length) {
        phases.push(0);
      }

      configProperty.value = {
        ...config,
        amplitudes,
        phases,
      };
    };

    amplitudeProperties.forEach((prop) => {
      prop.link(updateConfig);
    });

    const slidersVBox = new VBox({
      spacing: 10,
      align: "left",
      children: sliderNodes,
    });

    // Normalization info
    const normalizationText = new Text("", {
      font: new PhetFont(12),
      fill: QPPWColors.textFillProperty,
    });

    const updateNormalization = () => {
      const sumSquared = amplitudeProperties.reduce((sum, prop) => sum + prop.value * prop.value, 0);
      normalizationText.string = `Sum of |cᵢ|² = ${sumSquared.toFixed(3)}`;

      // Change color if not normalized
      if (Math.abs(sumSquared - 1.0) > 0.01) {
        normalizationText.fill = "orange";
      } else {
        normalizationText.fill = QPPWColors.textFillProperty.value;
      }
    };

    amplitudeProperties.forEach((prop) => {
      prop.link(updateNormalization);
    });

    // Normalize button
    const normalizeButton = new RectangularPushButton({
      content: new Text("Normalize", {
        font: new PhetFont(12),
        fill: QPPWColors.textFillProperty,
      }),
      listener: () => {
        const sumSquared = amplitudeProperties.reduce((sum, prop) => sum + prop.value * prop.value, 0);
        const normFactor = Math.sqrt(sumSquared);

        if (normFactor > 0) {
          amplitudeProperties.forEach((prop) => {
            prop.value = prop.value / normFactor;
          });
        }
      },
      baseColor: QPPWColors.controlPanelBackgroundColorProperty,
    });

    const normalizationRow = new HBox({
      spacing: 10,
      children: [normalizationText, normalizeButton],
    });

    return new VBox({
      spacing: 15,
      align: "left",
      children: [titleText, descriptionText, slidersVBox, normalizationRow],
    });
  }
}
