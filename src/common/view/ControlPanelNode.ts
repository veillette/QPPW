/**
 * ControlPanelNode contains all the controls for the One Well screen.
 * This includes potential selection, display options, particle mass, and well parameters.
 */

import { Node, Text, VBox, HBox, HSeparator } from "scenerystack/scenery";
import { Panel, Checkbox, VerticalAquaRadioButtonGroup, ComboBox, HSlider } from "scenerystack/sun";
import { Dimension2 } from "scenerystack/dot";
import { OneWellModel } from "../../one-well/model/OneWellModel.js";
import { TwoWellsModel } from "../../two-wells/model/TwoWellsModel.js";
import { PotentialType } from "../model/PotentialFunction.js";
import QPPWColors from "../../QPPWColors.js";
import { PhetFont } from "scenerystack/scenery-phet";
import stringManager from "../../i18n/StringManager.js";

interface ComboBoxItem<T> {
  value: T;
  createNode: () => Node;
}

export interface ControlPanelNodeOptions {
  // Whether to show the particle mass control group
  showParticleMass?: boolean;
  // Filter which potential types to show (if undefined, shows all)
  allowedPotentialTypes?: PotentialType[];
}

interface ResolvedControlPanelNodeOptions {
  showParticleMass: boolean;
  allowedPotentialTypes: PotentialType[] | undefined;
}

export class ControlPanelNode extends Node {
  private readonly model: OneWellModel | TwoWellsModel;
  private readonly options: ResolvedControlPanelNodeOptions;

  public constructor(
    model: OneWellModel | TwoWellsModel,
    listBoxParent: Node,
    providedOptions?: ControlPanelNodeOptions,
  ) {
    super();

    this.model = model;

    // Default options
    this.options = {
      showParticleMass: true,
      allowedPotentialTypes: undefined,
      ...providedOptions,
    };

    // Create all control groups
    const energyChartGroup = this.createEnergyChartGroup(listBoxParent);
    const bottomChartGroup = this.createBottomChartGroup();
    const particleMassGroup = this.options.showParticleMass ? this.createParticleMassGroup() : null;
    const wellConfigGroup = this.createWellConfigurationGroup();

    // Arrange groups vertically (only include particle mass if enabled)
    const children: Node[] = [
      energyChartGroup,
      new HSeparator({ stroke: QPPWColors.gridLineProperty }),
      bottomChartGroup,
      new HSeparator({ stroke: QPPWColors.gridLineProperty }),
    ];

    if (particleMassGroup) {
      children.push(particleMassGroup);
      children.push(new HSeparator({ stroke: QPPWColors.gridLineProperty }));
    }

    children.push(wellConfigGroup);
    children.push(new HSeparator({ stroke: QPPWColors.gridLineProperty }));

    // Arrange groups vertically
    const content = new VBox({
      spacing: 15,
      align: "left",
      children: children,
    });

    const controlPanel = new Panel(content, {
      fill: QPPWColors.panelFillProperty,
      stroke: QPPWColors.panelStrokeProperty,
      xMargin: 15,
      yMargin: 15,
      cornerRadius: 5,
    });

    this.addChild(controlPanel);
  }

  /**
   * Creates the Energy Chart control group.
   */
  private createEnergyChartGroup(listBoxParent: Node): Node {
    const titleText = new Text(stringManager.energyChartStringProperty, {
      font: new PhetFont({ size: 16, weight: "bold" }),
      fill: QPPWColors.textFillProperty,
    });

    // Potential Well dropdown - all available options
    const allPotentialItems: ComboBoxItem<PotentialType>[] = [
      {
        value: PotentialType.INFINITE_WELL,
        createNode: () =>
          new Text(stringManager.squareInfiniteStringProperty, {
            font: new PhetFont(14),
            fill: QPPWColors.textFillProperty,
          }),
      },
      {
        value: PotentialType.FINITE_WELL,
        createNode: () =>
          new Text(stringManager.squareFiniteStringProperty, {
            font: new PhetFont(14),
            fill: QPPWColors.textFillProperty,
          }),
      },
      {
        value: PotentialType.HARMONIC_OSCILLATOR,
        createNode: () =>
          new Text(stringManager.harmonicOscillatorStringProperty, {
            font: new PhetFont(14),
            fill: QPPWColors.textFillProperty,
          }),
      },
      {
        value: PotentialType.ASYMMETRIC_TRIANGLE,
        createNode: () =>
          new Text(stringManager.asymmetricTriangleStringProperty, {
            font: new PhetFont(14),
            fill: QPPWColors.textFillProperty,
          }),
      },
      {
        value: PotentialType.COULOMB_1D,
        createNode: () =>
          new Text(stringManager.coulomb1DStringProperty, {
            font: new PhetFont(14),
            fill: QPPWColors.textFillProperty,
          }),
      },
      {
        value: PotentialType.COULOMB_3D,
        createNode: () =>
          new Text(stringManager.coulomb3DStringProperty, {
            font: new PhetFont(14),
            fill: QPPWColors.textFillProperty,
          }),
      },
    ];

    // Filter potential types if specified in options
    const potentialItems = this.options.allowedPotentialTypes
      ? allPotentialItems.filter((item) => this.options.allowedPotentialTypes!.includes(item.value))
      : allPotentialItems;

    const potentialComboBox = new ComboBox(
      this.model.potentialTypeProperty,
      potentialItems,
      listBoxParent,
      {
        xMargin: 8,
        yMargin: 6,
        cornerRadius: 4,
        buttonFill: QPPWColors.controlPanelBackgroundColorProperty,
        buttonStroke: QPPWColors.controlPanelStrokeColorProperty,
        listFill: QPPWColors.controlPanelBackgroundColorProperty,
        listStroke: QPPWColors.controlPanelStrokeColorProperty,
        highlightFill: QPPWColors.controlPanelStrokeColorProperty,
      },
    );

    const potentialLabelText = new Text(stringManager.potentialWellStringProperty, {
      font: new PhetFont(14),
      fill: QPPWColors.textFillProperty,
    });

    const potentialRowNode = new HBox({
      spacing: 10,
      children: [potentialLabelText, potentialComboBox],
    });

    // Note: Configure Potential and Superposition State buttons would open dialogs
    // For now, we'll add placeholders
    // const configurePotentialButton = new RectangularPushButton({
    //   content: new Text("Configure Potential...", { font: new PhetFont(12) }),
    //   listener: () => { /* Open dialog */ },
    // });

    // const superpositionButton = new RectangularPushButton({
    //   content: new Text("Superposition State...", { font: new PhetFont(12) }),
    //   listener: () => { /* Open dialog */ },
    // });

    return new VBox({
      spacing: 8,
      align: "left",
      children: [
        titleText,
        potentialRowNode,
        // configurePotentialButton,
        // superpositionButton,
      ],
    });
  }

  /**
   * Creates the Bottom Chart control group.
   */
  private createBottomChartGroup(): Node {
    // Display mode radio buttons
    const displayModeItems = [
      {
        value: "probabilityDensity" as const,
        createNode: () =>
          new Text(stringManager.probabilityDensityStringProperty, {
            font: new PhetFont(14),
            fill: QPPWColors.textFillProperty,
          }),
      },
      {
        value: "waveFunction" as const,
        createNode: () =>
          new Text(stringManager.wavefunctionStringProperty, {
            font: new PhetFont(14),
            fill: QPPWColors.textFillProperty,
          }),
      },
      {
        value: "phaseColor" as const,
        createNode: () =>
          new Text(stringManager.phaseColorStringProperty, {
            font: new PhetFont(14),
            fill: QPPWColors.textFillProperty,
          }),
      },
    ];

    const displayModeRadioButtonGroup = new VerticalAquaRadioButtonGroup(
      this.model.displayModeProperty,
      displayModeItems,
      {
        spacing: 8,
        radioButtonOptions: {
          radius: 8,
        },
      },
    );

    const displayLabel = new Text(stringManager.displayStringProperty, {
      font: new PhetFont(14),
      fill: QPPWColors.textFillProperty,
    });

    // Wave function views checkboxes
    const waveFunctionViewsLabel = new Text(stringManager.waveFunctionViewsStringProperty, {
      font: new PhetFont(14),
      fill: QPPWColors.textFillProperty,
    });

    const realPartCheckbox = new Checkbox(
      this.model.showRealPartProperty,
      new Text(stringManager.realPartStringProperty, {
        font: new PhetFont(12),
        fill: QPPWColors.textFillProperty,
      }),
      { boxWidth: 16 },
    );

    const imaginaryPartCheckbox = new Checkbox(
      this.model.showImaginaryPartProperty,
      new Text(stringManager.imaginaryPartStringProperty, {
        font: new PhetFont(12),
        fill: QPPWColors.textFillProperty,
      }),
      { boxWidth: 16 },
    );

    const magnitudeCheckbox = new Checkbox(
      this.model.showMagnitudeProperty,
      new Text(stringManager.magnitudeStringProperty, {
        font: new PhetFont(12),
        fill: QPPWColors.textFillProperty,
      }),
      { boxWidth: 16 },
    );

    const phaseCheckbox = new Checkbox(
      this.model.showPhaseProperty,
      new Text(stringManager.phaseStringProperty, {
        font: new PhetFont(12),
        fill: QPPWColors.textFillProperty,
      }),
      { boxWidth: 16 },
    );

    const waveFunctionCheckboxes = new VBox({
      spacing: 6,
      align: "left",
      children: [
        realPartCheckbox,
        imaginaryPartCheckbox,
        magnitudeCheckbox,
        phaseCheckbox,
      ],
      leftMargin: 20,
    });

    // Enable/disable wave function views based on display mode
    this.model.displayModeProperty.link((mode) => {
      const enabled = mode === "waveFunction";
      realPartCheckbox.enabled = enabled;
      imaginaryPartCheckbox.enabled = enabled;
      magnitudeCheckbox.enabled = enabled;
      phaseCheckbox.enabled = enabled;
    });

    return new VBox({
      spacing: 8,
      align: "left",
      children: [
        displayLabel,
        displayModeRadioButtonGroup,
        waveFunctionViewsLabel,
        waveFunctionCheckboxes,
      ],
    });
  }

  /**
   * Creates the Particle Mass control group.
   */
  private createParticleMassGroup(): Node {
    const titleText = new Text(stringManager.particleMassStringProperty, {
      font: new PhetFont({ size: 16, weight: "bold" }),
      fill: QPPWColors.textFillProperty,
    });

    const massValueText = new Text("", {
      font: new PhetFont(14),
      fill: QPPWColors.textFillProperty,
    });

    // Update text when mass changes
    this.model.particleMassProperty.link((mass) => {
      massValueText.string = `${mass.toFixed(2)} mₑ`;
    });

    const massSlider = new HSlider(
      this.model.particleMassProperty,
      this.model.particleMassProperty.range!,
      {
        trackSize: new Dimension2(150, 4),
        thumbSize: new Dimension2(15, 30),
      },
    );

    return new VBox({
      spacing: 8,
      align: "left",
      children: [
        titleText,
        new HBox({
          spacing: 10,
          children: [massSlider, massValueText],
        }),
      ],
    });
  }

  /**
   * Creates the Well Configuration control group.
   */
  private createWellConfigurationGroup(): Node {
    const titleText = new Text(stringManager.wellConfigurationStringProperty, {
      font: new PhetFont({ size: 16, weight: "bold" }),
      fill: QPPWColors.textFillProperty,
    });

    // Well Width slider
    const widthValueText = new Text("", {
      font: new PhetFont(12),
      fill: QPPWColors.textFillProperty,
    });

    this.model.wellWidthProperty.link((width) => {
      widthValueText.string = `${width.toFixed(2)} nm`;
    });

    const widthSlider = new HSlider(
      this.model.wellWidthProperty,
      this.model.wellWidthProperty.range!,
      {
        trackSize: new Dimension2(150, 4),
        thumbSize: new Dimension2(15, 30),
      },
    );

    const widthRow = new VBox({
      spacing: 4,
      align: "left",
      children: [
        new Text(stringManager.wellWidthStringProperty, {
          font: new PhetFont(12),
          fill: QPPWColors.textFillProperty,
        }),
        new HBox({
          spacing: 10,
          children: [widthSlider, widthValueText],
        }),
      ],
    });

    // Well Depth slider (only for finite wells)
    const depthValueText = new Text("", {
      font: new PhetFont(12),
      fill: QPPWColors.textFillProperty,
    });

    this.model.wellDepthProperty.link((depth) => {
      depthValueText.string = `${depth.toFixed(2)} eV`;
    });

    const depthSlider = new HSlider(
      this.model.wellDepthProperty,
      this.model.wellDepthProperty.range!,
      {
        trackSize: new Dimension2(150, 4),
        thumbSize: new Dimension2(15, 30),
      },
    );

    const depthRow = new VBox({
      spacing: 4,
      align: "left",
      children: [
        new Text(stringManager.wellDepthStringProperty, {
          font: new PhetFont(12),
          fill: QPPWColors.textFillProperty,
        }),
        new HBox({
          spacing: 10,
          children: [depthSlider, depthValueText],
        }),
      ],
    });

    // Enable/disable depth slider based on potential type
    this.model.potentialTypeProperty.link((type) => {
      const needsDepth =
        type === PotentialType.FINITE_WELL ||
        type === PotentialType.HARMONIC_OSCILLATOR ||
        type === PotentialType.ASYMMETRIC_TRIANGLE;
      depthRow.visible = needsDepth;
    });

    return new VBox({
      spacing: 10,
      align: "left",
      children: [titleText, widthRow, depthRow],
    });
  }
}