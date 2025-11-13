/**
 * ControlPanelNode contains all the controls for the One Well screen.
 * This includes potential selection, display options, particle mass, and well parameters.
 */

import { Node, Text, VBox, HBox, HSeparator } from "scenerystack/scenery";
import { Panel, Checkbox, VerticalAquaRadioButtonGroup, ComboBox, HSlider } from "scenerystack/sun";
import { Dimension2 } from "scenerystack/dot";
import { OneWellModel } from "../../one-well/model/OneWellModel.js";
import { PotentialType } from "../model/PotentialFunction.js";
import QPPWColors from "../../QPPWColors.js";
import { ResetAllButton, PhetFont } from "scenerystack/scenery-phet";

interface ComboBoxItem<T> {
  value: T;
  createNode: () => Node;
}

export class ControlPanelNode extends Node {
  private readonly model: OneWellModel;
  private readonly resetCallback: () => void;

  public constructor(model: OneWellModel, resetCallback: () => void, listBoxParent: Node) {
    super();

    this.model = model;
    this.resetCallback = resetCallback;

    // Create all control groups
    const energyChartGroup = this.createEnergyChartGroup(listBoxParent);
    const bottomChartGroup = this.createBottomChartGroup();
    const particleMassGroup = this.createParticleMassGroup();
    const wellConfigGroup = this.createWellConfigurationGroup();
    const actionsGroup = this.createActionsGroup();

    // Arrange groups vertically
    const content = new VBox({
      spacing: 15,
      align: "left",
      children: [
        energyChartGroup,
        new HSeparator({ stroke: QPPWColors.gridLineProperty }),
        bottomChartGroup,
        new HSeparator({ stroke: QPPWColors.gridLineProperty }),
        particleMassGroup,
        new HSeparator({ stroke: QPPWColors.gridLineProperty }),
        wellConfigGroup,
        new HSeparator({ stroke: QPPWColors.gridLineProperty }),
        actionsGroup,
      ],
    });

    const panel = new Panel(content, {
      fill: QPPWColors.panelFillProperty,
      stroke: QPPWColors.panelStrokeProperty,
      xMargin: 15,
      yMargin: 15,
      cornerRadius: 5,
    });

    this.addChild(panel);
  }

  /**
   * Creates the Energy Chart control group.
   */
  private createEnergyChartGroup(listBoxParent: Node): Node {
    const titleText = new Text("Energy Chart", {
      font: new PhetFont({ size: 16, weight: "bold" }),
      fill: QPPWColors.textFillProperty,
    });

    // Potential Well dropdown
    const potentialItems: ComboBoxItem<PotentialType>[] = [
      {
        value: PotentialType.INFINITE_WELL,
        createNode: () => new Text("Square (Infinite)", { font: new PhetFont(14), fill: QPPWColors.textFillProperty }),
      },
      {
        value: PotentialType.FINITE_WELL,
        createNode: () => new Text("Square (Finite)", { font: new PhetFont(14), fill: QPPWColors.textFillProperty }),
      },
      {
        value: PotentialType.HARMONIC_OSCILLATOR,
        createNode: () => new Text("Harmonic Oscillator", { font: new PhetFont(14), fill: QPPWColors.textFillProperty }),
      },
      {
        value: PotentialType.ASYMMETRIC_TRIANGLE,
        createNode: () => new Text("Asymmetric Triangle", { font: new PhetFont(14), fill: QPPWColors.textFillProperty }),
      },
      {
        value: PotentialType.COULOMB_1D,
        createNode: () => new Text("1D Coulomb", { font: new PhetFont(14), fill: QPPWColors.textFillProperty }),
      },
      {
        value: PotentialType.COULOMB_3D,
        createNode: () => new Text("3D Coulomb", { font: new PhetFont(14), fill: QPPWColors.textFillProperty }),
      },
    ];

    const potentialComboBox = new ComboBox(this.model.potentialTypeProperty, potentialItems, listBoxParent, {
      xMargin: 8,
      yMargin: 6,
      cornerRadius: 4,
    });

    const potentialLabel = new Text("Potential Well:", {
      font: new PhetFont(14),
      fill: QPPWColors.textFillProperty,
    });

    const potentialRow = new HBox({
      spacing: 10,
      children: [potentialLabel, potentialComboBox],
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
        potentialRow,
        // configurePotentialButton,
        // superpositionButton,
      ],
    });
  }

  /**
   * Creates the Bottom Chart control group.
   */
  private createBottomChartGroup(): Node {
    const titleText = new Text("Bottom Chart", {
      font: new PhetFont({ size: 16, weight: "bold" }),
      fill: QPPWColors.textFillProperty,
    });

    // Display mode radio buttons
    const displayModeItems = [
      {
        value: "probabilityDensity" as const,
        createNode: () => new Text("Probability Density", { font: new PhetFont(14), fill: QPPWColors.textFillProperty }),
      },
      {
        value: "waveFunction" as const,
        createNode: () => new Text("Wave Function", { font: new PhetFont(14), fill: QPPWColors.textFillProperty }),
      },
    ];

    const displayModeGroup = new VerticalAquaRadioButtonGroup(this.model.displayModeProperty, displayModeItems, {
      spacing: 8,
      radioButtonOptions: {
        radius: 8,
      },
    });

    const displayLabel = new Text("Display:", {
      font: new PhetFont(14),
      fill: QPPWColors.textFillProperty,
    });

    // Wave function views checkboxes
    const waveFunctionViewsLabel = new Text("Wave Function views:", {
      font: new PhetFont(14),
      fill: QPPWColors.textFillProperty,
    });

    const realPartCheckbox = new Checkbox(
      this.model.showRealPartProperty,
      new Text("real part", { font: new PhetFont(12), fill: QPPWColors.textFillProperty }),
      { boxWidth: 16 },
    );

    const imaginaryPartCheckbox = new Checkbox(
      this.model.showImaginaryPartProperty,
      new Text("imaginary part", { font: new PhetFont(12), fill: QPPWColors.textFillProperty }),
      { boxWidth: 16 },
    );

    const magnitudeCheckbox = new Checkbox(
      this.model.showMagnitudeProperty,
      new Text("magnitude", { font: new PhetFont(12), fill: QPPWColors.textFillProperty }),
      { boxWidth: 16 },
    );

    const phaseCheckbox = new Checkbox(
      this.model.showPhaseProperty,
      new Text("phase", { font: new PhetFont(12), fill: QPPWColors.textFillProperty }),
      { boxWidth: 16 },
    );

    const waveFunctionCheckboxes = new VBox({
      spacing: 6,
      align: "left",
      children: [realPartCheckbox, imaginaryPartCheckbox, magnitudeCheckbox, phaseCheckbox],
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
      children: [titleText, displayLabel, displayModeGroup, waveFunctionViewsLabel, waveFunctionCheckboxes],
    });
  }

  /**
   * Creates the Particle Mass control group.
   */
  private createParticleMassGroup(): Node {
    const titleText = new Text("Particle Mass", {
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

    const massSlider = new HSlider(this.model.particleMassProperty, this.model.particleMassProperty.range!, {
      trackSize: new Dimension2(150, 4),
      thumbSize: new Dimension2(15, 30),
    });

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
    const titleText = new Text("Well Configuration", {
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

    const widthSlider = new HSlider(this.model.wellWidthProperty, this.model.wellWidthProperty.range!, {
      trackSize: new Dimension2(150, 4),
      thumbSize: new Dimension2(15, 30),
    });

    const widthRow = new VBox({
      spacing: 4,
      align: "left",
      children: [
        new Text("Well Width", { font: new PhetFont(12), fill: QPPWColors.textFillProperty }),
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

    const depthSlider = new HSlider(this.model.wellDepthProperty, this.model.wellDepthProperty.range!, {
      trackSize: new Dimension2(150, 4),
      thumbSize: new Dimension2(15, 30),
    });

    const depthRow = new VBox({
      spacing: 4,
      align: "left",
      children: [
        new Text("Well Depth", { font: new PhetFont(12), fill: QPPWColors.textFillProperty }),
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

  /**
   * Creates the Actions control group.
   */
  private createActionsGroup(): Node {
    const resetButton = new ResetAllButton({
      listener: this.resetCallback,
      radius: 20,
    });

    return new VBox({
      spacing: 8,
      align: "center",
      children: [resetButton],
    });
  }
}
