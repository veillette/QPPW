/**
 * IntroControlPanelNode is a simplified control panel for the intro screen.
 * It excludes the superposition combo box and phase color display mode.
 */

import { Node, Text, VBox, HBox, HSeparator } from "scenerystack/scenery";
import { Panel, Checkbox, ComboBox, HSlider } from "scenerystack/sun";
import { Dimension2 } from "scenerystack/dot";
import { IntroModel } from "../model/IntroModel.js";
import { PotentialType } from "../../common/model/PotentialFunction.js";
import QPPWColors from "../../QPPWColors.js";
import { PhetFont } from "scenerystack/scenery-phet";
import stringManager from "../../i18n/StringManager.js";

interface ComboBoxItem<T> {
  value: T;
  createNode: () => Node;
}

export class IntroControlPanelNode extends Node {
  private readonly model: IntroModel;
  private readonly probabilityChartNode?: import("../../common/view/WaveFunctionChartNode.js").WaveFunctionChartNode;
  private readonly waveFunctionChartNode?: import("../../common/view/WaveFunctionChartNode.js").WaveFunctionChartNode;

  public constructor(
    model: IntroModel,
    listBoxParent: Node,
    probabilityChartNode?: import("../../common/view/WaveFunctionChartNode.js").WaveFunctionChartNode,
    waveFunctionChartNode?: import("../../common/view/WaveFunctionChartNode.js").WaveFunctionChartNode,
  ) {
    super();

    this.model = model;
    this.probabilityChartNode = probabilityChartNode;
    this.waveFunctionChartNode = waveFunctionChartNode;

    // Create all control groups
    const energyChartGroup = this.createEnergyChartGroup(listBoxParent);
    const bottomChartGroup = this.createBottomChartGroup();
    const wellConfigGroup = this.createWellConfigurationGroup();

    // Arrange groups vertically
    const children: Node[] = [
      energyChartGroup,
      new HSeparator({ stroke: QPPWColors.gridLineProperty }),
      bottomChartGroup,
      new HSeparator({ stroke: QPPWColors.gridLineProperty }),
      wellConfigGroup,
      new HSeparator({ stroke: QPPWColors.gridLineProperty }),
    ];

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
   * Creates the Energy Chart control group (without superposition).
   */
  private createEnergyChartGroup(listBoxParent: Node): Node {
    const titleText = new Text(stringManager.energyChartStringProperty, {
      font: new PhetFont({ size: 16, weight: "bold" }),
      fill: QPPWColors.textFillProperty,
    });

    // Potential Well dropdown - limited to intro-friendly potentials
    const potentialItems: ComboBoxItem<PotentialType>[] = [
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
        value: PotentialType.MORSE,
        createNode: () =>
          new Text(stringManager.morseStringProperty, {
            font: new PhetFont(14),
            fill: QPPWColors.textFillProperty,
          }),
      },
      {
        value: PotentialType.POSCHL_TELLER,
        createNode: () =>
          new Text(stringManager.poschlTellerStringProperty, {
            font: new PhetFont(14),
            fill: QPPWColors.textFillProperty,
          }),
      },
      {
        value: PotentialType.ROSEN_MORSE,
        createNode: () =>
          new Text(stringManager.rosenMorseStringProperty, {
            font: new PhetFont(14),
            fill: QPPWColors.textFillProperty,
          }),
      },
      {
        value: PotentialType.ECKART,
        createNode: () =>
          new Text(stringManager.eckartStringProperty, {
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
        value: PotentialType.TRIANGULAR,
        createNode: () =>
          new Text(stringManager.triangularStringProperty, {
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

    const potentialLabelText = new Text(
      stringManager.potentialWellStringProperty,
      {
        font: new PhetFont(14),
        fill: QPPWColors.textFillProperty,
      },
    );

    const potentialRowNode = new HBox({
      spacing: 10,
      children: [potentialLabelText, potentialComboBox],
    });

    // No superposition dropdown in intro screen

    const children: Node[] = [titleText, potentialRowNode];

    return new VBox({
      spacing: 8,
      align: "left",
      children: children,
    });
  }

  /**
   * Creates the Bottom Chart control group (for probability density chart only).
   * Display mode controls removed since intro screen shows both charts separately.
   */
  private createBottomChartGroup(): Node {
    // Classical probability checkbox
    const classicalProbabilityCheckboxContent = new Checkbox(
      this.model.showClassicalProbabilityProperty,
      new Text(stringManager.classicalProbabilityDensityStringProperty, {
        font: new PhetFont(12),
        fill: QPPWColors.textFillProperty,
      }),
      { boxWidth: 16 },
    );

    const classicalProbabilityCheckbox = new Node({
      children: [classicalProbabilityCheckboxContent],
      x: 20,
    });

    // Show Zeros checkbox
    const showZerosCheckboxContent = new Checkbox(
      this.model.showZerosProperty,
      new Text(stringManager.showZerosStringProperty, {
        font: new PhetFont(12),
        fill: QPPWColors.textFillProperty,
      }),
      { boxWidth: 16 },
    );

    const showZerosCheckbox = new Node({
      children: [showZerosCheckboxContent],
      x: 20,
    });

    // Area measurement tool checkbox (for probability density chart)
    const areaToolCheckboxContent = this.probabilityChartNode
      ? new Checkbox(
          this.probabilityChartNode.showAreaToolProperty,
          new Text("Measure Area", {
            font: new PhetFont(12),
            fill: QPPWColors.textFillProperty,
          }),
          { boxWidth: 16 },
        )
      : null;

    const areaToolCheckbox = areaToolCheckboxContent
      ? new Node({
          children: [areaToolCheckboxContent],
          x: 20,
        })
      : null;

    // Curvature visualization tool checkbox (for wavefunction chart)
    const curvatureToolCheckboxContent = this.waveFunctionChartNode
      ? new Checkbox(
          this.waveFunctionChartNode.showCurvatureToolProperty,
          new Text("Show Curvature", {
            font: new PhetFont(12),
            fill: QPPWColors.textFillProperty,
          }),
          { boxWidth: 16 },
        )
      : null;

    const curvatureToolCheckbox = curvatureToolCheckboxContent
      ? new Node({
          children: [curvatureToolCheckboxContent],
          x: 20,
        })
      : null;

    // Build children array (no display mode controls since intro screen shows separate charts)
    const children: Node[] = [
      classicalProbabilityCheckbox,
      showZerosCheckbox,
    ];

    // Add area tool checkbox if it exists
    if (areaToolCheckbox) {
      children.push(areaToolCheckbox);
    }

    // Add curvature tool checkbox if it exists
    if (curvatureToolCheckbox) {
      children.push(curvatureToolCheckbox);
    }

    return new VBox({
      spacing: 8,
      align: "left",
      children,
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

    this.model.wellWidthProperty.link((width: number) => {
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

    // Well Depth slider
    const depthValueText = new Text("", {
      font: new PhetFont(12),
      fill: QPPWColors.textFillProperty,
    });

    this.model.wellDepthProperty.link((depth: number) => {
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

    // Barrier Height slider (only for Rosen-Morse and Eckart potentials)
    const barrierHeightValueText = new Text("", {
      font: new PhetFont(12),
      fill: QPPWColors.textFillProperty,
    });

    this.model.barrierHeightProperty.link((height: number) => {
      barrierHeightValueText.string = `${height.toFixed(2)} eV`;
    });

    const barrierHeightSlider = new HSlider(
      this.model.barrierHeightProperty,
      this.model.barrierHeightProperty.range!,
      {
        trackSize: new Dimension2(150, 4),
        thumbSize: new Dimension2(15, 30),
      },
    );

    const barrierHeightRow = new VBox({
      spacing: 4,
      align: "left",
      children: [
        new Text(stringManager.barrierHeightStringProperty, {
          font: new PhetFont(12),
          fill: QPPWColors.textFillProperty,
        }),
        new HBox({
          spacing: 10,
          children: [barrierHeightSlider, barrierHeightValueText],
        }),
      ],
    });

    // Potential Offset slider (only for triangular potential)
    const offsetValueText = new Text("", {
      font: new PhetFont(12),
      fill: QPPWColors.textFillProperty,
    });

    this.model.potentialOffsetProperty.link((offset: number) => {
      offsetValueText.string = `${offset.toFixed(2)} eV`;
    });

    const offsetSlider = new HSlider(
      this.model.potentialOffsetProperty,
      this.model.potentialOffsetProperty.range!,
      {
        trackSize: new Dimension2(150, 4),
        thumbSize: new Dimension2(15, 30),
      },
    );

    const offsetRow = new VBox({
      spacing: 4,
      align: "left",
      children: [
        new Text(stringManager.potentialOffsetStringProperty, {
          font: new PhetFont(12),
          fill: QPPWColors.textFillProperty,
        }),
        new HBox({
          spacing: 10,
          children: [offsetSlider, offsetValueText],
        }),
      ],
    });

    // Show/hide sliders based on potential type
    this.model.potentialTypeProperty.link((type: PotentialType) => {
      // Depth slider visibility
      const showDepth =
        type === PotentialType.FINITE_WELL ||
        type === PotentialType.HARMONIC_OSCILLATOR ||
        type === PotentialType.MORSE ||
        type === PotentialType.POSCHL_TELLER ||
        type === PotentialType.ROSEN_MORSE ||
        type === PotentialType.ECKART ||
        type === PotentialType.ASYMMETRIC_TRIANGLE ||
        type === PotentialType.TRIANGULAR ||
        type === PotentialType.COULOMB_1D ||
        type === PotentialType.COULOMB_3D;
      depthRow.visible = showDepth;

      // Barrier height slider visibility (only for Rosen-Morse and Eckart)
      const showBarrierHeight =
        type === PotentialType.ROSEN_MORSE || type === PotentialType.ECKART;
      barrierHeightRow.visible = showBarrierHeight;

      // Offset slider visibility (only for triangular potential)
      const showOffset = type === PotentialType.TRIANGULAR;
      offsetRow.visible = showOffset;
    });

    return new VBox({
      spacing: 8,
      align: "left",
      children: [titleText, widthRow, depthRow, barrierHeightRow, offsetRow],
    });
  }
}
