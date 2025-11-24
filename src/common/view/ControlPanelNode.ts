/**
 * ControlPanelNode contains all the controls for the One Well screen.
 * This includes potential selection, display options, particle mass, and well parameters.
 */

import {
  Node,
  Text,
  VBox,
  HBox,
  HSeparator,
  RichText,
} from "scenerystack/scenery";
import {
  Panel,
  Checkbox,
  VerticalAquaRadioButtonGroup,
  ComboBox,
  HSlider,
} from "scenerystack/sun";
import { Dimension2 } from "scenerystack/dot";
import { OneWellModel } from "../../one-well/model/OneWellModel.js";
import { TwoWellsModel } from "../../two-wells/model/TwoWellsModel.js";
import { ManyWellsModel } from "../../many-wells/model/ManyWellsModel.js";
import { PotentialType } from "../model/PotentialFunction.js";
import { SuperpositionType } from "../model/SuperpositionType.js";
import { SuperpositionDialog } from "./SuperpositionDialog.js";
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
  private readonly model: OneWellModel | TwoWellsModel | ManyWellsModel;
  private readonly options: ResolvedControlPanelNodeOptions;

  public constructor(
    model: OneWellModel | TwoWellsModel | ManyWellsModel,
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
    const particleMassGroup = this.options.showParticleMass
      ? this.createParticleMassGroup()
      : null;
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
      {
        value: PotentialType.DOUBLE_SQUARE_WELL,
        createNode: () =>
          new Text(stringManager.doubleSquareWellStringProperty, {
            font: new PhetFont(14),
            fill: QPPWColors.textFillProperty,
          }),
      },
      {
        value: PotentialType.MULTI_SQUARE_WELL,
        createNode: () =>
          new Text(stringManager.multiSquareWellStringProperty, {
            font: new PhetFont(14),
            fill: QPPWColors.textFillProperty,
          }),
      },
      {
        value: PotentialType.MULTI_COULOMB_1D,
        createNode: () =>
          new Text(stringManager.multiCoulomb1DStringProperty, {
            font: new PhetFont(14),
            fill: QPPWColors.textFillProperty,
          }),
      },
    ];

    // Filter potential types if specified in options
    const potentialItems = this.options.allowedPotentialTypes
      ? allPotentialItems.filter((item) =>
          this.options.allowedPotentialTypes!.includes(item.value),
        )
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

    // Superposition State dropdown
    const superpositionItems: ComboBoxItem<SuperpositionType>[] = [
      {
        value: SuperpositionType.PSI_I_PSI_J,
        createNode: () =>
          new RichText(stringManager.psiIPsiJStringProperty, {
            font: new PhetFont(14),
            fill: QPPWColors.textFillProperty,
          }),
      },
      {
        value: SuperpositionType.SINGLE,
        createNode: () =>
          new RichText(stringManager.psiKStringProperty, {
            font: new PhetFont(14),
            fill: QPPWColors.textFillProperty,
          }),
      },
      {
        value: SuperpositionType.LOCALIZED_NARROW,
        createNode: () =>
          new Text(stringManager.localizedNarrowStringProperty, {
            font: new PhetFont(14),
            fill: QPPWColors.textFillProperty,
          }),
      },
      {
        value: SuperpositionType.LOCALIZED_WIDE,
        createNode: () =>
          new Text(stringManager.localizedWideStringProperty, {
            font: new PhetFont(14),
            fill: QPPWColors.textFillProperty,
          }),
      },
      {
        value: SuperpositionType.COHERENT,
        createNode: () =>
          new Text(stringManager.coherentStateStringProperty, {
            font: new PhetFont(14),
            fill: QPPWColors.textFillProperty,
          }),
      },
      {
        value: SuperpositionType.CUSTOM,
        createNode: () =>
          new Text(stringManager.customStringProperty, {
            font: new PhetFont(14),
            fill: QPPWColors.textFillProperty,
          }),
      },
    ];

    const superpositionComboBox = new ComboBox(
      this.model.superpositionTypeProperty,
      superpositionItems,
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

    const superpositionLabelText = new Text(
      stringManager.superpositionStringProperty,
      {
        font: new PhetFont(14),
        fill: QPPWColors.textFillProperty,
      },
    );

    const superpositionRowNode = new HBox({
      spacing: 10,
      children: [superpositionLabelText, superpositionComboBox],
    });

    // Coherent state displacement slider (OneWellModel only)
    let displacementRow: Node | null = null;
    if ("coherentDisplacementProperty" in this.model) {
      const oneWellModel = this.model as OneWellModel;

      const displacementValueText = new Text("", {
        font: new PhetFont(12),
        fill: QPPWColors.textFillProperty,
      });

      oneWellModel.coherentDisplacementProperty.link((displacement: number) => {
        displacementValueText.string = `${displacement.toFixed(2)} nm`;
      });

      const displacementSlider = new HSlider(
        oneWellModel.coherentDisplacementProperty,
        oneWellModel.coherentDisplacementProperty.range!,
        {
          trackSize: new Dimension2(120, 4),
          thumbSize: new Dimension2(15, 30),
        },
      );

      displacementRow = new VBox({
        spacing: 4,
        align: "left",
        children: [
          new Text(stringManager.displacementStringProperty, {
            font: new PhetFont(12),
            fill: QPPWColors.textFillProperty,
          }),
          new HBox({
            spacing: 10,
            children: [displacementSlider, displacementValueText],
          }),
        ],
        visible: false, // Initially hidden
      });

      // Show/hide displacement slider based on superposition type
      this.model.superpositionTypeProperty.link((type: SuperpositionType) => {
        displacementRow!.visible = type === SuperpositionType.COHERENT;
      });
    }

    // Track the previous superposition type to revert if dialog is cancelled
    let previousSuperpositionType: SuperpositionType =
      this.model.superpositionTypeProperty.value;
    let isHandlingDialogResult = false;

    // Open dialog when "Custom..." is selected
    this.model.superpositionTypeProperty.link((type: SuperpositionType) => {
      // Skip if we're handling dialog result to avoid recursion
      if (isHandlingDialogResult) {
        return;
      }

      if (type === SuperpositionType.CUSTOM) {
        const dialog = new SuperpositionDialog(
          this.model.superpositionConfigProperty,
          this.model.getBoundStates(),
          () => {
            // OK button pressed - keep CUSTOM selection
            isHandlingDialogResult = true;
            previousSuperpositionType = SuperpositionType.CUSTOM;
            isHandlingDialogResult = false;
          },
          () => {
            // Cancel button pressed - revert to previous selection
            isHandlingDialogResult = true;
            this.model.superpositionTypeProperty.value =
              previousSuperpositionType;
            isHandlingDialogResult = false;
          },
        );
        dialog.show();
      } else {
        // Update previous type when user selects a non-CUSTOM option
        previousSuperpositionType = type;
      }
    });

    const children: Node[] = [
      titleText,
      potentialRowNode,
      superpositionRowNode,
    ];

    if (displacementRow) {
      children.push(displacementRow);
    }

    return new VBox({
      spacing: 8,
      align: "left",
      children: children,
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
    const waveFunctionViewsLabel = new Text(
      stringManager.waveFunctionViewsStringProperty,
      {
        font: new PhetFont(14),
        fill: QPPWColors.textFillProperty,
      },
    );

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
    this.model.displayModeProperty.link((mode: string) => {
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
    this.model.particleMassProperty.link((mass: number) => {
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

    // Well Depth slider (only for finite wells)
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

    let barrierHeightRow: Node | null = null;
    if ("barrierHeightProperty" in this.model) {
      const oneWellModel = this.model as OneWellModel;

      oneWellModel.barrierHeightProperty.link((height: number) => {
        barrierHeightValueText.string = `${height.toFixed(2)} eV`;
      });

      const barrierHeightSlider = new HSlider(
        oneWellModel.barrierHeightProperty,
        oneWellModel.barrierHeightProperty.range!,
        {
          trackSize: new Dimension2(150, 4),
          thumbSize: new Dimension2(15, 30),
        },
      );

      barrierHeightRow = new VBox({
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
    }

    // Potential Offset slider (only for triangular potential)
    const offsetValueText = new Text("", {
      font: new PhetFont(12),
      fill: QPPWColors.textFillProperty,
    });

    let offsetRow: Node | null = null;
    if ("potentialOffsetProperty" in this.model) {
      const oneWellModel = this.model as OneWellModel;

      oneWellModel.potentialOffsetProperty.link((offset: number) => {
        offsetValueText.string = `${offset.toFixed(2)} eV`;
      });

      const offsetSlider = new HSlider(
        oneWellModel.potentialOffsetProperty,
        oneWellModel.potentialOffsetProperty.range!,
        {
          trackSize: new Dimension2(150, 4),
          thumbSize: new Dimension2(15, 30),
        },
      );

      offsetRow = new VBox({
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
    }

    // Well Separation slider (only for double/multi square well and multi-Coulomb 1D)
    const separationValueText = new Text("", {
      font: new PhetFont(12),
      fill: QPPWColors.textFillProperty,
    });

    // Check if the model has wellSeparationProperty (TwoWellsModel or ManyWellsModel)
    let separationRow: Node | null = null;
    if ("wellSeparationProperty" in this.model) {
      const wellModel = this.model as TwoWellsModel | ManyWellsModel;

      wellModel.wellSeparationProperty.link((separation: number) => {
        separationValueText.string = `${separation.toFixed(2)} nm`;
      });

      const separationSlider = new HSlider(
        wellModel.wellSeparationProperty,
        wellModel.wellSeparationProperty.range!,
        {
          trackSize: new Dimension2(150, 4),
          thumbSize: new Dimension2(15, 30),
        },
      );

      separationRow = new VBox({
        spacing: 4,
        align: "left",
        children: [
          new Text(stringManager.wellSeparationStringProperty, {
            font: new PhetFont(12),
            fill: QPPWColors.textFillProperty,
          }),
          new HBox({
            spacing: 10,
            children: [separationSlider, separationValueText],
          }),
        ],
      });
    }

    // Number of Wells slider (only for multi-well potentials)
    const numberOfWellsValueText = new Text("", {
      font: new PhetFont(12),
      fill: QPPWColors.textFillProperty,
    });

    // Check if the model has numberOfWellsProperty (ManyWellsModel)
    let numberOfWellsRow: Node | null = null;
    if ("numberOfWellsProperty" in this.model) {
      const manyWellsModel = this.model as ManyWellsModel;

      manyWellsModel.numberOfWellsProperty.link((numberOfWells: number) => {
        numberOfWellsValueText.string = `${numberOfWells}`;
      });

      const numberOfWellsSlider = new HSlider(
        manyWellsModel.numberOfWellsProperty,
        manyWellsModel.numberOfWellsProperty.range!,
        {
          trackSize: new Dimension2(150, 4),
          thumbSize: new Dimension2(15, 30),
          constrainValue: (value: number) => Math.round(value), // Integer values only
        },
      );

      numberOfWellsRow = new VBox({
        spacing: 4,
        align: "left",
        children: [
          new Text(stringManager.numberOfWellsStringProperty, {
            font: new PhetFont(12),
            fill: QPPWColors.textFillProperty,
          }),
          new HBox({
            spacing: 10,
            children: [numberOfWellsSlider, numberOfWellsValueText],
          }),
        ],
      });
    }

    // Enable/disable width slider based on potential type
    // Coulomb potentials have fixed spatial extent and don't use well width
    this.model.potentialTypeProperty.link((type: PotentialType) => {
      const needsWidth =
        type !== PotentialType.COULOMB_1D && type !== PotentialType.COULOMB_3D;
      widthRow.visible = needsWidth;
    });

    // Enable/disable depth slider based on potential type
    this.model.potentialTypeProperty.link((type: PotentialType) => {
      const needsDepth =
        type === PotentialType.FINITE_WELL ||
        type === PotentialType.HARMONIC_OSCILLATOR ||
        type === PotentialType.MORSE ||
        type === PotentialType.POSCHL_TELLER ||
        type === PotentialType.ROSEN_MORSE ||
        type === PotentialType.ECKART ||
        type === PotentialType.ASYMMETRIC_TRIANGLE ||
        type === PotentialType.TRIANGULAR ||
        type === PotentialType.DOUBLE_SQUARE_WELL ||
        type === PotentialType.MULTI_SQUARE_WELL;
      depthRow.visible = needsDepth;
    });

    // Enable/disable barrier height slider based on potential type (only for Rosen-Morse and Eckart)
    if (barrierHeightRow) {
      this.model.potentialTypeProperty.link((type: PotentialType) => {
        const needsBarrierHeight =
          type === PotentialType.ROSEN_MORSE || type === PotentialType.ECKART;
        barrierHeightRow!.visible = needsBarrierHeight;
      });
    }

    // Enable/disable offset slider based on potential type (only for triangular potential)
    if (offsetRow) {
      this.model.potentialTypeProperty.link((type: PotentialType) => {
        const needsOffset = type === PotentialType.TRIANGULAR;
        offsetRow!.visible = needsOffset;
      });
    }

    // Enable/disable separation slider based on potential type
    if (separationRow) {
      this.model.potentialTypeProperty.link((type: PotentialType) => {
        const needsSeparation =
          type === PotentialType.DOUBLE_SQUARE_WELL ||
          type === PotentialType.MULTI_SQUARE_WELL ||
          type === PotentialType.MULTI_COULOMB_1D;
        separationRow!.visible = needsSeparation;
      });
    }

    // Enable/disable number of wells slider based on potential type
    if (numberOfWellsRow) {
      this.model.potentialTypeProperty.link((type: PotentialType) => {
        const needsNumberOfWells =
          type === PotentialType.MULTI_SQUARE_WELL ||
          type === PotentialType.MULTI_COULOMB_1D;
        numberOfWellsRow!.visible = needsNumberOfWells;
      });
    }

    // Electric Field slider (only for ManyWellsModel)
    const electricFieldValueText = new Text("", {
      font: new PhetFont(12),
      fill: QPPWColors.textFillProperty,
    });

    let electricFieldRow: Node | null = null;
    if ("electricFieldProperty" in this.model) {
      const manyWellsModel = this.model as ManyWellsModel;

      manyWellsModel.electricFieldProperty.link((field: number) => {
        electricFieldValueText.string = `${field.toFixed(3)} eV/nm`;
      });

      const electricFieldSlider = new HSlider(
        manyWellsModel.electricFieldProperty,
        manyWellsModel.electricFieldProperty.range!,
        {
          trackSize: new Dimension2(150, 4),
          thumbSize: new Dimension2(15, 30),
        },
      );

      electricFieldRow = new VBox({
        spacing: 4,
        align: "left",
        children: [
          new Text(stringManager.electricFieldStringProperty, {
            font: new PhetFont(12),
            fill: QPPWColors.textFillProperty,
          }),
          new HBox({
            spacing: 10,
            children: [electricFieldSlider, electricFieldValueText],
          }),
        ],
      });
    }

    // Enable/disable electric field slider based on potential type
    if (electricFieldRow) {
      this.model.potentialTypeProperty.link((type: PotentialType) => {
        const needsElectricField =
          type === PotentialType.MULTI_SQUARE_WELL ||
          type === PotentialType.MULTI_COULOMB_1D;
        electricFieldRow!.visible = needsElectricField;
      });
    }

    const children: Node[] = [titleText, widthRow, depthRow];
    if (barrierHeightRow) {
      children.push(barrierHeightRow);
    }
    if (offsetRow) {
      children.push(offsetRow);
    }
    if (numberOfWellsRow) {
      children.push(numberOfWellsRow);
    }
    if (separationRow) {
      children.push(separationRow);
    }
    if (electricFieldRow) {
      children.push(electricFieldRow);
    }

    return new VBox({
      spacing: 10,
      align: "left",
      children: children,
    });
  }
}
