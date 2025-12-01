// NOTE: brand.js needs to be the first import. This is because SceneryStack for sims needs a very specific loading
// order: init.ts => assert.ts => splash.ts => brand.ts => everything else (here)
import "./brand.js";

import { onReadyToLaunch, Sim, PreferencesModel } from "scenerystack/sim";
import { Tandem } from "scenerystack/tandem";
import { VBox, Text, HStrut, HBox } from "scenerystack/scenery";
import {
  Checkbox,
  VerticalAquaRadioButtonGroup,
  HSlider,
} from "scenerystack/sun";
import { PhetFont } from "scenerystack/scenery-phet";
import { Dimension2, Range } from "scenerystack";
import { NumberProperty } from "scenerystack/axon";
import stringManager from "./i18n/StringManager.js";
import QPPWColors from "./QPPWColors.js";
import QPPWPreferences from "./QPPWPreferences.js";
import { NumericalMethod } from "./common/model/Schrodinger1DSolver.js";
import { OneWellScreen } from "./one-well/OneWellScreen.js";
import { TwoWellsScreen } from "./two-wells/TwoWellsScreen.js";
import { ManyWellsScreen } from "./many-wells/ManyWellsScreen.js";
import { IntroScreen } from "./intro/IntroScreen.js";
import { OneWellScreenIcon } from "./one-well/view/OneWellScreenIcon.js";
import { TwoWellsScreenIcon } from "./two-wells/view/TwoWellsScreenIcon.js";
import { ManyWellsScreenIcon } from "./many-wells/view/ManyWellsScreenIcon.js";
import { IntroScreenIcon } from "./intro/view/IntroScreenIcon.js";

onReadyToLaunch(() => {
  const screenNames = stringManager.getScreenNames();

  // Get preferences strings
  const preferencesLabels = stringManager.getPreferencesLabels();
  const numericalMethodNames = stringManager.getNumericalMethodNames();
  const numericalMethodDescriptions =
    stringManager.getNumericalMethodDescriptions();

  const simOptions = {
    webgl: true,
    preferencesModel: new PreferencesModel({
      visualOptions: {
        supportsProjectorMode: true,
        supportsInteractiveHighlights: true,
      },
      audioOptions: {
        supportsVoicing: true,
        supportsSound: true,
      },
      simulationOptions: {
        customPreferences: [
          {
            createContent: (_tandem: Tandem) => {
              // Auto-pause preference
              const autoPauseSectionVBox = new VBox({
                align: "left",
                spacing: 8,
                children: [
                  new Checkbox(
                    QPPWPreferences.autoPauseWhenTabHiddenProperty,
                    new Text(
                      preferencesLabels.autoPauseWhenTabHiddenStringProperty,
                      {
                        font: new PhetFont(16),
                      },
                    ),
                    {
                      boxWidth: 16,
                    },
                  ),
                  new Text(
                    preferencesLabels.autoPauseDescriptionStringProperty,
                    {
                      font: new PhetFont(12),
                      maxWidth: 600,
                    },
                  ),
                ],
              });

              // Numerical method section
              const numericalMethodRadioButtonGroup =
                new VerticalAquaRadioButtonGroup(
                  QPPWPreferences.numericalMethodProperty,
                  [
                    {
                      value: NumericalMethod.NUMEROV,
                      createNode: () =>
                        new VBox({
                          align: "left",
                          spacing: 4,
                          children: [
                            new Text(
                              numericalMethodNames.numerovStringProperty,
                              {
                                font: new PhetFont(14),
                              },
                            ),
                            new Text(
                              numericalMethodDescriptions.numerovStringProperty,
                              {
                                font: new PhetFont(11),
                                fill: "rgb(80,80,80)",
                                maxWidth: 550,
                              },
                            ),
                          ],
                        }),
                      tandemName: "numerovRadioButton",
                    },
                    {
                      value: NumericalMethod.MATRIX_NUMEROV,
                      createNode: () =>
                        new VBox({
                          align: "left",
                          spacing: 4,
                          children: [
                            new Text(
                              numericalMethodNames.matrixNumerovStringProperty,
                              {
                                font: new PhetFont(14),
                              },
                            ),
                            new Text(
                              numericalMethodDescriptions.matrixNumerovStringProperty,
                              {
                                font: new PhetFont(11),
                                fill: "rgb(80,80,80)",
                                maxWidth: 550,
                              },
                            ),
                          ],
                        }),
                      tandemName: "matrixNumerovRadioButton",
                    },
                    {
                      value: NumericalMethod.DVR,
                      createNode: () =>
                        new VBox({
                          align: "left",
                          spacing: 4,
                          children: [
                            new Text(numericalMethodNames.dvrStringProperty, {
                              font: new PhetFont(14),
                            }),
                            new Text(
                              numericalMethodDescriptions.dvrStringProperty,
                              {
                                font: new PhetFont(11),
                                fill: "rgb(80,80,80)",
                                maxWidth: 550,
                              },
                            ),
                          ],
                        }),
                      tandemName: "dvrRadioButton",
                    },
                    {
                      value: NumericalMethod.FGH,
                      createNode: () =>
                        new VBox({
                          align: "left",
                          spacing: 4,
                          children: [
                            new Text(numericalMethodNames.fghStringProperty, {
                              font: new PhetFont(14),
                            }),
                            new Text(
                              numericalMethodDescriptions.fghStringProperty,
                              {
                                font: new PhetFont(11),
                                fill: "rgb(80,80,80)",
                                maxWidth: 550,
                              },
                            ),
                          ],
                        }),
                      tandemName: "fghRadioButton",
                    },
                    {
                      value: NumericalMethod.SPECTRAL,
                      createNode: () =>
                        new VBox({
                          align: "left",
                          spacing: 4,
                          children: [
                            new Text(
                              numericalMethodNames.spectralStringProperty,
                              {
                                font: new PhetFont(14),
                              },
                            ),
                            new Text(
                              numericalMethodDescriptions.spectralStringProperty,
                              {
                                font: new PhetFont(11),
                                fill: "rgb(80,80,80)",
                                maxWidth: 550,
                              },
                            ),
                          ],
                        }),
                      tandemName: "spectralRadioButton",
                    },
                    {
                      value: NumericalMethod.QUANTUM_BOUND,
                      createNode: () =>
                        new VBox({
                          align: "left",
                          spacing: 4,
                          children: [
                            new Text(
                              numericalMethodNames.quantumBoundStringProperty,
                              {
                                font: new PhetFont(14),
                              },
                            ),
                            new Text(
                              numericalMethodDescriptions.quantumBoundStringProperty,
                              {
                                font: new PhetFont(11),
                                fill: "rgb(80,80,80)",
                                maxWidth: 550,
                              },
                            ),
                          ],
                        }),
                      tandemName: "quantumBoundRadioButton",
                    },
                  ],
                  {
                    spacing: 12,
                    radioButtonOptions: {
                      radius: 8,
                    },
                  },
                );

              const numericalMethodSectionVBox = new VBox({
                align: "left",
                spacing: 12,
                children: [
                  new Text(preferencesLabels.numericalMethodStringProperty, {
                    font: new PhetFont({ size: 16, weight: "bold" }),
                  }),
                  new Text(
                    preferencesLabels.numericalMethodDescriptionStringProperty,
                    {
                      font: new PhetFont(12),
                      maxWidth: 600,
                    },
                  ),
                  numericalMethodRadioButtonGroup,
                ],
              });

              // Grid points slider section - uses exponent (6-11) for powers of 2
              // Ensure initial value is a power of 2
              const initialGridPoints =
                QPPWPreferences.gridPointsProperty.value;
              const initialExponent = Math.round(Math.log2(initialGridPoints));
              const correctedGridPoints = Math.pow(
                2,
                Math.max(5, Math.min(9, initialExponent)),
              );

              // Update gridPointsProperty if it wasn't a power of 2
              if (
                QPPWPreferences.gridPointsProperty.value !== correctedGridPoints
              ) {
                QPPWPreferences.gridPointsProperty.value = correctedGridPoints;
              }

              // Create a NumberProperty for the exponent that derives from gridPointsProperty
              const exponentProperty = new NumberProperty(
                Math.log2(correctedGridPoints),
                {
                  range: new Range(5, 9),
                },
              );

              // Track whether we're in the middle of updating to prevent infinite loops
              let isUpdating = false;

              // Bidirectional sync between exponent and grid points
              // Use lazyLink to avoid firing during initialization
              exponentProperty.lazyLink((exponent: number) => {
                if (isUpdating) return;
                isUpdating = true;
                const gridPoints = Math.pow(2, Math.round(exponent));
                if (QPPWPreferences.gridPointsProperty.value !== gridPoints) {
                  QPPWPreferences.gridPointsProperty.value = gridPoints;
                }
                isUpdating = false;
              });

              QPPWPreferences.gridPointsProperty.lazyLink(
                (gridPoints: number) => {
                  if (isUpdating) return;
                  isUpdating = true;
                  const exponent = Math.log2(gridPoints);
                  if (!isNaN(exponent) && exponentProperty.value !== exponent) {
                    exponentProperty.value = exponent;
                  }
                  isUpdating = false;
                },
              );

              const gridPointsSlider = new HSlider(
                exponentProperty,
                exponentProperty.range,
                {
                  trackSize: new Dimension2(400, 5),
                  thumbSize: new Dimension2(20, 40),
                  majorTickLength: 15,
                  minorTickLength: 10,
                  constrainValue: (value: number) => Math.round(value),
                },
              );

              // Add major ticks with actual grid point values
              gridPointsSlider.addMajorTick(
                5,
                new Text("32", { font: new PhetFont(12) }),
              );
              gridPointsSlider.addMajorTick(
                6,
                new Text("64", { font: new PhetFont(12) }),
              );
              gridPointsSlider.addMajorTick(
                7,
                new Text("128", { font: new PhetFont(12) }),
              );
              gridPointsSlider.addMajorTick(
                8,
                new Text("256", { font: new PhetFont(12) }),
              );
              gridPointsSlider.addMajorTick(
                9,
                new Text("512", { font: new PhetFont(12) }),
              );

              const gridPointsValueText = new Text("", {
                font: new PhetFont({ size: 14, weight: "bold" }),
              });

              QPPWPreferences.gridPointsProperty.link((value: number) => {
                gridPointsValueText.string = `${value} points`;
              });

              const gridPointsSectionVBox = new VBox({
                align: "left",
                spacing: 12,
                children: [
                  new Text(preferencesLabels.gridPointsStringProperty, {
                    font: new PhetFont({ size: 16, weight: "bold" }),
                  }),
                  new Text(
                    preferencesLabels.gridPointsDescriptionStringProperty,
                    {
                      font: new PhetFont(12),
                      maxWidth: 600,
                    },
                  ),
                  new HBox({
                    spacing: 15,
                    children: [gridPointsSlider, gridPointsValueText],
                  }),
                ],
              });

              return new VBox({
                align: "left",
                spacing: 20,
                children: [
                  autoPauseSectionVBox,
                  new HStrut(650), // Set minimum width
                  numericalMethodSectionVBox,
                  gridPointsSectionVBox,
                ],
              });
            },
          },
        ],
      },
    }),
  };

  const screens = [
    new IntroScreen({
      name: screenNames.introStringProperty,
      tandem: Tandem.ROOT.createTandem("introScreen"),
      backgroundColorProperty: QPPWColors.backgroundColorProperty,
      homeScreenIcon: new IntroScreenIcon(),
    }),
    new OneWellScreen({
      name: screenNames.oneWellStringProperty,
      tandem: Tandem.ROOT.createTandem("oneWellScreen"),
      backgroundColorProperty: QPPWColors.backgroundColorProperty,
      homeScreenIcon: new OneWellScreenIcon(),
    }),
    new TwoWellsScreen({
      name: screenNames.twoWellsStringProperty,
      tandem: Tandem.ROOT.createTandem("twoWellsScreen"),
      backgroundColorProperty: QPPWColors.backgroundColorProperty,
      homeScreenIcon: new TwoWellsScreenIcon(),
    }),
    new ManyWellsScreen({
      name: screenNames.manyWellsStringProperty,
      tandem: Tandem.ROOT.createTandem("manyWellsScreen"),
      backgroundColorProperty: QPPWColors.backgroundColorProperty,
      homeScreenIcon: new ManyWellsScreenIcon(),
    }),
  ];

  const sim = new Sim(
    stringManager.getTitleStringProperty(),
    screens,
    simOptions,
  );
  sim.start();
});
