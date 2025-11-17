// NOTE: brand.js needs to be the first import. This is because SceneryStack for sims needs a very specific loading
// order: init.ts => assert.ts => splash.ts => brand.ts => everything else (here)
import "./brand.js";

import { onReadyToLaunch, Sim, PreferencesModel } from "scenerystack/sim";
import { Tandem } from "scenerystack/tandem";
import { VBox, Text, HStrut, HBox } from "scenerystack/scenery";
import { Checkbox, VerticalAquaRadioButtonGroup, HSlider } from "scenerystack/sun";
import { PhetFont } from "scenerystack/scenery-phet";
import { Dimension2 } from "scenerystack";
import stringManager from "./i18n/StringManager.js";
import QPPWColors from "./QPPWColors.js";
import QPPWPreferences from "./QPPWPreferences.js";
import { NumericalMethod } from "./common/model/Schrodinger1DSolver.js";
import { OneWellScreen } from "./one-well/OneWellScreen.js";
import { TwoWellsScreen } from "./two-wells/TwoWellsScreen.js";
import { ManyWellsScreen } from "./many-wells/ManyWellsScreen.js";
import { OneWellScreenIcon } from "./one-well/view/OneWellScreenIcon.js";
import { TwoWellsScreenIcon } from "./two-wells/view/TwoWellsScreenIcon.js";
import { ManyWellsScreenIcon } from "./many-wells/view/ManyWellsScreenIcon.js";

onReadyToLaunch(() => {
  const screenNames = stringManager.getScreenNames();

  // Get preferences strings
  const preferencesLabels = stringManager.getPreferencesLabels();
  const numericalMethodNames = stringManager.getNumericalMethodNames();
  const numericalMethodDescriptions = stringManager.getNumericalMethodDescriptions();

  const simOptions = {
    webgl: true,
    preferencesModel: new PreferencesModel({
      visualOptions: {
        supportsProjectorMode: true,
        supportsInteractiveHighlights: true,
      },
      simulationOptions: {
        customPreferences: [
          {
            // eslint-disable-next-line @typescript-eslint/no-unused-vars
            createContent: (_tandem: Tandem) => {
              // Auto-pause preference
              const autoPauseSection = new VBox({
                align: "left",
                spacing: 8,
                children: [
                  new Checkbox(
                    QPPWPreferences.autoPauseWhenTabHiddenProperty,
                    new Text(preferencesLabels.autoPauseWhenTabHiddenStringProperty, {
                      font: new PhetFont(16),
                    }),
                    {
                      boxWidth: 16,
                    },
                  ),
                  new Text(preferencesLabels.autoPauseDescriptionStringProperty, {
                    font: new PhetFont(12),
                    maxWidth: 600,
                  }),
                ],
              });

              // Numerical method section
              const numericalMethodRadioButtonGroup = new VerticalAquaRadioButtonGroup(
                QPPWPreferences.numericalMethodProperty,
                [
                  {
                    value: NumericalMethod.NUMEROV,
                    createNode: () => new VBox({
                      align: "left",
                      spacing: 4,
                      children: [
                        new Text(numericalMethodNames.numerovStringProperty, {
                          font: new PhetFont(14),
                        }),
                        new Text(numericalMethodDescriptions.numerovStringProperty, {
                          font: new PhetFont(11),
                          fill: "rgb(80,80,80)",
                          maxWidth: 550,
                        }),
                      ],
                    }),
                    tandemName: "numerovRadioButton",
                  },
                  {
                    value: NumericalMethod.MATRIX_NUMEROV,
                    createNode: () => new VBox({
                      align: "left",
                      spacing: 4,
                      children: [
                        new Text(numericalMethodNames.matrixNumerovStringProperty, {
                          font: new PhetFont(14),
                        }),
                        new Text(numericalMethodDescriptions.matrixNumerovStringProperty, {
                          font: new PhetFont(11),
                          fill: "rgb(80,80,80)",
                          maxWidth: 550,
                        }),
                      ],
                    }),
                    tandemName: "matrixNumerovRadioButton",
                  },
                  {
                    value: NumericalMethod.DVR,
                    createNode: () => new VBox({
                      align: "left",
                      spacing: 4,
                      children: [
                        new Text(numericalMethodNames.dvrStringProperty, {
                          font: new PhetFont(14),
                        }),
                        new Text(numericalMethodDescriptions.dvrStringProperty, {
                          font: new PhetFont(11),
                          fill: "rgb(80,80,80)",
                          maxWidth: 550,
                        }),
                      ],
                    }),
                    tandemName: "dvrRadioButton",
                  },
                  {
                    value: NumericalMethod.FGH,
                    createNode: () => new VBox({
                      align: "left",
                      spacing: 4,
                      children: [
                        new Text(numericalMethodNames.fghStringProperty, {
                          font: new PhetFont(14),
                        }),
                        new Text(numericalMethodDescriptions.fghStringProperty, {
                          font: new PhetFont(11),
                          fill: "rgb(80,80,80)",
                          maxWidth: 550,
                        }),
                      ],
                    }),
                    tandemName: "fghRadioButton",
                  },
                  {
                    value: NumericalMethod.SPECTRAL,
                    createNode: () => new VBox({
                      align: "left",
                      spacing: 4,
                      children: [
                        new Text(numericalMethodNames.spectralStringProperty, {
                          font: new PhetFont(14),
                        }),
                        new Text(numericalMethodDescriptions.spectralStringProperty, {
                          font: new PhetFont(11),
                          fill: "rgb(80,80,80)",
                          maxWidth: 550,
                        }),
                      ],
                    }),
                    tandemName: "spectralRadioButton",
                  },
                ],
                {
                  spacing: 12,
                  radioButtonOptions: {
                    radius: 8,
                  },
                },
              );

              const numericalMethodSection = new VBox({
                align: "left",
                spacing: 12,
                children: [
                  new Text(preferencesLabels.numericalMethodStringProperty, {
                    font: new PhetFont({ size: 16, weight: "bold" }),
                  }),
                  new Text(preferencesLabels.numericalMethodDescriptionStringProperty, {
                    font: new PhetFont(12),
                    maxWidth: 600,
                  }),
                  numericalMethodRadioButtonGroup,
                ],
              });

              // Grid points slider section
              const gridPointsSlider = new HSlider(
                QPPWPreferences.gridPointsProperty,
                QPPWPreferences.gridPointsProperty.range,
                {
                  trackSize: new Dimension2(400, 5),
                  thumbSize: new Dimension2(20, 40),
                  majorTickLength: 15,
                  minorTickLength: 10,
                },
              );

              // Add major ticks
              gridPointsSlider.addMajorTick(256, new Text("256", { font: new PhetFont(12) }));
              gridPointsSlider.addMajorTick(512, new Text("512", { font: new PhetFont(12) }));
              gridPointsSlider.addMajorTick(1024, new Text("1024", { font: new PhetFont(12) }));
              gridPointsSlider.addMajorTick(2000, new Text("2000", { font: new PhetFont(12) }));

              // Add minor ticks
              for (let i = 300; i <= 2000; i += 100) {
                if (i % 500 !== 0) {
                  gridPointsSlider.addMinorTick(i);
                }
              }

              const gridPointsValueText = new Text("", {
                font: new PhetFont({ size: 14, weight: "bold" }),
              });

              QPPWPreferences.gridPointsProperty.link((value) => {
                gridPointsValueText.string = `${value} points`;
              });

              const gridPointsSection = new VBox({
                align: "left",
                spacing: 12,
                children: [
                  new Text(preferencesLabels.gridPointsStringProperty, {
                    font: new PhetFont({ size: 16, weight: "bold" }),
                  }),
                  new Text(preferencesLabels.gridPointsDescriptionStringProperty, {
                    font: new PhetFont(12),
                    maxWidth: 600,
                  }),
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
                  autoPauseSection,
                  new HStrut(650), // Set minimum width
                  numericalMethodSection,
                  gridPointsSection,
                ],
              });
            },
          },
        ],
      },
    }),
  };

  const screens = [
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
