// NOTE: brand.js needs to be the first import. This is because SceneryStack for sims needs a very specific loading
// order: init.ts => assert.ts => splash.ts => brand.ts => everything else (here)
import "./brand.js";

import { onReadyToLaunch, Sim, PreferencesModel } from "scenerystack/sim";
import { Tandem } from "scenerystack/tandem";
import { VBox, Text, HStrut } from "scenerystack/scenery";
import { Checkbox, VerticalAquaRadioButtonGroup } from "scenerystack/sun";
import { PhetFont } from "scenerystack/scenery-phet";
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

              return new VBox({
                align: "left",
                spacing: 20,
                children: [
                  autoPauseSection,
                  new HStrut(650), // Set minimum width
                  numericalMethodSection,
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
