// NOTE: brand.js needs to be the first import. This is because SceneryStack for sims needs a very specific loading
// order: init.ts => assert.ts => splash.ts => brand.ts => everything else (here)
import "./brand.js";

import { onReadyToLaunch, Sim, PreferencesModel } from "scenerystack/sim";
import { Tandem } from "scenerystack/tandem";
import { SimScreen } from "./screen-name/SimScreen.js";
import { StringManager } from "./i18n/StringManager.js";
import QPPWPreferences from "./QPPWPreferences.js";
import QPPWColors from "./QPPWColors.js";

onReadyToLaunch(() => {
  // Get the string manager instance
  const stringManager = StringManager.getInstance();
  const screenNames = stringManager.getScreenNames();

  // The title uses localized strings
  const titleStringProperty = screenNames.oneWellScreen; // Will be "QPPW" or full title

  // TODO: Create three specific screen classes:
  // - OneWellScreen (extending Screen with OneWellModel and OneWellScreenView)
  // - TwoWellsScreen (extending Screen with TwoWellsModel and TwoWellsScreenView)
  // - ManyWellsScreen (extending Screen with ManyWellsModel and ManyWellsScreenView)
  // Each ScreenView should extend BaseScreenView from common/view/BaseScreenView.ts

  const screens = [
    new SimScreen({ tandem: Tandem.ROOT.createTandem("oneWellScreen") }),
    // new TwoWellsScreen({ tandem: Tandem.ROOT.createTandem("twoWellsScreen") }),
    // new ManyWellsScreen({ tandem: Tandem.ROOT.createTandem("manyWellsScreen") }),
  ];

  // Create preferences model
  const preferencesModel = new PreferencesModel({
    simulationOptions: {
      customPreferences: [
        {
          createContent: () => {
            // TODO: Create preferences panel UI using QPPWPreferences
            // This should include toggles for:
            // - autoPauseWhenTabHiddenProperty
            // - announceParameterChangesProperty
            // - announceStateChangesProperty
            // - announceDragInteractionsProperty
            return null;
          },
        },
      ],
    },
  });

  const sim = new Sim(titleStringProperty, screens, {
    preferencesModel: preferencesModel,
    // Enable WebGL for better performance
    webgl: true,
  });

  sim.start();
});
