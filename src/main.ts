// NOTE: brand.js needs to be the first import. This is because SceneryStack for sims needs a very specific loading
// order: init.ts => assert.ts => splash.ts => brand.ts => everything else (here)
import "./brand.js";

import { onReadyToLaunch, Sim } from "scenerystack/sim";
import { Tandem } from "scenerystack/tandem";
import stringManager from "./i18n/StringManager.js";
import QPPWColors from "./QPPWColors.js";
import { OneWellScreen } from "./one-well/OneWellScreen.js";
import { TwoWellsScreen } from "./two-wells/TwoWellsScreen.js";
import { ManyWellsScreen } from "./many-wells/ManyWellsScreen.js";
import { OneWellScreenIcon } from "./one-well/view/OneWellScreenIcon.js";
import { TwoWellsScreenIcon } from "./two-wells/view/TwoWellsScreenIcon.js";
import { ManyWellsScreenIcon } from "./many-wells/view/ManyWellsScreenIcon.js";

onReadyToLaunch(() => {
  const screenNames = stringManager.getScreenNames();

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
  );
  sim.start();
});
