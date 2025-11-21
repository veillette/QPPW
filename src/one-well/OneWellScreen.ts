/**
 * OneWellScreen represents the screen for exploring a single quantum potential well.
 */

import { Screen, ScreenOptions } from "scenerystack/sim";
import { OneWellModel } from "./model/OneWellModel.js";
import { OneWellScreenView } from "./view/OneWellScreenView.js";

export class OneWellScreen extends Screen<OneWellModel, OneWellScreenView> {
  public constructor(options: ScreenOptions) {
    super(
      () => new OneWellModel(),
      (model: OneWellModel) => new OneWellScreenView(model),
      options,
    );
  }
}
