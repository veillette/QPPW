/**
 * TwoWellsScreen represents the screen for exploring quantum tunneling in a double potential well.
 */

import { Screen, ScreenOptions } from "scenerystack/sim";
import { TwoWellsModel } from "./model/TwoWellsModel.js";
import { TwoWellsScreenView } from "./view/TwoWellsScreenView.js";

export class TwoWellsScreen extends Screen<TwoWellsModel, TwoWellsScreenView> {
  public constructor(options: ScreenOptions) {
    super(
      () => new TwoWellsModel(),
      (model: TwoWellsModel) => new TwoWellsScreenView(model),
      options,
    );
  }
}
