/**
 * ManyWellsScreen represents the screen for exploring energy bands in periodic potential wells.
 */

import { Screen, ScreenOptions } from "scenerystack/sim";
import { ManyWellsModel } from "./model/ManyWellsModel.js";
import { ManyWellsScreenView } from "./view/ManyWellsScreenView.js";

export class ManyWellsScreen extends Screen<
  ManyWellsModel,
  ManyWellsScreenView
> {
  public constructor(options: ScreenOptions) {
    super(
      () => new ManyWellsModel(),
      (model: ManyWellsModel) => new ManyWellsScreenView(model),
      options,
    );
  }
}
