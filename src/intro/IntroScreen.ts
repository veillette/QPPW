/**
 * IntroScreen represents the introductory screen for exploring quantum potential wells.
 * This screen provides a simplified interface without play/pause controls or superposition options.
 */

import { Screen, ScreenOptions } from "scenerystack/sim";
import { IntroModel } from "./model/IntroModel.js";
import { IntroScreenView } from "./view/IntroScreenView.js";

export class IntroScreen extends Screen<IntroModel, IntroScreenView> {
    public constructor(options: ScreenOptions) {
        super(
            () => new IntroModel(),
            (model: IntroModel) => new IntroScreenView(model),
            options,
        );
    }
}
