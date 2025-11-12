/**
 * OneWellScreenIcon provides the icon for the One Well screen.
 */

import { Rectangle } from "scenerystack/scenery";
import { ScreenIcon } from "scenerystack/sim";
import QPPWColors from "../../QPPWColors.js";

export class OneWellScreenIcon extends ScreenIcon {
  public constructor() {
    // Create a simple icon showing a single potential well
    const iconNode = new Rectangle(0, 0, 50, 50, {
      fill: QPPWColors.backgroundColorProperty,
      stroke: QPPWColors.potentialWellProperty,
      lineWidth: 3,
    });

    super(iconNode, {
      fill: QPPWColors.backgroundColorProperty,
    });
  }
}
