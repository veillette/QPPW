/**
 * TwoWellsScreenIcon provides the icon for the Two Wells screen.
 */

import { Rectangle, HBox } from "scenerystack/scenery";
import { ScreenIcon } from "scenerystack/sim";
import QPPWColors from "../../QPPWColors.js";

export class TwoWellsScreenIcon extends ScreenIcon {
  public constructor() {
    // Create a simple icon showing two potential wells with a barrier
    const leftWell = new Rectangle(0, 0, 20, 50, {
      fill: QPPWColors.backgroundColorProperty,
      stroke: QPPWColors.potentialWellProperty,
      lineWidth: 2,
    });

    const barrier = new Rectangle(0, 0, 8, 50, {
      fill: QPPWColors.potentialBarrierProperty,
      opacity: 0.6,
    });

    const rightWell = new Rectangle(0, 0, 20, 50, {
      fill: QPPWColors.backgroundColorProperty,
      stroke: QPPWColors.potentialWellProperty,
      lineWidth: 2,
    });

    const iconNode = new HBox({
      spacing: 0,
      children: [leftWell, barrier, rightWell],
    });

    super(iconNode, {
      fill: QPPWColors.backgroundColorProperty,
    });
  }
}
