/**
 * ManyWellsScreenIcon provides the icon for the Many Wells screen.
 */

import { Rectangle, Node, HBox } from "scenerystack/scenery";
import { ScreenIcon } from "scenerystack/sim";
import QPPWColors from "../../QPPWColors.js";

export class ManyWellsScreenIcon extends ScreenIcon {
  public constructor() {
    // Create a simple icon showing multiple potential wells
    const wells: Node[] = [];
    const numWells = 4;
    const wellWidth = 10;
    const barrierWidth = 4;

    for (let i = 0; i < numWells; i++) {
      // Add well
      wells.push(
        new Rectangle(0, 0, wellWidth, 50, {
          fill: QPPWColors.backgroundColorProperty,
          stroke: QPPWColors.potentialWellProperty,
          lineWidth: 2,
        }),
      );

      // Add barrier (except after last well)
      if (i < numWells - 1) {
        wells.push(
          new Rectangle(0, 0, barrierWidth, 50, {
            fill: QPPWColors.potentialBarrierProperty,
            opacity: 0.5,
          }),
        );
      }
    }

    const iconNode = new HBox({
      spacing: 0,
      children: wells,
    });

    super(iconNode, {
      fill: QPPWColors.backgroundColorProperty,
    });
  }
}
