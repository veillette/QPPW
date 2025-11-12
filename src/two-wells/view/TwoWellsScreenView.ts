/**
 * TwoWellsScreenView is the main view for the Two Wells screen.
 * It displays a double quantum potential well and demonstrates quantum tunneling.
 */

import { BaseScreenView } from "../../common/view/BaseScreenView.js";
import { TwoWellsModel } from "../model/TwoWellsModel.js";
import { ScreenViewOptions } from "scenerystack/sim";
import { Node, Rectangle, Text, VBox } from "scenerystack/scenery";
import { Panel } from "scenerystack/sun";
import QPPWColors from "../../QPPWColors.js";
import stringManager from "../../i18n/StringManager.js";

export class TwoWellsScreenView extends BaseScreenView {
  private readonly leftWellNode: Rectangle;
  private readonly rightWellNode: Rectangle;
  private readonly barrierNode: Rectangle;
  private readonly controlPanel: Panel;

  public constructor(model: TwoWellsModel, options?: ScreenViewOptions) {
    super(
      () => {
        model.reset();
        this.reset();
      },
      options,
    );

    const wellWidth = 150;
    const wellHeight = 300;
    const barrierWidth = 50;
    const centerX = this.layoutBounds.centerX - 100;
    const centerY = this.layoutBounds.centerY;

    // Create the left potential well
    this.leftWellNode = new Rectangle(0, 0, wellWidth, wellHeight, {
      fill: QPPWColors.backgroundColorProperty,
      stroke: QPPWColors.potentialWellProperty,
      lineWidth: 3,
      centerX: centerX - barrierWidth / 2 - wellWidth / 2,
      centerY: centerY,
    });
    this.addChild(this.leftWellNode);

    // Create the barrier between wells
    this.barrierNode = new Rectangle(0, 0, barrierWidth, wellHeight, {
      fill: QPPWColors.potentialBarrierProperty,
      opacity: 0.5,
      centerX: centerX,
      centerY: centerY,
    });
    this.addChild(this.barrierNode);

    // Create the right potential well
    this.rightWellNode = new Rectangle(0, 0, wellWidth, wellHeight, {
      fill: QPPWColors.backgroundColorProperty,
      stroke: QPPWColors.potentialWellProperty,
      lineWidth: 3,
      centerX: centerX + barrierWidth / 2 + wellWidth / 2,
      centerY: centerY,
    });
    this.addChild(this.rightWellNode);

    // Add title text
    const titleText = new Text(stringManager.twoWellsStringProperty, {
      font: "24px sans-serif",
      fill: QPPWColors.textFillProperty,
      centerX: centerX,
      top: 20,
    });
    this.addChild(titleText);

    // Create control panel
    this.controlPanel = this.createControlPanel();
    this.addChild(this.controlPanel);

    // Add placeholder content text
    const contentText = new Text("Double well with quantum tunneling", {
      font: "16px sans-serif",
      fill: QPPWColors.labelFillProperty,
      centerX: centerX,
      top: titleText.bottom + 20,
    });
    this.addChild(contentText);
  }

  /**
   * Creates the control panel for adjusting well and barrier parameters.
   */
  private createControlPanel(): Panel {
    const content = new VBox({
      spacing: 10,
      align: "left",
      children: [
        new Text(stringManager.wellWidthStringProperty, {
          font: "14px sans-serif",
          fill: QPPWColors.textFillProperty,
        }),
        new Text(stringManager.barrierHeightStringProperty, {
          font: "14px sans-serif",
          fill: QPPWColors.textFillProperty,
        }),
        new Text(stringManager.barrierWidthStringProperty, {
          font: "14px sans-serif",
          fill: QPPWColors.textFillProperty,
        }),
        new Text(stringManager.tunnelingStringProperty, {
          font: "14px sans-serif",
          fill: QPPWColors.textFillProperty,
        }),
      ],
    });

    return new Panel(content, {
      fill: QPPWColors.panelFillProperty,
      stroke: QPPWColors.panelStrokeProperty,
      xMargin: 10,
      yMargin: 10,
      right: this.layoutBounds.maxX - 10,
      top: 20,
    });
  }

  /**
   * Creates the content for the info dialog.
   */
  public createInfoDialogContent(): Node {
    return new Text(
      "Explore quantum tunneling in a double potential well.\n" +
        "Adjust barrier parameters to see how tunneling probability changes.\n" +
        "Watch particles tunnel through classically forbidden regions!",
      {
        font: "14px sans-serif",
        fill: QPPWColors.textFillProperty,
      },
    );
  }

  /**
   * Creates the screen summary content for accessibility.
   */
  public createScreenSummaryContent(): Node {
    return new Text(
      "Two Wells screen demonstrates quantum tunneling between two potential wells.",
      {
        font: "14px sans-serif",
        fill: QPPWColors.textFillProperty,
      },
    );
  }

  /**
   * Resets the screen view to its initial state.
   */
  public override reset(): void {
    super.reset();
    // Add screen-specific reset logic here
  }

  /**
   * Steps the screen view forward in time.
   * @param dt - The time step in seconds
   */
  public override step(dt: number): void {
    super.step(dt);
    // Add animation/update logic here
  }
}
