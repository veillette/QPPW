/**
 * OneWellScreenView is the main view for the One Well screen.
 * It displays a single quantum potential well and its wave functions.
 */

import { BaseScreenView } from "../../common/view/BaseScreenView.js";
import { OneWellModel } from "../model/OneWellModel.js";
import { ScreenViewOptions } from "scenerystack/sim";
import { Node, Rectangle, Text, VBox } from "scenerystack/scenery";
import { Panel } from "scenerystack/sun";
import QPPWColors from "../../QPPWColors.js";
import stringManager from "../../i18n/StringManager.js";

export class OneWellScreenView extends BaseScreenView {
  private readonly potentialWellNode: Rectangle;
  private readonly controlPanel: Panel;

  public constructor(model: OneWellModel, options?: ScreenViewOptions) {
    super(
      () => {
        model.reset();
        this.reset();
      },
      options,
    );

    // Create the potential well visualization
    this.potentialWellNode = new Rectangle(0, 0, 400, 300, {
      fill: QPPWColors.backgroundColorProperty,
      stroke: QPPWColors.potentialWellProperty,
      lineWidth: 3,
      centerX: this.layoutBounds.centerX - 100,
      centerY: this.layoutBounds.centerY,
    });
    this.addChild(this.potentialWellNode);

    // Add title text
    const titleText = new Text(stringManager.oneWellStringProperty, {
      font: "24px sans-serif",
      fill: QPPWColors.textFillProperty,
      centerX: this.layoutBounds.centerX - 100,
      top: 20,
    });
    this.addChild(titleText);

    // Create control panel
    this.controlPanel = this.createControlPanel();
    this.addChild(this.controlPanel);

    // Add placeholder content text
    const contentText = new Text("Single potential well visualization", {
      font: "16px sans-serif",
      fill: QPPWColors.labelFillProperty,
      centerX: this.potentialWellNode.centerX,
      centerY: this.potentialWellNode.centerY,
    });
    this.addChild(contentText);
  }

  /**
   * Creates the control panel for adjusting well parameters.
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
        new Text(stringManager.wellDepthStringProperty, {
          font: "14px sans-serif",
          fill: QPPWColors.textFillProperty,
        }),
        new Text(stringManager.energyStringProperty, {
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
      "Explore quantum mechanics in a single potential well.\n" +
        "Adjust the well parameters to see how energy levels change.",
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
      "One Well screen shows a single quantum potential well with adjustable parameters.",
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
