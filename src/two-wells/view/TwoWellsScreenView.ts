/**
 * TwoWellsScreenView is the main view for the Two Wells screen.
 * It displays a double quantum potential well and demonstrates quantum tunneling.
 */

import { BaseScreenView } from "../../common/view/BaseScreenView.js";
import { TwoWellsModel } from "../model/TwoWellsModel.js";
import { ScreenViewOptions } from "scenerystack/sim";
import { Node, Rectangle, Text, VBox } from "scenerystack/scenery";
import { Panel } from "scenerystack/sun";
import { PhetFont } from "scenerystack/scenery-phet";
import QPPWColors from "../../QPPWColors.js";
import stringManager from "../../i18n/StringManager.js";

export class TwoWellsScreenView extends BaseScreenView {
  private readonly leftWellNode: Rectangle;
  private readonly rightWellNode: Rectangle;
  private readonly barrierNode: Rectangle;
  private readonly customControlPanel: Panel;

  public constructor(model: TwoWellsModel, options?: ScreenViewOptions) {
    super(model, options);

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
      font: new PhetFont(24),
      fill: QPPWColors.textFillProperty,
      centerX: centerX,
      top: 20,
    });
    this.addChild(titleText);

    // Create control panel
    this.customControlPanel = this.createControlPanel();
    this.addChild(this.customControlPanel);

    // Add placeholder content text
    const contentText = new Text(stringManager.doubleWellStringProperty, {
      font: new PhetFont(16),
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
          font: new PhetFont(14),
          fill: QPPWColors.textFillProperty,
        }),
        new Text(stringManager.barrierHeightStringProperty, {
          font: new PhetFont(14),
          fill: QPPWColors.textFillProperty,
        }),
        new Text(stringManager.barrierWidthStringProperty, {
          font: new PhetFont(14),
          fill: QPPWColors.textFillProperty,
        }),
        new Text(stringManager.tunnelingStringProperty, {
          font: new PhetFont(14),
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
    return new Text(stringManager.twoWellsDescriptionStringProperty, {
      font: new PhetFont(14),
      fill: QPPWColors.textFillProperty,
    });
  }

  /**
   * Creates the screen summary content for accessibility.
   */
  public createScreenSummaryContent(): Node {
    return new Text(stringManager.twoWellsSummaryStringProperty, {
      font: new PhetFont(14),
      fill: QPPWColors.textFillProperty,
    });
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
    (this.model as TwoWellsModel).step(dt);
  }
}
