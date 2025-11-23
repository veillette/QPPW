/**
 * ManyWellsScreenView is the main view for the Many Wells screen.
 * It displays multiple quantum potential wells and demonstrates energy band formation.
 */

import { BaseScreenView } from "../../common/view/BaseScreenView.js";
import { ManyWellsModel } from "../model/ManyWellsModel.js";
import { ScreenViewOptions } from "scenerystack/sim";
import { Node, Text, VBox } from "scenerystack/scenery";
import { Panel } from "scenerystack/sun";
import QPPWColors from "../../QPPWColors.js";
import stringManager from "../../i18n/StringManager.js";

export class ManyWellsScreenView extends BaseScreenView {
  private readonly wellsContainer: Node;
  private readonly customControlPanel: Panel;

  public constructor(model: ManyWellsModel, options?: ScreenViewOptions) {
    super(model, options);

    // Create container for the wells
    this.wellsContainer = new Node();
    this.addChild(this.wellsContainer);

    // Add title text
    const titleText = new Text(stringManager.manyWellsStringProperty, {
      font: "24px sans-serif",
      fill: QPPWColors.textFillProperty,
      centerX: this.layoutBounds.centerX,
      top: 20,
    });
    this.addChild(titleText);

    // Create control panel
    this.customControlPanel = this.createControlPanel();
    this.addChild(this.customControlPanel);
  }

  /**
   * Creates the control panel for adjusting well parameters.
   */
  private createControlPanel(): Panel {
    const content = new VBox({
      spacing: 10,
      align: "left",
      children: [
        new Text(stringManager.numberOfWellsStringProperty, {
          font: "14px sans-serif",
          fill: QPPWColors.textFillProperty,
        }),
        new Text(stringManager.latticeConstantStringProperty, {
          font: "14px sans-serif",
          fill: QPPWColors.textFillProperty,
        }),
        new Text(stringManager.energyBandsStringProperty, {
          font: "14px sans-serif",
          fill: QPPWColors.textFillProperty,
        }),
        new Text(stringManager.wellWidthStringProperty, {
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
      "Explore energy bands in a periodic potential.\n" +
        "Add or remove wells to see how energy bands form.\n" +
        "This demonstrates the foundation of solid-state physics!",
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
      "Many Wells screen demonstrates energy band formation in periodic potentials.",
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
