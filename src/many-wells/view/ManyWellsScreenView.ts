/**
 * ManyWellsScreenView is the main view for the Many Wells screen.
 * It displays multiple quantum potential wells and demonstrates energy band formation.
 */

import { BaseScreenView } from "../../common/view/BaseScreenView.js";
import { ManyWellsModel } from "../model/ManyWellsModel.js";
import { ScreenViewOptions } from "scenerystack/sim";
import { Node, Rectangle, Text, VBox, RichText } from "scenerystack/scenery";
import { Panel } from "scenerystack/sun";
import { PhetFont } from "scenerystack/scenery-phet";
import QPPWColors from "../../QPPWColors.js";
import stringManager from "../../i18n/StringManager.js";

export class ManyWellsScreenView extends BaseScreenView {
  private readonly wellsContainer: Node;
  private readonly customControlPanel: Panel;
  private readonly _model: ManyWellsModel;

  public constructor(model: ManyWellsModel, options?: ScreenViewOptions) {
    super(model, options);

    this._model = model;

    // Create container for the wells
    this.wellsContainer = new Node();
    this.addChild(this.wellsContainer);

    // Create the initial wells array
    this.createWellsVisualization();

    // Add title text
    const titleText = new Text(stringManager.manyWellsStringProperty, {
      font: new PhetFont(24),
      fill: QPPWColors.textFillProperty,
      centerX: this.layoutBounds.centerX,
      top: 20,
    });
    this.addChild(titleText);

    // Create control panel
    this.customControlPanel = this.createControlPanel();
    this.addChild(this.customControlPanel);

    // Add placeholder content text
    const contentText = new Text(stringManager.multipleWellsStringProperty, {
      font: new PhetFont(16),
      fill: QPPWColors.labelFillProperty,
      centerX: this.layoutBounds.centerX - 100,
      top: titleText.bottom + 20,
    });
    this.addChild(contentText);

    // Listen to numberOfWellsProperty changes
    model.numberOfWellsProperty.link(() => {
      this.updateWellsVisualization();
    });
  }

  /**
   * Creates the visualization of multiple potential wells.
   */
  private createWellsVisualization(): void {
    const numWells = this._model.numberOfWellsProperty.value;
    const wellWidth = 40;
    const wellHeight = 200;
    const barrierWidth = 15;
    const spacing = wellWidth + barrierWidth;

    const totalWidth = numWells * wellWidth + (numWells - 1) * barrierWidth;
    const startX = this.layoutBounds.centerX - totalWidth / 2 - 100;
    const centerY = this.layoutBounds.centerY;

    for (let i = 0; i < numWells; i++) {
      // Create well
      const well = new Rectangle(0, 0, wellWidth, wellHeight, {
        fill: QPPWColors.backgroundColorProperty,
        stroke: QPPWColors.potentialWellProperty,
        lineWidth: 2,
        left: startX + i * spacing,
        centerY: centerY,
      });
      this.wellsContainer.addChild(well);

      // Create barrier (except after the last well)
      if (i < numWells - 1) {
        const barrier = new Rectangle(0, 0, barrierWidth, wellHeight, {
          fill: QPPWColors.potentialBarrierProperty,
          opacity: 0.4,
          left: startX + i * spacing + wellWidth,
          centerY: centerY,
        });
        this.wellsContainer.addChild(barrier);
      }
    }
  }

  /**
   * Updates the wells visualization when parameters change.
   */
  private updateWellsVisualization(): void {
    this.wellsContainer.removeAllChildren();
    this.createWellsVisualization();
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
          font: new PhetFont(14),
          fill: QPPWColors.textFillProperty,
        }),
        new Text(stringManager.latticeConstantStringProperty, {
          font: new PhetFont(14),
          fill: QPPWColors.textFillProperty,
        }),
        new Text(stringManager.energyBandsStringProperty, {
          font: new PhetFont(14),
          fill: QPPWColors.textFillProperty,
        }),
        new Text(stringManager.wellWidthStringProperty, {
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
    const titleText = new Text(stringManager.manyWellsStringProperty, {
      font: new PhetFont({ size: 18, weight: "bold" }),
      fill: QPPWColors.textFillProperty,
    });

    const descriptionText = new RichText(stringManager.manyWellsDescriptionStringProperty, {
      font: new PhetFont(14),
      fill: QPPWColors.textFillProperty,
      maxWidth: 500,
    });

    const keyConceptsTitle = new Text(stringManager.keyConceptsTitleStringProperty, {
      font: new PhetFont({ size: 14, weight: "bold" }),
      fill: QPPWColors.textFillProperty,
    });

    const keyConceptsList = new RichText(
      stringManager.manyWellsKeyConceptsStringProperty,
      {
        font: new PhetFont(13),
        fill: QPPWColors.textFillProperty,
        maxWidth: 500,
      }
    );

    const interactionTitle = new Text(stringManager.interactionsTitleStringProperty, {
      font: new PhetFont({ size: 14, weight: "bold" }),
      fill: QPPWColors.textFillProperty,
    });

    const interactionsList = new RichText(
      stringManager.manyWellsInteractionsStringProperty,
      {
        font: new PhetFont(13),
        fill: QPPWColors.textFillProperty,
        maxWidth: 500,
      }
    );

    return new VBox({
      spacing: 12,
      align: "left",
      children: [
        titleText,
        descriptionText,
        keyConceptsTitle,
        keyConceptsList,
        interactionTitle,
        interactionsList,
      ],
    });
  }

  /**
   * Creates the screen summary content for accessibility.
   */
  public createScreenSummaryContent(): Node {
    const summaryText = new RichText(
      stringManager.manyWellsEducationalContentStringProperty,
      {
        font: new PhetFont(13),
        fill: QPPWColors.textFillProperty,
        maxWidth: 600,
      }
    );

    return new VBox({
      spacing: 10,
      align: "left",
      children: [summaryText],
    });
  }

  /**
   * Resets the screen view to its initial state.
   */
  public override reset(): void {
    super.reset();
    this.updateWellsVisualization();
  }

  /**
   * Steps the screen view forward in time.
   * @param dt - The time step in seconds
   */
  public override step(dt: number): void {
    super.step(dt);
    this._model.step(dt);
  }
}
