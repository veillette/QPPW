/**
 * SimulationControlBar provides playback controls and time display.
 * This is displayed at the bottom of quantum physics screens.
 */

import { Node, Text, HBox, VBox, Rectangle } from "scenerystack/scenery";
import {
  TimeControlNode,
  PhetFont,
} from "scenerystack/scenery-phet";

import { DerivedProperty } from "scenerystack/axon";
import { BaseModel } from "../model/BaseModel.js";
import QPPWColors from "../../QPPWColors.js";
import stringManager from "../../i18n/StringManager.js";

export class SimulationControlBar extends Node {
  private readonly model: BaseModel;
  private readonly timeText: Text;
  private readonly barWidth: number;

  public constructor(model: BaseModel, options?: { width?: number }) {
    super();

    this.model = model;
    this.barWidth = options?.width ?? 800;

    // Create background bar
    const background = new Rectangle(0, 0, this.barWidth, 80, {
      fill: QPPWColors.panelFillProperty,
      stroke: null, // Remove border to avoid line appearing below chart
      lineWidth: 0,
    });
    this.addChild(background);

    // Time display
    const timeLabel = new Text(stringManager.timeStringProperty, {
      font: new PhetFont(14),
      fill: QPPWColors.textFillProperty,
    });

    this.timeText = new Text("0.00 fs", {
      font: new PhetFont({ size: 16, weight: "bold" }),
      fill: QPPWColors.textFillProperty,
    });

    this.model.timeProperty.link((time) => {
      this.timeText.string = `${time.toFixed(2)} fs`;
    });

    const timeDisplay = new VBox({
      spacing: 5,
      align: "center",
      children: [timeLabel, this.timeText],
    });

    // Time controls (play/pause and step buttons)
    const timeControlNode = new TimeControlNode(this.model.isPlayingProperty, {
      timeSpeedProperty: this.model.timeSpeedProperty,
      playPauseStepButtonOptions: {
        includeStepForwardButton: true,
        includeStepBackwardButton: true,
        stepForwardButtonOptions: {
          listener: () => {
            // Step forward by one frame (forced even when paused)
            this.model.step(BaseModel.MANUAL_STEP_SIZE, true);
          },
          enabledProperty: DerivedProperty.not(this.model.isPlayingProperty),
          radius: 15, // Smaller than play/pause button
        },
        stepBackwardButtonOptions: {
          listener: () => {
            // Step backward by one frame (negative time step, forced even when paused)
            this.model.step(-BaseModel.MANUAL_STEP_SIZE, true);
          },
          enabledProperty: DerivedProperty.not(this.model.isPlayingProperty),
          radius: 15, // Smaller than play/pause button
        },
      },
      speedRadioButtonGroupPlacement: "left",
      speedRadioButtonGroupOptions: {
        labelOptions: {
          fill: QPPWColors.textFillProperty,
        },
      },
    })


    const playbackButtons = new HBox({
      spacing: 10,
      children: [timeControlNode],
    });

    const playbackSection = new VBox({
      spacing: 8,
      align: "center",
      children: [ playbackButtons],
    });

    // Arrange all sections horizontally
    const content = new HBox({
      spacing: 40,
      align: "center",
      children: [timeDisplay, playbackSection],
      centerX: this.barWidth / 2,
      centerY: 40,
    });

    this.addChild(content);
  }
}
