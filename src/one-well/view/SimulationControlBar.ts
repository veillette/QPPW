/**
 * SimulationControlBar provides playback controls and time display.
 * This is displayed at the bottom of the One Well screen.
 */

import { Node, Text, HBox, VBox, Rectangle } from "scenerystack/scenery";
import {
  TimeControlNode,
  PhetFont,
  ResetAllButton,
} from "scenerystack/scenery-phet";

import { DerivedProperty } from "scenerystack/axon";
import { OneWellModel } from "../model/OneWellModel.js";
import QPPWColors from "../../QPPWColors.js";

export class SimulationControlBar extends Node {
  private readonly model: OneWellModel;
  private readonly timeText: Text;
  private readonly barWidth: number;

  public constructor(model: OneWellModel, options?: { width?: number }) {
    super();

    this.model = model;
    this.barWidth = options?.width ?? 800;

    // Create background bar
    const background = new Rectangle(0, 0, this.barWidth, 80, {
      fill: QPPWColors.panelFillProperty,
      stroke: QPPWColors.panelStrokeProperty,
      lineWidth: 1,
    });
    this.addChild(background);

    // Time display
    const timeLabel = new Text("Time:", {
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



     // Default time step for manual stepping (in seconds)
    const manualStepSize = 0.016; // ~1 frame at 60 FPS

    // Create derived property: stepper buttons enabled only when paused
    const stepperEnabledProperty = new DerivedProperty(
      [this.model.isPlayingProperty],
      (isPlaying) => !isPlaying,
    );

    // Time controls (play/pause and step buttons)
    const timeControlNode = new TimeControlNode(this.model.isPlayingProperty, {
      timeSpeedProperty: this.model.timeSpeedProperty,
      playPauseStepButtonOptions: {
        includeStepForwardButton: true,
        includeStepBackwardButton: true,
        stepForwardButtonOptions: {
          listener: () => {
            // Step forward by one frame (forced even when paused)
            this.model.step(manualStepSize, true);
          },
          enabledProperty: stepperEnabledProperty,
          radius: 15, // Smaller than play/pause button
        },
        stepBackwardButtonOptions: {
          listener: () => {
            // Step backward by one frame (negative time step, forced even when paused)
            this.model.step(-manualStepSize, true);
          },
          enabledProperty: stepperEnabledProperty,
          radius: 15, // Smaller than play/pause button
        },
      },
      speedRadioButtonGroupPlacement: "left",
      speedRadioButtonGroupOptions: {
        labelOptions: {
          fill: QPPWColors.textFillProperty,
        },
      },
    });


    // Reset button
    const resetButton = new ResetAllButton({
      listener: () => {
        this.model.reset();
      }
    });
    this.addChild(resetButton);


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
