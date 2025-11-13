/**
 * SimulationControlBar provides playback controls and time display.
 * This is displayed at the bottom of the One Well screen.
 */

import { Node, Text, HBox, VBox, Rectangle } from "scenerystack/scenery";
import { PlayPauseButton, StepForwardButton, PhetFont } from "scenerystack/scenery-phet";
import { RectangularPushButton, VerticalAquaRadioButtonGroup } from "scenerystack/sun";
import { Path } from "scenerystack/scenery";
import { Shape } from "scenerystack/kite";
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

    // Speed control
    const speedLabel = new Text("Speed:", {
      font: new PhetFont(14),
      fill: QPPWColors.textFillProperty,
    });

    const speedItems = [
      {
        value: "normal" as const,
        createNode: () => new Text("Normal", { font: new PhetFont(12), fill: QPPWColors.textFillProperty }),
      },
      {
        value: "fast" as const,
        createNode: () => new Text("Fast", { font: new PhetFont(12), fill: QPPWColors.textFillProperty }),
      },
    ];

    const speedControl = new VerticalAquaRadioButtonGroup(this.model.simulationSpeedProperty, speedItems, {
      spacing: 6,
      radioButtonOptions: {
        radius: 6,
      },
    });

    const speedSection = new VBox({
      spacing: 8,
      align: "center",
      children: [speedLabel, speedControl],
    });

    // Playback controls
    const restartButton = new RectangularPushButton({
      content: this.createRestartIcon(),
      baseColor: "#e0e0e0",
      listener: () => {
        this.model.restart();
      },
      xMargin: 8,
      yMargin: 8,
    });

    const playPauseButton = new PlayPauseButton(this.model.isPlayingProperty, {
      radius: 20,
    });

    const stepButton = new StepForwardButton({
      listener: () => {
        // Step forward by a small time increment
        const dt = 0.016; // ~60 FPS
        this.model.step(dt);
      },
      radius: 15,
      enabledProperty: new DerivedProperty([this.model.isPlayingProperty], (isPlaying) => !isPlaying),
    });

    const playbackLabel = new Text("Playback:", {
      font: new PhetFont(14),
      fill: QPPWColors.textFillProperty,
    });

    const playbackButtons = new HBox({
      spacing: 10,
      children: [restartButton, playPauseButton, stepButton],
    });

    const playbackSection = new VBox({
      spacing: 8,
      align: "center",
      children: [playbackLabel, playbackButtons],
    });

    // Arrange all sections horizontally
    const content = new HBox({
      spacing: 40,
      align: "center",
      children: [timeDisplay, speedSection, playbackSection],
      centerX: this.barWidth / 2,
      centerY: 40,
    });

    this.addChild(content);
  }

  /**
   * Creates a restart/rewind icon.
   */
  private createRestartIcon(): Node {
    const size = 20;
    const shape = new Shape();

    // Draw a circular arrow (restart symbol)
    shape.arc(0, 0, size / 2, Math.PI * 0.7, Math.PI * 2.3, false);

    // Add arrow head
    shape.moveTo(size / 2, -size / 4);
    shape.lineTo(size / 2 + 4, -size / 4 - 6);
    shape.lineTo(size / 2 + 6, -size / 4 + 2);

    const path = new Path(shape, {
      stroke: "black",
      lineWidth: 2,
      lineCap: "round",
      lineJoin: "round",
    });

    // Add a small triangle at the center
    const triangleShape = new Shape();
    triangleShape.moveTo(-size / 3, size / 3);
    triangleShape.lineTo(-size / 3, -size / 3);
    triangleShape.lineTo(0, 0);
    triangleShape.close();

    const triangle = new Path(triangleShape, {
      fill: "black",
    });

    return new Node({
      children: [path, triangle],
    });
  }
}
