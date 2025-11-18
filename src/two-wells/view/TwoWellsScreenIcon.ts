/**
 * TwoWellsScreenIcon provides the icon for the Two Wells screen.
 * Features a colorful double potential well with tunneling visualization.
 */

import { Rectangle, Node, Path, LinearGradient } from "scenerystack/scenery";
import { Shape } from "scenerystack/kite";
import { ScreenIcon } from "scenerystack/sim";
import QPPWColors from "../../QPPWColors.js";

export class TwoWellsScreenIcon extends ScreenIcon {
  public constructor() {
    const width = 70;
    const height = 50;

    // Create background with gradient
    const backgroundGradient = new LinearGradient(0, 0, 0, height)
      .addColorStop(0, '#1a1a3a')
      .addColorStop(1, '#0a0a1f');

    const background = new Rectangle(0, 0, width, height, {
      fill: backgroundGradient,
      cornerRadius: 4,
    });

    // Create double well shape (W-shaped)
    const wellShape = new Shape()
      .moveTo(5, 8)
      .lineTo(5, 40)
      .lineTo(27, 40)
      .lineTo(27, 15)
      .lineTo(43, 15)
      .lineTo(43, 40)
      .lineTo(65, 40)
      .lineTo(65, 8);

    const well = new Path(wellShape, {
      stroke: '#9696c8',  // Light purple for well outline
      lineWidth: 2.5,
    });

    // Central barrier with gradient (red to orange)
    const barrierGradient = new LinearGradient(27, 0, 43, 0)
      .addColorStop(0, '#b43232')
      .addColorStop(0.5, '#ff6b3d')
      .addColorStop(1, '#b43232');

    const barrier = new Rectangle(27, 15, 16, 25, {
      fill: barrierGradient,
      opacity: 0.7,
    });

    // Left well wave function (bright blue)
    const leftWaveShape = new Shape().moveTo(7, 30);
    for (let x = 7; x <= 25; x += 0.5) {
      const normalizedX = (x - 7) / 18;
      const amplitude = 8 * Math.sin(normalizedX * Math.PI);
      leftWaveShape.lineTo(x, 30 - amplitude);
    }

    const leftWave = new Path(leftWaveShape, {
      stroke: '#00c8ff',  // Bright cyan/blue
      lineWidth: 2,
    });

    // Right well wave function (bright orange - showing tunneled state)
    const rightWaveShape = new Shape().moveTo(45, 30);
    for (let x = 45; x <= 63; x += 0.5) {
      const normalizedX = (x - 45) / 18;
      const amplitude = 6 * Math.sin(normalizedX * Math.PI);  // Slightly smaller
      rightWaveShape.lineTo(x, 30 - amplitude);
    }

    const rightWave = new Path(rightWaveShape, {
      stroke: '#ff9632',  // Bright orange
      lineWidth: 2,
    });

    // Tunneling effect - dashed line through barrier (purple)
    const tunnelShape = new Shape()
      .moveTo(25, 26)
      .lineTo(45, 26);

    const tunnelEffect = new Path(tunnelShape, {
      stroke: '#a020f0',  // Purple
      lineWidth: 1.5,
      lineDash: [2, 2],
      opacity: 0.9,
    });

    // Energy levels (bright green)
    const leftEnergy = new Rectangle(8, 22, 17, 1.5, {
      fill: '#00ff96',  // Bright green
      opacity: 0.8,
    });

    const rightEnergy = new Rectangle(45, 22, 17, 1.5, {
      fill: '#00ff96',  // Bright green
      opacity: 0.8,
    });

    const iconNode = new Node({
      children: [background, barrier, well, leftEnergy, rightEnergy, leftWave, rightWave, tunnelEffect],
    });

    super(iconNode, {
      fill: QPPWColors.backgroundColorProperty,
    });
  }
}
