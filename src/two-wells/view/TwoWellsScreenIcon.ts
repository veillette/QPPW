/**
 * TwoWellsScreenIcon provides the icon for the Two Wells screen.
 * Features a colorful double potential well with tunneling visualization.
 */

import { Rectangle, Node, Path, LinearGradient } from "scenerystack/scenery";
import { Shape } from "scenerystack/kite";
import { ScreenIcon } from "scenerystack/sim";
import QPPWColors from "../../QPPWColors.js";

// Dimensions
const ICON_WIDTH = 70;
const ICON_HEIGHT = 50;
const CORNER_RADIUS = 4;

// Colors
const BACKGROUND_GRADIENT_TOP = '#1a1a3a';
const BACKGROUND_GRADIENT_BOTTOM = '#0a0a1f';
const WELL_STROKE_COLOR = '#9696c8';
const BARRIER_GRADIENT_EDGE = '#b43232';
const BARRIER_GRADIENT_CENTER = '#ff6b3d';
const LEFT_WAVE_COLOR = '#00c8ff';
const RIGHT_WAVE_COLOR = '#ff9632';
const TUNNEL_EFFECT_COLOR = '#a020f0';
const ENERGY_LEVEL_COLOR = '#00ff96';

// Line widths
const WELL_LINE_WIDTH = 2.5;
const WAVE_LINE_WIDTH = 2;
const TUNNEL_LINE_WIDTH = 1.5;

// Opacity
const BARRIER_OPACITY = 0.7;
const TUNNEL_OPACITY = 0.9;
const ENERGY_LEVEL_OPACITY = 0.8;

export class TwoWellsScreenIcon extends ScreenIcon {
  public constructor() {
    // Create background with gradient
    const backgroundGradient = new LinearGradient(0, 0, 0, ICON_HEIGHT)
      .addColorStop(0, BACKGROUND_GRADIENT_TOP)
      .addColorStop(1, BACKGROUND_GRADIENT_BOTTOM);

    const background = new Rectangle(0, 0, ICON_WIDTH, ICON_HEIGHT, {
      fill: backgroundGradient,
      cornerRadius: CORNER_RADIUS,
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
      stroke: WELL_STROKE_COLOR,
      lineWidth: WELL_LINE_WIDTH,
    });

    // Central barrier with gradient
    const barrierGradient = new LinearGradient(27, 0, 43, 0)
      .addColorStop(0, BARRIER_GRADIENT_EDGE)
      .addColorStop(0.5, BARRIER_GRADIENT_CENTER)
      .addColorStop(1, BARRIER_GRADIENT_EDGE);

    const barrier = new Rectangle(27, 15, 16, 25, {
      fill: barrierGradient,
      opacity: BARRIER_OPACITY,
    });

    // Left well wave function
    const leftWaveShape = new Shape().moveTo(7, 30);
    for (let x = 7; x <= 25; x += 0.5) {
      const normalizedX = (x - 7) / 18;
      const amplitude = 8 * Math.sin(normalizedX * Math.PI);
      leftWaveShape.lineTo(x, 30 - amplitude);
    }

    const leftWave = new Path(leftWaveShape, {
      stroke: LEFT_WAVE_COLOR,
      lineWidth: WAVE_LINE_WIDTH,
    });

    // Right well wave function (showing tunneled state)
    const rightWaveShape = new Shape().moveTo(45, 30);
    for (let x = 45; x <= 63; x += 0.5) {
      const normalizedX = (x - 45) / 18;
      const amplitude = 6 * Math.sin(normalizedX * Math.PI);
      rightWaveShape.lineTo(x, 30 - amplitude);
    }

    const rightWave = new Path(rightWaveShape, {
      stroke: RIGHT_WAVE_COLOR,
      lineWidth: WAVE_LINE_WIDTH,
    });

    // Tunneling effect - dashed line through barrier
    const tunnelShape = new Shape()
      .moveTo(25, 26)
      .lineTo(45, 26);

    const tunnelEffect = new Path(tunnelShape, {
      stroke: TUNNEL_EFFECT_COLOR,
      lineWidth: TUNNEL_LINE_WIDTH,
      lineDash: [2, 2],
      opacity: TUNNEL_OPACITY,
    });

    // Energy levels
    const leftEnergy = new Rectangle(8, 22, 17, 1.5, {
      fill: ENERGY_LEVEL_COLOR,
      opacity: ENERGY_LEVEL_OPACITY,
    });

    const rightEnergy = new Rectangle(45, 22, 17, 1.5, {
      fill: ENERGY_LEVEL_COLOR,
      opacity: ENERGY_LEVEL_OPACITY,
    });

    const iconNode = new Node({
      children: [background, barrier, well, leftEnergy, rightEnergy, leftWave, rightWave, tunnelEffect],
    });

    super(iconNode, {
      fill: QPPWColors.backgroundColorProperty,
    });
  }
}
