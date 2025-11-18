/**
 * OneWellScreenIcon provides the icon for the One Well screen.
 * Features a colorful single potential well with a wave function visualization.
 */

import { Rectangle, Node, Path, LinearGradient } from "scenerystack/scenery";
import { Shape } from "scenerystack/kite";
import { ScreenIcon } from "scenerystack/sim";
import QPPWColors from "../../QPPWColors.js";

// Dimensions
const ICON_WIDTH = 60;
const ICON_HEIGHT = 50;
const CORNER_RADIUS = 4;

// Colors
const BACKGROUND_GRADIENT_TOP = '#1a1a3a';
const BACKGROUND_GRADIENT_BOTTOM = '#0a0a1f';
const WELL_STROKE_COLOR = '#9696c8';
const WAVE_FUNCTION_COLOR = '#00c8ff';
const PROBABILITY_FILL_COLOR = 'rgba(255, 200, 0, 0.4)';
const ENERGY_LEVEL_COLOR = '#00ff96';

// Line widths
const WELL_LINE_WIDTH = 3;
const WAVE_LINE_WIDTH = 2.5;

// Opacity
const ENERGY_LEVEL_OPACITY = 0.8;

export class OneWellScreenIcon extends ScreenIcon {
  public constructor() {
    // Create background with gradient
    const backgroundGradient = new LinearGradient(0, 0, 0, ICON_HEIGHT)
      .addColorStop(0, BACKGROUND_GRADIENT_TOP)
      .addColorStop(1, BACKGROUND_GRADIENT_BOTTOM);

    const background = new Rectangle(0, 0, ICON_WIDTH, ICON_HEIGHT, {
      fill: backgroundGradient,
      cornerRadius: CORNER_RADIUS,
    });

    // Create potential well shape (U-shaped)
    const wellShape = new Shape()
      .moveTo(5, 10)
      .lineTo(5, 40)
      .lineTo(55, 40)
      .lineTo(55, 10);

    const well = new Path(wellShape, {
      stroke: WELL_STROKE_COLOR,
      lineWidth: WELL_LINE_WIDTH,
    });

    // Create wave function (sinusoidal curve)
    const waveShape = new Shape().moveTo(8, 25);
    for (let x = 8; x <= 52; x += 1) {
      const normalizedX = (x - 8) / 44;
      const amplitude = 10 * Math.sin(normalizedX * Math.PI);
      const y = 25 - amplitude * Math.sin(normalizedX * Math.PI);
      waveShape.lineTo(x, y);
    }

    const waveFunction = new Path(waveShape, {
      stroke: WAVE_FUNCTION_COLOR,
      lineWidth: WAVE_LINE_WIDTH,
    });

    // Create probability density fill
    const probabilityShape = new Shape().moveTo(8, 38);
    for (let x = 8; x <= 52; x += 1) {
      const normalizedX = (x - 8) / 44;
      const amplitude = Math.pow(Math.sin(normalizedX * Math.PI), 2) * 18;
      probabilityShape.lineTo(x, 38 - amplitude);
    }
    probabilityShape.lineTo(52, 38).close();

    const probabilityFill = new Path(probabilityShape, {
      fill: PROBABILITY_FILL_COLOR,
    });

    // Energy level indicator
    const energyLevel = new Rectangle(10, 20, 40, 2, {
      fill: ENERGY_LEVEL_COLOR,
      opacity: ENERGY_LEVEL_OPACITY,
    });

    const iconNode = new Node({
      children: [background, probabilityFill, well, energyLevel, waveFunction],
    });

    super(iconNode, {
      fill: QPPWColors.backgroundColorProperty,
    });
  }
}
