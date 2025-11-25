/**
 * IntroScreenIcon provides the icon for the Intro screen.
 * Features a simple potential well with a wave function visualization.
 */

import { Rectangle, Node, Path, LinearGradient } from "scenerystack/scenery";
import { Shape } from "scenerystack/kite";
import { ScreenIcon } from "scenerystack/sim";
import QPPWColors from "../../QPPWColors.js";

// Dimensions
const ICON_WIDTH = 60;
const ICON_HEIGHT = 50;
const CORNER_RADIUS = 4;

// Layout proportions
const PADDING = 5;
const WELL_TOP = ICON_HEIGHT * 0.2;
const WELL_BOTTOM = ICON_HEIGHT - PADDING * 2;
const WAVE_PADDING = PADDING + 3;
const WAVE_CENTER_Y = ICON_HEIGHT / 2;
const WAVE_WIDTH = ICON_WIDTH - 2 * WAVE_PADDING;
const PROBABILITY_BASELINE = WELL_BOTTOM - 2;

// Wave parameters
const WAVE_AMPLITUDE = 8;
const PROBABILITY_AMPLITUDE = 15;

// Colors - using a welcoming intro theme
const BACKGROUND_GRADIENT_TOP = "#2a2a4a";
const BACKGROUND_GRADIENT_BOTTOM = "#1a1a2f";
const WELL_STROKE_COLOR = "#a8a8d8";
const WAVE_FUNCTION_COLOR = "#00d4ff";
const PROBABILITY_FILL_COLOR = "rgba(100, 200, 255, 0.5)";

// Line widths
const WELL_LINE_WIDTH = 3;
const WAVE_LINE_WIDTH = 2.5;

export class IntroScreenIcon extends ScreenIcon {
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
      .moveTo(PADDING, WELL_TOP)
      .lineTo(PADDING, WELL_BOTTOM)
      .lineTo(ICON_WIDTH - PADDING, WELL_BOTTOM)
      .lineTo(ICON_WIDTH - PADDING, WELL_TOP);

    const well = new Path(wellShape, {
      stroke: WELL_STROKE_COLOR,
      lineWidth: WELL_LINE_WIDTH,
    });

    // Create wave function (simple sinusoidal curve for intro)
    const waveShape = new Shape().moveTo(WAVE_PADDING, WAVE_CENTER_Y);
    for (let x = WAVE_PADDING; x <= WAVE_PADDING + WAVE_WIDTH; x += 1) {
      const normalizedX = (x - WAVE_PADDING) / WAVE_WIDTH;
      const amplitude = WAVE_AMPLITUDE * Math.sin(normalizedX * Math.PI);
      const y = WAVE_CENTER_Y - amplitude * Math.sin(normalizedX * Math.PI);
      waveShape.lineTo(x, y);
    }

    const waveFunction = new Path(waveShape, {
      stroke: WAVE_FUNCTION_COLOR,
      lineWidth: WAVE_LINE_WIDTH,
    });

    // Create probability density fill
    const probabilityShape = new Shape().moveTo(
      WAVE_PADDING,
      PROBABILITY_BASELINE,
    );
    for (let x = WAVE_PADDING; x <= WAVE_PADDING + WAVE_WIDTH; x += 1) {
      const normalizedX = (x - WAVE_PADDING) / WAVE_WIDTH;
      const amplitude =
        Math.pow(Math.sin(normalizedX * Math.PI), 2) * PROBABILITY_AMPLITUDE;
      probabilityShape.lineTo(x, PROBABILITY_BASELINE - amplitude);
    }
    probabilityShape
      .lineTo(WAVE_PADDING + WAVE_WIDTH, PROBABILITY_BASELINE)
      .close();

    const probabilityFill = new Path(probabilityShape, {
      fill: PROBABILITY_FILL_COLOR,
    });

    const iconNode = new Node({
      children: [background, probabilityFill, well, waveFunction],
    });

    super(iconNode, {
      fill: QPPWColors.backgroundColorProperty,
    });
  }
}
