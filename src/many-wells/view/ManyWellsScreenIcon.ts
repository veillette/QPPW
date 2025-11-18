/**
 * ManyWellsScreenIcon provides the icon for the Many Wells screen.
 * Features a colorful periodic potential with energy band visualization.
 */

import { Rectangle, Node, Path, LinearGradient } from "scenerystack/scenery";
import { Shape } from "scenerystack/kite";
import { ScreenIcon } from "scenerystack/sim";
import QPPWColors from "../../QPPWColors.js";

// Dimensions
const ICON_WIDTH = 80;
const ICON_HEIGHT = 50;
const CORNER_RADIUS = 4;
const NUM_WELLS = 4;
const WELL_WIDTH = 12;
const BARRIER_WIDTH = 6;

// Colors
const BACKGROUND_GRADIENT_TOP = '#1a1a3a';
const BACKGROUND_GRADIENT_BOTTOM = '#0a0a1f';
const POTENTIAL_STROKE_COLOR = '#9696c8';
const BARRIER_GRADIENT_EDGE = '#b43232';
const BARRIER_GRADIENT_CENTER = '#ff6b3d';
const BLOCH_WAVE_COLOR = '#00c8ff';

// Energy band colors
const BAND_COLORS = ['#00ff96', '#ffff64', '#00c8ff', '#ff96ff'];
const BAND_Y_POSITIONS = [15, 22, 29, 36];

// Line widths
const POTENTIAL_LINE_WIDTH = 2;
const BLOCH_WAVE_LINE_WIDTH = 1.5;

// Opacity
const BARRIER_OPACITY = 0.6;
const BASE_BAND_OPACITY = 0.8;
const BAND_OPACITY_DECREMENT = 0.1;

export class ManyWellsScreenIcon extends ScreenIcon {
  public constructor() {
    // Create background with gradient
    const backgroundGradient = new LinearGradient(0, 0, 0, ICON_HEIGHT)
      .addColorStop(0, BACKGROUND_GRADIENT_TOP)
      .addColorStop(1, BACKGROUND_GRADIENT_BOTTOM);

    const background = new Rectangle(0, 0, ICON_WIDTH, ICON_HEIGHT, {
      fill: backgroundGradient,
      cornerRadius: CORNER_RADIUS,
    });

    // Create periodic potential shape
    const totalWidth = NUM_WELLS * WELL_WIDTH + (NUM_WELLS - 1) * BARRIER_WIDTH;
    const startX = (ICON_WIDTH - totalWidth) / 2;

    const potentialShape = new Shape().moveTo(startX, 8);

    for (let i = 0; i < NUM_WELLS; i++) {
      const wellX = startX + i * (WELL_WIDTH + BARRIER_WIDTH);

      // Down into well
      potentialShape.lineTo(wellX, 8);
      potentialShape.lineTo(wellX, 40);

      // Across bottom of well
      potentialShape.lineTo(wellX + WELL_WIDTH, 40);

      // Up out of well (if not last)
      if (i < NUM_WELLS - 1) {
        potentialShape.lineTo(wellX + WELL_WIDTH, 8);
        // Across top of barrier
        potentialShape.lineTo(wellX + WELL_WIDTH + BARRIER_WIDTH, 8);
      } else {
        potentialShape.lineTo(wellX + WELL_WIDTH, 8);
      }
    }

    const potential = new Path(potentialShape, {
      stroke: POTENTIAL_STROKE_COLOR,
      lineWidth: POTENTIAL_LINE_WIDTH,
    });

    // Create barriers with gradient fills
    const barriers: Node[] = [];
    for (let i = 0; i < NUM_WELLS - 1; i++) {
      const barrierX = startX + (i + 1) * WELL_WIDTH + i * BARRIER_WIDTH;

      const barrierGradient = new LinearGradient(barrierX, 0, barrierX + BARRIER_WIDTH, 0)
        .addColorStop(0, BARRIER_GRADIENT_EDGE)
        .addColorStop(0.5, BARRIER_GRADIENT_CENTER)
        .addColorStop(1, BARRIER_GRADIENT_EDGE);

      barriers.push(new Rectangle(barrierX, 8, BARRIER_WIDTH, 32, {
        fill: barrierGradient,
        opacity: BARRIER_OPACITY,
      }));
    }

    // Energy bands (multiple horizontal lines with gradient colors)
    const energyBands: Node[] = [];

    for (let band = 0; band < BAND_COLORS.length; band++) {
      // Draw energy band across all wells
      for (let i = 0; i < NUM_WELLS; i++) {
        const wellX = startX + i * (WELL_WIDTH + BARRIER_WIDTH);
        energyBands.push(new Rectangle(wellX + 1, BAND_Y_POSITIONS[band], WELL_WIDTH - 2, 1.5, {
          fill: BAND_COLORS[band],
          opacity: BASE_BAND_OPACITY - band * BAND_OPACITY_DECREMENT,
        }));
      }
    }

    // Bloch wave function (sinusoidal across all wells)
    const waveShape = new Shape();
    let firstPoint = true;

    for (let i = 0; i < NUM_WELLS; i++) {
      const wellX = startX + i * (WELL_WIDTH + BARRIER_WIDTH);

      for (let x = wellX + 1; x <= wellX + WELL_WIDTH - 1; x += 0.5) {
        const localX = (x - wellX) / WELL_WIDTH;
        const globalPhase = i * Math.PI;
        const amplitude = 4 * Math.sin(localX * Math.PI);
        const y = 26 - amplitude * Math.cos(globalPhase);

        if (firstPoint) {
          waveShape.moveTo(x, y);
          firstPoint = false;
        } else {
          waveShape.lineTo(x, y);
        }
      }
    }

    const blochWave = new Path(waveShape, {
      stroke: BLOCH_WAVE_COLOR,
      lineWidth: BLOCH_WAVE_LINE_WIDTH,
    });

    const iconNode = new Node({
      children: [background, ...barriers, potential, ...energyBands, blochWave],
    });

    super(iconNode, {
      fill: QPPWColors.backgroundColorProperty,
    });
  }
}
