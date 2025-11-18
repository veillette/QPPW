/**
 * ManyWellsScreenIcon provides the icon for the Many Wells screen.
 * Features a colorful periodic potential with energy band visualization.
 */

import { Rectangle, Node, Path, LinearGradient } from "scenerystack/scenery";
import { Shape } from "scenerystack/kite";
import { ScreenIcon } from "scenerystack/sim";
import QPPWColors from "../../QPPWColors.js";

export class ManyWellsScreenIcon extends ScreenIcon {
  public constructor() {
    const width = 80;
    const height = 50;
    const numWells = 4;

    // Create background with gradient
    const backgroundGradient = new LinearGradient(0, 0, 0, height)
      .addColorStop(0, '#1a1a3a')
      .addColorStop(1, '#0a0a1f');

    const background = new Rectangle(0, 0, width, height, {
      fill: backgroundGradient,
      cornerRadius: 4,
    });

    // Create periodic potential shape
    const wellWidth = 12;
    const barrierWidth = 6;
    const totalWidth = numWells * wellWidth + (numWells - 1) * barrierWidth;
    const startX = (width - totalWidth) / 2;

    const potentialShape = new Shape().moveTo(startX, 8);

    for (let i = 0; i < numWells; i++) {
      const wellX = startX + i * (wellWidth + barrierWidth);

      // Down into well
      potentialShape.lineTo(wellX, 8);
      potentialShape.lineTo(wellX, 40);

      // Across bottom of well
      potentialShape.lineTo(wellX + wellWidth, 40);

      // Up out of well (if not last)
      if (i < numWells - 1) {
        potentialShape.lineTo(wellX + wellWidth, 8);
        // Across top of barrier
        potentialShape.lineTo(wellX + wellWidth + barrierWidth, 8);
      } else {
        potentialShape.lineTo(wellX + wellWidth, 8);
      }
    }

    const potential = new Path(potentialShape, {
      stroke: '#9696c8',  // Light purple
      lineWidth: 2,
    });

    // Create barriers with gradient fills
    const barriers: Node[] = [];
    for (let i = 0; i < numWells - 1; i++) {
      const barrierX = startX + (i + 1) * wellWidth + i * barrierWidth;

      const barrierGradient = new LinearGradient(barrierX, 0, barrierX + barrierWidth, 0)
        .addColorStop(0, '#b43232')
        .addColorStop(0.5, '#ff6b3d')
        .addColorStop(1, '#b43232');

      barriers.push(new Rectangle(barrierX, 8, barrierWidth, 32, {
        fill: barrierGradient,
        opacity: 0.6,
      }));
    }

    // Energy bands (multiple horizontal lines with gradient colors)
    const bandColors = ['#00ff96', '#ffff64', '#00c8ff', '#ff96ff'];
    const bandYPositions = [15, 22, 29, 36];
    const energyBands: Node[] = [];

    for (let band = 0; band < bandColors.length; band++) {
      // Draw energy band across all wells
      for (let i = 0; i < numWells; i++) {
        const wellX = startX + i * (wellWidth + barrierWidth);
        energyBands.push(new Rectangle(wellX + 1, bandYPositions[band], wellWidth - 2, 1.5, {
          fill: bandColors[band],
          opacity: 0.8 - band * 0.1,
        }));
      }
    }

    // Bloch wave function (sinusoidal across all wells) - bright blue
    const waveShape = new Shape();
    let firstPoint = true;

    for (let i = 0; i < numWells; i++) {
      const wellX = startX + i * (wellWidth + barrierWidth);

      for (let x = wellX + 1; x <= wellX + wellWidth - 1; x += 0.5) {
        const localX = (x - wellX) / wellWidth;
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
      stroke: '#00c8ff',  // Bright cyan
      lineWidth: 1.5,
    });

    const iconNode = new Node({
      children: [background, ...barriers, potential, ...energyBands, blochWave],
    });

    super(iconNode, {
      fill: QPPWColors.backgroundColorProperty,
    });
  }
}
