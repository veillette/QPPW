/**
 * OneWellScreenIcon provides the icon for the One Well screen.
 * Features a colorful single potential well with a wave function visualization.
 */

import { Rectangle, Node, Path, LinearGradient } from "scenerystack/scenery";
import { Shape } from "scenerystack/kite";
import { ScreenIcon } from "scenerystack/sim";
import QPPWColors from "../../QPPWColors.js";

export class OneWellScreenIcon extends ScreenIcon {
  public constructor() {
    const width = 60;
    const height = 50;

    // Create background with gradient
    const backgroundGradient = new LinearGradient(0, 0, 0, height)
      .addColorStop(0, '#1a1a3a')
      .addColorStop(1, '#0a0a1f');

    const background = new Rectangle(0, 0, width, height, {
      fill: backgroundGradient,
      cornerRadius: 4,
    });

    // Create potential well shape (U-shaped)
    const wellShape = new Shape()
      .moveTo(5, 10)
      .lineTo(5, 40)
      .lineTo(55, 40)
      .lineTo(55, 10);

    const well = new Path(wellShape, {
      stroke: '#9696c8',  // Light purple for well outline
      lineWidth: 3,
    });

    // Create wave function (sinusoidal curve) - bright blue
    const waveShape = new Shape().moveTo(8, 25);
    for (let x = 8; x <= 52; x += 1) {
      const normalizedX = (x - 8) / 44;
      const amplitude = 10 * Math.sin(normalizedX * Math.PI);
      const y = 25 - amplitude * Math.sin(normalizedX * Math.PI);
      waveShape.lineTo(x, y);
    }

    const waveFunction = new Path(waveShape, {
      stroke: '#00c8ff',  // Bright cyan/blue
      lineWidth: 2.5,
    });

    // Create probability density fill (semi-transparent gold)
    const probabilityShape = new Shape().moveTo(8, 38);
    for (let x = 8; x <= 52; x += 1) {
      const normalizedX = (x - 8) / 44;
      const amplitude = Math.pow(Math.sin(normalizedX * Math.PI), 2) * 18;
      probabilityShape.lineTo(x, 38 - amplitude);
    }
    probabilityShape.lineTo(52, 38).close();

    const probabilityFill = new Path(probabilityShape, {
      fill: 'rgba(255, 200, 0, 0.4)',  // Semi-transparent gold
    });

    // Energy level indicator (bright green horizontal line)
    const energyLevel = new Rectangle(10, 20, 40, 2, {
      fill: '#00ff96',  // Bright green
      opacity: 0.8,
    });

    const iconNode = new Node({
      children: [background, probabilityFill, well, energyLevel, waveFunction],
    });

    super(iconNode, {
      fill: QPPWColors.backgroundColorProperty,
    });
  }
}
