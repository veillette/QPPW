/**
 * QPPWColors defines the color palette for the Quantum Physics Potential Wells simulation.
 * Colors are defined as ProfileColorProperty instances which can adapt to different color profiles.
 */

import { ProfileColorProperty } from "scenerystack/scenery";
import { Color } from "scenerystack/scenery";

const QPPWColors = {
  // Background colors
  backgroundColorProperty: new ProfileColorProperty("qppw", "background", {
    default: new Color(10, 10, 30), // Deep blue-black for quantum theme
  }),

  // Panel and control colors
  panelFillProperty: new ProfileColorProperty("qppw", "panelFill", {
    default: new Color(240, 240, 245),
  }),

  panelStrokeProperty: new ProfileColorProperty("qppw", "panelStroke", {
    default: new Color(100, 100, 120),
  }),

  // Wavefunction colors
  wavefunctionRealProperty: new ProfileColorProperty("qppw", "wavefunctionReal", {
    default: new Color(0, 150, 255), // Bright blue
  }),

  wavefunctionImaginaryProperty: new ProfileColorProperty("qppw", "wavefunctionImaginary", {
    default: new Color(255, 100, 0), // Orange
  }),

  wavefunctionProbabilityProperty: new ProfileColorProperty("qppw", "wavefunctionProbability", {
    default: new Color(255, 200, 0), // Gold/yellow
  }),

  // Potential well colors
  potentialWellProperty: new ProfileColorProperty("qppw", "potentialWell", {
    default: new Color(150, 150, 200), // Light purple
  }),

  potentialBarrierProperty: new ProfileColorProperty("qppw", "potentialBarrier", {
    default: new Color(180, 50, 50), // Red
  }),

  // Energy level colors
  energyLevelProperty: new ProfileColorProperty("qppw", "energyLevel", {
    default: new Color(0, 255, 150), // Bright green
  }),

  energyLevelSelectedProperty: new ProfileColorProperty("qppw", "energyLevelSelected", {
    default: new Color(255, 255, 100), // Bright yellow
  }),

  // Grid and axis colors
  gridLineProperty: new ProfileColorProperty("qppw", "gridLine", {
    default: new Color(80, 80, 100, 0.3),
  }),

  axisProperty: new ProfileColorProperty("qppw", "axis", {
    default: new Color(200, 200, 220),
  }),

  // Text colors
  textFillProperty: new ProfileColorProperty("qppw", "textFill", {
    default: new Color(255, 255, 255),
  }),

  labelFillProperty: new ProfileColorProperty("qppw", "labelFill", {
    default: new Color(200, 200, 220),
  }),
};

export default QPPWColors;
