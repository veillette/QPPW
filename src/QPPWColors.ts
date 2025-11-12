/**
 * Color profiles for the Quantum Bound States (QPPW) Simulation
 * Provides support for different color profiles (default and projector mode)
 */

import { ProfileColorProperty } from "scenerystack/scenery";
import qppw from "./QPPWNamespace.js";

const QPPWColors = {
  // Backgrounds
  screenBackgroundColorProperty: new ProfileColorProperty(qppw, "screenBackground", {
    default: "#000000",
    projector: "#ffffff",
  }),

  graphBackgroundColorProperty: new ProfileColorProperty(qppw, "graphBackground", {
    default: "#111111",
    projector: "#f5f5f5",
  }),

  // Control Panels
  controlPanelBackgroundColorProperty: new ProfileColorProperty(qppw, "controlPanelBackground", {
    default: "rgba(40, 40, 40, 0.9)",
    projector: "rgba(245, 245, 245, 0.9)",
  }),

  controlPanelStrokeColorProperty: new ProfileColorProperty(qppw, "controlPanelStroke", {
    default: "#666666",
    projector: "#999999",
  }),

  // Text
  textColorProperty: new ProfileColorProperty(qppw, "text", {
    default: "#ffffff",
    projector: "#000000",
  }),

  labelColorProperty: new ProfileColorProperty(qppw, "label", {
    default: "#e0e0e0",
    projector: "#333333",
  }),

  // Graph Elements
  gridLineColorProperty: new ProfileColorProperty(qppw, "gridLine", {
    default: "#444444",
    projector: "#cccccc",
  }),

  axisColorProperty: new ProfileColorProperty(qppw, "axis", {
    default: "#888888",
    projector: "#555555",
  }),

  // Wave Function Colors
  realPartColorProperty: new ProfileColorProperty(qppw, "realPart", {
    default: "#ff4444",
    projector: "#cc0000",
  }),

  imaginaryPartColorProperty: new ProfileColorProperty(qppw, "imaginaryPart", {
    default: "#4444ff",
    projector: "#0000cc",
  }),

  probabilityDensityColorProperty: new ProfileColorProperty(qppw, "probabilityDensity", {
    default: "#ffaa00",
    projector: "#dd8800",
  }),

  magnitudeColorProperty: new ProfileColorProperty(qppw, "magnitude", {
    default: "#00ff88",
    projector: "#00cc66",
  }),

  phaseColorProperty: new ProfileColorProperty(qppw, "phase", {
    default: "#ff00ff",
    projector: "#cc00cc",
  }),

  // Potential Well
  potentialWellColorProperty: new ProfileColorProperty(qppw, "potentialWell", {
    default: "#888888",
    projector: "#666666",
  }),

  potentialWellFillColorProperty: new ProfileColorProperty(qppw, "potentialWellFill", {
    default: "rgba(100, 100, 100, 0.3)",
    projector: "rgba(150, 150, 150, 0.3)",
  }),

  // Energy Levels
  energyLevelColorProperty: new ProfileColorProperty(qppw, "energyLevel", {
    default: "#00ffff",
    projector: "#0099cc",
  }),

  selectedEnergyLevelColorProperty: new ProfileColorProperty(qppw, "selectedEnergyLevel", {
    default: "#ffff00",
    projector: "#cccc00",
  }),

  // Interactive Elements
  hoverColorProperty: new ProfileColorProperty(qppw, "hover", {
    default: "#ffaa00",
    projector: "#ff8800",
  }),

  selectedColorProperty: new ProfileColorProperty(qppw, "selected", {
    default: "#00ff00",
    projector: "#00cc00",
  }),

  // Focus Indicators
  focusColorProperty: new ProfileColorProperty(qppw, "focus", {
    default: "#44aaff",
    projector: "#0066cc",
  }),

  highContrastFocusColorProperty: new ProfileColorProperty(qppw, "highContrastFocus", {
    default: "#ffff00",
    projector: "#000000",
  }),

  // UI Elements
  buttonBackgroundColorProperty: new ProfileColorProperty(qppw, "buttonBackground", {
    default: "#555555",
    projector: "#dddddd",
  }),

  disabledColorProperty: new ProfileColorProperty(qppw, "disabled", {
    default: "#555555",
    projector: "#aaaaaa",
  }),

  // Measurement Tools
  measuringTapeColorProperty: new ProfileColorProperty(qppw, "measuringTape", {
    default: "#ffff00",
    projector: "#cccc00",
  }),

  // Scene Grid
  sceneGridColorProperty: new ProfileColorProperty(qppw, "sceneGrid", {
    default: "rgba(255, 255, 255, 0.2)",
    projector: "rgba(0, 0, 0, 0.2)",
  }),

  // Info Elements
  infoButtonColorProperty: new ProfileColorProperty(qppw, "infoButton", {
    default: "#4488ff",
    projector: "#3366cc",
  }),
};

qppw.register("QPPWColors", QPPWColors);

export default QPPWColors;
