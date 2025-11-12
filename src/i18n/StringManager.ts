/**
 * StringManager handles internationalization for the QPPW simulation.
 * It provides StringProperty instances for all translatable strings in the application.
 */

import { StringProperty } from "scenerystack/axon";

export class StringManager {
  // Simulation title
  public readonly titleStringProperty: StringProperty;

  // Screen names
  public readonly oneWellStringProperty: StringProperty;
  public readonly twoWellsStringProperty: StringProperty;
  public readonly manyWellsStringProperty: StringProperty;

  // Common labels
  public readonly energyStringProperty: StringProperty;
  public readonly positionStringProperty: StringProperty;
  public readonly wavefunctionStringProperty: StringProperty;
  public readonly probabilityStringProperty: StringProperty;
  public readonly potentialStringProperty: StringProperty;

  // Controls
  public readonly playStringProperty: StringProperty;
  public readonly pauseStringProperty: StringProperty;
  public readonly resetStringProperty: StringProperty;
  public readonly stepStringProperty: StringProperty;

  // One Well screen strings
  public readonly singleWellStringProperty: StringProperty;
  public readonly wellWidthStringProperty: StringProperty;
  public readonly wellDepthStringProperty: StringProperty;

  // Two Wells screen strings
  public readonly doubleWellStringProperty: StringProperty;
  public readonly barrierHeightStringProperty: StringProperty;
  public readonly barrierWidthStringProperty: StringProperty;
  public readonly tunnelingStringProperty: StringProperty;

  // Many Wells screen strings
  public readonly multipleWellsStringProperty: StringProperty;
  public readonly numberOfWellsStringProperty: StringProperty;
  public readonly latticeConstantStringProperty: StringProperty;
  public readonly energyBandsStringProperty: StringProperty;

  // Preferences strings
  public readonly preferencesStringProperty: StringProperty;
  public readonly numericalMethodStringProperty: StringProperty;
  public readonly numericalMethodDescriptionStringProperty: StringProperty;
  public readonly numerovStringProperty: StringProperty;
  public readonly numerovDescriptionStringProperty: StringProperty;
  public readonly dvrStringProperty: StringProperty;
  public readonly dvrDescriptionStringProperty: StringProperty;
  public readonly autoPauseWhenTabHiddenStringProperty: StringProperty;
  public readonly autoPauseDescriptionStringProperty: StringProperty;

  public constructor() {
    // Initialize all string properties with English defaults
    this.titleStringProperty = new StringProperty("Quantum Physics: Potential Wells");

    // Screen names
    this.oneWellStringProperty = new StringProperty("One Well");
    this.twoWellsStringProperty = new StringProperty("Two Wells");
    this.manyWellsStringProperty = new StringProperty("Many Wells");

    // Common labels
    this.energyStringProperty = new StringProperty("Energy");
    this.positionStringProperty = new StringProperty("Position");
    this.wavefunctionStringProperty = new StringProperty("Wavefunction");
    this.probabilityStringProperty = new StringProperty("Probability");
    this.potentialStringProperty = new StringProperty("Potential");

    // Controls
    this.playStringProperty = new StringProperty("Play");
    this.pauseStringProperty = new StringProperty("Pause");
    this.resetStringProperty = new StringProperty("Reset");
    this.stepStringProperty = new StringProperty("Step");

    // One Well screen strings
    this.singleWellStringProperty = new StringProperty("Single Well");
    this.wellWidthStringProperty = new StringProperty("Well Width");
    this.wellDepthStringProperty = new StringProperty("Well Depth");

    // Two Wells screen strings
    this.doubleWellStringProperty = new StringProperty("Double Well");
    this.barrierHeightStringProperty = new StringProperty("Barrier Height");
    this.barrierWidthStringProperty = new StringProperty("Barrier Width");
    this.tunnelingStringProperty = new StringProperty("Tunneling");

    // Many Wells screen strings
    this.multipleWellsStringProperty = new StringProperty("Multiple Wells");
    this.numberOfWellsStringProperty = new StringProperty("Number of Wells");
    this.latticeConstantStringProperty = new StringProperty("Lattice Constant");
    this.energyBandsStringProperty = new StringProperty("Energy Bands");

    // Preferences strings
    this.preferencesStringProperty = new StringProperty("Preferences");
    this.numericalMethodStringProperty = new StringProperty("Numerical Method");
    this.numericalMethodDescriptionStringProperty = new StringProperty("Choose the numerical method for solving the Schrödinger equation.");
    this.numerovStringProperty = new StringProperty("Numerov Method");
    this.numerovDescriptionStringProperty = new StringProperty("Traditional shooting method - accurate and stable for most potentials");
    this.dvrStringProperty = new StringProperty("DVR (Discrete Variable Representation)");
    this.dvrDescriptionStringProperty = new StringProperty("Matrix diagonalization method - faster and more robust for complex potentials");
    this.autoPauseWhenTabHiddenStringProperty = new StringProperty("Auto-pause when tab is hidden");
    this.autoPauseDescriptionStringProperty = new StringProperty("Automatically pause the simulation when the browser tab is not visible");
  }

  /**
   * Gets the title string property for the simulation.
   */
  public getTitleStringProperty(): StringProperty {
    return this.titleStringProperty;
  }

  /**
   * Gets an object containing the screen name properties.
   */
  public getScreenNames() {
    return {
      oneWellStringProperty: this.oneWellStringProperty,
      twoWellsStringProperty: this.twoWellsStringProperty,
      manyWellsStringProperty: this.manyWellsStringProperty,
    };
  }

  /**
   * Gets an object containing the preferences string properties.
   */
  public getPreferencesLabels() {
    return {
      preferencesStringProperty: this.preferencesStringProperty,
      numericalMethodStringProperty: this.numericalMethodStringProperty,
      numericalMethodDescriptionStringProperty: this.numericalMethodDescriptionStringProperty,
      autoPauseWhenTabHiddenStringProperty: this.autoPauseWhenTabHiddenStringProperty,
      autoPauseDescriptionStringProperty: this.autoPauseDescriptionStringProperty,
    };
  }

  /**
   * Gets an object containing the numerical method name properties.
   */
  public getNumericalMethodNames() {
    return {
      numerovStringProperty: this.numerovStringProperty,
      dvrStringProperty: this.dvrStringProperty,
    };
  }

  /**
   * Gets an object containing the numerical method description properties.
   */
  public getNumericalMethodDescriptions() {
    return {
      numerovStringProperty: this.numerovDescriptionStringProperty,
      dvrStringProperty: this.dvrDescriptionStringProperty,
    };
  }
}

// Create and export a singleton instance
const stringManager = new StringManager();
export default stringManager;
