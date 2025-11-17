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
  public readonly potentialEnergyStringProperty: StringProperty;

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
  public readonly wellSeparationStringProperty: StringProperty;

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
  public readonly matrixNumerovStringProperty: StringProperty;
  public readonly matrixNumerovDescriptionStringProperty: StringProperty;
  public readonly dvrStringProperty: StringProperty;
  public readonly dvrDescriptionStringProperty: StringProperty;
  public readonly fghStringProperty: StringProperty;
  public readonly fghDescriptionStringProperty: StringProperty;
  public readonly spectralStringProperty: StringProperty;
  public readonly spectralDescriptionStringProperty: StringProperty;
  public readonly autoPauseWhenTabHiddenStringProperty: StringProperty;
  public readonly autoPauseDescriptionStringProperty: StringProperty;

  // Control Panel strings
  public readonly energyChartStringProperty: StringProperty;
  public readonly bottomChartStringProperty: StringProperty;
  public readonly particleMassStringProperty: StringProperty;
  public readonly wellConfigurationStringProperty: StringProperty;
  public readonly potentialWellStringProperty: StringProperty;
  public readonly displayStringProperty: StringProperty;
  public readonly waveFunctionViewsStringProperty: StringProperty;

  // Potential type strings
  public readonly squareInfiniteStringProperty: StringProperty;
  public readonly squareFiniteStringProperty: StringProperty;
  public readonly harmonicOscillatorStringProperty: StringProperty;
  public readonly morseStringProperty: StringProperty;
  public readonly poschlTellerStringProperty: StringProperty;
  public readonly rosenMorseStringProperty: StringProperty;
  public readonly eckartStringProperty: StringProperty;
  public readonly asymmetricTriangleStringProperty: StringProperty;
  public readonly coulomb1DStringProperty: StringProperty;
  public readonly coulomb3DStringProperty: StringProperty;
  public readonly doubleSquareWellStringProperty: StringProperty;

  // Display mode strings
  public readonly probabilityDensityStringProperty: StringProperty;
  public readonly phaseColorStringProperty: StringProperty;
  public readonly realPartStringProperty: StringProperty;
  public readonly imaginaryPartStringProperty: StringProperty;
  public readonly magnitudeStringProperty: StringProperty;
  public readonly phaseStringProperty: StringProperty;

  // Axis labels
  public readonly energyEvStringProperty: StringProperty;
  public readonly positionNmStringProperty: StringProperty;
  public readonly totalEnergyStringProperty: StringProperty;
  public readonly waveFunctionMagnitudeStringProperty: StringProperty;
  public readonly energyLevelLabelStringProperty: StringProperty;

  // Units
  public readonly electronVoltsStringProperty: StringProperty;
  public readonly nanometersStringProperty: StringProperty;
  public readonly electronMassStringProperty: StringProperty;

  // Time display
  public readonly timeStringProperty: StringProperty;

  // Screen descriptions
  public readonly oneWellDescriptionStringProperty: StringProperty;
  public readonly twoWellsDescriptionStringProperty: StringProperty;
  public readonly manyWellsDescriptionStringProperty: StringProperty;
  public readonly oneWellSummaryStringProperty: StringProperty;
  public readonly twoWellsSummaryStringProperty: StringProperty;
  public readonly manyWellsSummaryStringProperty: StringProperty;

  // Superposition strings
  public readonly superpositionStringProperty: StringProperty;
  public readonly psiIPsiJStringProperty: StringProperty;
  public readonly psiKStringProperty: StringProperty;
  public readonly localizedNarrowStringProperty: StringProperty;
  public readonly localizedWideStringProperty: StringProperty;
  public readonly coherentStateStringProperty: StringProperty;
  public readonly customStringProperty: StringProperty;
  public readonly configureSuperpositionStringProperty: StringProperty;
  public readonly amplitudeStringProperty: StringProperty;
  public readonly superpositionDialogTitleStringProperty: StringProperty;
  public readonly closeStringProperty: StringProperty;
  public readonly displacementStringProperty: StringProperty;

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
    this.potentialEnergyStringProperty = new StringProperty("Potential Energy");

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
    this.wellSeparationStringProperty = new StringProperty("Well Separation");

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
    this.matrixNumerovStringProperty = new StringProperty("Matrix Numerov");
    this.matrixNumerovDescriptionStringProperty = new StringProperty("Matrix diagonalization using Numerov formula - combines O(h⁴) accuracy with robustness");
    this.dvrStringProperty = new StringProperty("DVR (Discrete Variable Representation)");
    this.dvrDescriptionStringProperty = new StringProperty("Matrix diagonalization method - faster and more robust for complex potentials");
    this.fghStringProperty = new StringProperty("FGH (Fourier Grid Hamiltonian)");
    this.fghDescriptionStringProperty = new StringProperty("Plane wave basis method - natural for periodic systems with spectral accuracy");
    this.spectralStringProperty = new StringProperty("Spectral (Chebyshev)");
    this.spectralDescriptionStringProperty = new StringProperty("Chebyshev polynomial method - exponential convergence for smooth functions");
    this.autoPauseWhenTabHiddenStringProperty = new StringProperty("Auto-pause when tab is hidden");
    this.autoPauseDescriptionStringProperty = new StringProperty("Automatically pause the simulation when the browser tab is not visible");

    // Control Panel strings
    this.energyChartStringProperty = new StringProperty("Energy Chart");
    this.bottomChartStringProperty = new StringProperty("Bottom Chart");
    this.particleMassStringProperty = new StringProperty("Particle Mass");
    this.wellConfigurationStringProperty = new StringProperty("Well Configuration");
    this.potentialWellStringProperty = new StringProperty("Potential Well:");
    this.displayStringProperty = new StringProperty("Display:");
    this.waveFunctionViewsStringProperty = new StringProperty("Wave Function views:");

    // Potential type strings
    this.squareInfiniteStringProperty = new StringProperty("Square (Infinite)");
    this.squareFiniteStringProperty = new StringProperty("Square (Finite)");
    this.harmonicOscillatorStringProperty = new StringProperty("Harmonic Oscillator");
    this.morseStringProperty = new StringProperty("Morse");
    this.poschlTellerStringProperty = new StringProperty("Pöschl-Teller");
    this.rosenMorseStringProperty = new StringProperty("Rosen-Morse");
    this.eckartStringProperty = new StringProperty("Eckart");
    this.asymmetricTriangleStringProperty = new StringProperty("Asymmetric Triangle");
    this.coulomb1DStringProperty = new StringProperty("1D Coulomb");
    this.coulomb3DStringProperty = new StringProperty("3D Coulomb");
    this.doubleSquareWellStringProperty = new StringProperty("Double Square Well");

    // Display mode strings
    this.probabilityDensityStringProperty = new StringProperty("Probability Density");
    this.phaseColorStringProperty = new StringProperty("Phase (Color)");
    this.realPartStringProperty = new StringProperty("real part");
    this.imaginaryPartStringProperty = new StringProperty("imaginary part");
    this.magnitudeStringProperty = new StringProperty("magnitude");
    this.phaseStringProperty = new StringProperty("phase");

    // Axis labels
    this.energyEvStringProperty = new StringProperty("Energy (eV)");
    this.positionNmStringProperty = new StringProperty("Position (nm)");
    this.totalEnergyStringProperty = new StringProperty("Total Energy");
    this.waveFunctionMagnitudeStringProperty = new StringProperty("Wave Function Magnitude");
    this.energyLevelLabelStringProperty = new StringProperty("E{{level}} = {{value}} eV");

    // Units
    this.electronVoltsStringProperty = new StringProperty("eV");
    this.nanometersStringProperty = new StringProperty("nm");
    this.electronMassStringProperty = new StringProperty("m_e");

    // Time display
    this.timeStringProperty = new StringProperty("Time:");

    // Screen descriptions
    this.oneWellDescriptionStringProperty = new StringProperty(
      "Explore quantum mechanics in a single potential well.\nAdjust the well parameters to see how energy levels change."
    );
    this.twoWellsDescriptionStringProperty = new StringProperty(
      "Explore quantum tunneling in a double potential well.\nAdjust barrier parameters to see how tunneling probability changes.\nWatch particles tunnel through classically forbidden regions!"
    );
    this.manyWellsDescriptionStringProperty = new StringProperty(
      "Explore energy bands in a periodic potential.\nAdd or remove wells to see how energy bands form.\nThis demonstrates the foundation of solid-state physics!"
    );
    this.oneWellSummaryStringProperty = new StringProperty("One Well screen shows a single quantum potential well with adjustable parameters.");
    this.twoWellsSummaryStringProperty = new StringProperty("Two Wells screen demonstrates quantum tunneling between two potential wells.");
    this.manyWellsSummaryStringProperty = new StringProperty("Many Wells screen demonstrates energy band formation in periodic potentials.");

    // Superposition strings
    this.superpositionStringProperty = new StringProperty("Superposition:");
    this.psiIPsiJStringProperty = new StringProperty("ψₙ, ψₘ");
    this.psiKStringProperty = new StringProperty("ψₖ");
    this.localizedNarrowStringProperty = new StringProperty("Localized narrow");
    this.localizedWideStringProperty = new StringProperty("Localized wide");
    this.coherentStateStringProperty = new StringProperty("Coherent state");
    this.customStringProperty = new StringProperty("Custom...");
    this.configureSuperpositionStringProperty = new StringProperty("Configure Superposition");
    this.amplitudeStringProperty = new StringProperty("Amplitude");
    this.superpositionDialogTitleStringProperty = new StringProperty("Superposition Configuration");
    this.closeStringProperty = new StringProperty("Close");
    this.displacementStringProperty = new StringProperty("Displacement");
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
      matrixNumerovStringProperty: this.matrixNumerovStringProperty,
      dvrStringProperty: this.dvrStringProperty,
      fghStringProperty: this.fghStringProperty,
      spectralStringProperty: this.spectralStringProperty,
    };
  }

  /**
   * Gets an object containing the numerical method description properties.
   */
  public getNumericalMethodDescriptions() {
    return {
      numerovStringProperty: this.numerovDescriptionStringProperty,
      matrixNumerovStringProperty: this.matrixNumerovDescriptionStringProperty,
      dvrStringProperty: this.dvrDescriptionStringProperty,
      fghStringProperty: this.fghDescriptionStringProperty,
      spectralStringProperty: this.spectralDescriptionStringProperty,
    };
  }
}

// Create and export a singleton instance
const stringManager = new StringManager();
export default stringManager;
