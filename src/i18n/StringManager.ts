/**
 * StringManager - Centralized string management for the Quantum Bound States (QPPW) Simulation
 * Provides localized strings for the application
 */

import { DerivedProperty, LocalizedStringProperty, TReadOnlyProperty } from "scenerystack/axon";
import { localeProperty } from "scenerystack/joist";
import stringsEn from "./strings_en.json";
import stringsFr from "./strings_fr.json";

/**
 * StringManager provides a singleton pattern for accessing localized strings
 */
export class StringManager {
  private static instance: StringManager;

  private readonly stringsMap: Map<string, Record<string, string>>;

  private constructor() {
    this.stringsMap = new Map();
    this.stringsMap.set("en", stringsEn);
    this.stringsMap.set("fr", stringsFr);
  }

  /**
   * Get the singleton instance of StringManager
   */
  public static getInstance(): StringManager {
    if (!StringManager.instance) {
      StringManager.instance = new StringManager();
    }
    return StringManager.instance;
  }

  /**
   * Create a localized string property for a given key
   */
  private createStringProperty(key: string): TReadOnlyProperty<string> {
    return new DerivedProperty([localeProperty], (locale) => {
      const strings = this.stringsMap.get(locale) || this.stringsMap.get("en")!;
      return strings[key] || key;
    });
  }

  /**
   * Get screen name strings
   */
  public getScreenNames() {
    return {
      oneWellScreen: this.createStringProperty("oneWellScreen"),
      twoWellsScreen: this.createStringProperty("twoWellsScreen"),
      manyWellsScreen: this.createStringProperty("manyWellsScreen"),
    };
  }

  /**
   * Get control label strings
   */
  public getControlLabels() {
    return {
      wellDepth: this.createStringProperty("wellDepth"),
      wellWidth: this.createStringProperty("wellWidth"),
      wellSeparation: this.createStringProperty("wellSeparation"),
      numberOfWells: this.createStringProperty("numberOfWells"),
      barrierHeight: this.createStringProperty("barrierHeight"),
      barrierWidth: this.createStringProperty("barrierWidth"),
      mass: this.createStringProperty("mass"),
      energy: this.createStringProperty("energy"),
      eigenstate: this.createStringProperty("eigenstate"),
    };
  }

  /**
   * Get unit strings
   */
  public getUnits() {
    return {
      electronVolts: this.createStringProperty("electronVolts"),
      nanometers: this.createStringProperty("nanometers"),
      electronMass: this.createStringProperty("electronMass"),
    };
  }

  /**
   * Get graph property strings
   */
  public getGraphProperties() {
    return {
      position: this.createStringProperty("position"),
      waveFunction: this.createStringProperty("waveFunction"),
      probabilityDensity: this.createStringProperty("probabilityDensity"),
      realPart: this.createStringProperty("realPart"),
      imaginaryPart: this.createStringProperty("imaginaryPart"),
      magnitude: this.createStringProperty("magnitude"),
      phase: this.createStringProperty("phase"),
      potentialEnergy: this.createStringProperty("potentialEnergy"),
      energyLevel: this.createStringProperty("energyLevel"),
    };
  }

  /**
   * Get visualization label strings
   */
  public getVisualizationLabels() {
    return {
      showWaveFunction: this.createStringProperty("showWaveFunction"),
      showProbabilityDensity: this.createStringProperty("showProbabilityDensity"),
      showMagnitude: this.createStringProperty("showMagnitude"),
      showPhase: this.createStringProperty("showPhase"),
      showGrid: this.createStringProperty("showGrid"),
      showEnergyLevels: this.createStringProperty("showEnergyLevels"),
      showMeasuringTape: this.createStringProperty("showMeasuringTape"),
    };
  }

  /**
   * Get control button strings
   */
  public getControlButtons() {
    return {
      play: this.createStringProperty("play"),
      pause: this.createStringProperty("pause"),
      step: this.createStringProperty("step"),
      reset: this.createStringProperty("reset"),
      resetAll: this.createStringProperty("resetAll"),
    };
  }

  /**
   * Get preference strings
   */
  public getPreferences() {
    return {
      autoPauseWhenTabHidden: this.createStringProperty("autoPauseWhenTabHidden"),
      reducedMotion: this.createStringProperty("reducedMotion"),
      highContrastMode: this.createStringProperty("highContrastMode"),
      announceParameterChanges: this.createStringProperty("announceParameterChanges"),
      announceStateChanges: this.createStringProperty("announceStateChanges"),
      announceDragInteractions: this.createStringProperty("announceDragInteractions"),
    };
  }

  /**
   * Get accessibility announcement strings
   */
  public getAnnouncements() {
    return {
      simulationPaused: this.createStringProperty("simulationPaused"),
      simulationPlaying: this.createStringProperty("simulationPlaying"),
      simulationReset: this.createStringProperty("simulationReset"),
      wellDepthChanged: this.createStringProperty("wellDepthChanged"),
      wellWidthChanged: this.createStringProperty("wellWidthChanged"),
      energyLevelSelected: this.createStringProperty("energyLevelSelected"),
    };
  }

  /**
   * Get keyboard shortcut strings
   */
  public getKeyboardShortcuts() {
    return {
      playPause: this.createStringProperty("playPause"),
      stepForward: this.createStringProperty("stepForward"),
      resetAll: this.createStringProperty("resetAll"),
    };
  }

  /**
   * Get screen summary strings (for accessibility)
   */
  public getScreenSummaries() {
    return {
      oneWellScreenSummary: this.createStringProperty("oneWellScreenSummary"),
      twoWellsScreenSummary: this.createStringProperty("twoWellsScreenSummary"),
      manyWellsScreenSummary: this.createStringProperty("manyWellsScreenSummary"),
    };
  }
}
