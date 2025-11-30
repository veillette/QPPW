/**
 * ScreenSummaryNode provides the accessible screen summary for QPPW screens.
 * It implements the first section of the PhET Parallel DOM (PDOM) three-section layout.
 *
 * This component creates a dynamic, accessible description of the current simulation state
 * that updates automatically when model properties change. Screen readers announce this
 * content when users first navigate to a screen.
 */

import { Node } from "scenerystack/scenery";
import { DerivedProperty } from "scenerystack/axon";
import type { BaseModel } from "../../model/BaseModel.js";
import type { OneWellModel } from "../../../one-well/model/OneWellModel.js";
import type { TwoWellsModel } from "../../../two-wells/model/TwoWellsModel.js";
import type { ManyWellsModel } from "../../../many-wells/model/ManyWellsModel.js";

export type ScreenSummaryOptions = {
  screenName: string;
  screenDescription: string;
};

/**
 * ScreenSummaryNode creates an accessible screen summary following PhET PDOM patterns.
 * The summary includes:
 * - Screen title and description
 * - Current potential type and configuration
 * - Selected energy level information
 * - Key simulation parameters
 * - Statistical properties of the wavefunction
 */
export class ScreenSummaryNode extends Node {
  private readonly model: BaseModel | OneWellModel | TwoWellsModel | ManyWellsModel;
  private readonly screenName: string;

  public constructor(
    model: BaseModel | OneWellModel | TwoWellsModel | ManyWellsModel,
    options: ScreenSummaryOptions
  ) {
    super({
      tagName: "div",

      // PDOM - Accessible label for the screen summary section
      labelTagName: "h2",
      labelContent: "Screen Summary",

      // PDOM - Description of what this section contains
      descriptionTagName: "p",
      descriptionContent: `${options.screenDescription}`,
    });

    this.model = model;
    this.screenName = options.screenName;

    // Add current state description that updates dynamically
    this.addChild(this.createCurrentStateNode());

    // Add particle and potential parameters description
    this.addChild(this.createParametersNode());

    // Add statistical properties description
    this.addChild(this.createStatisticsNode());
  }

  /**
   * Creates a node describing the current simulation state.
   * This updates when potential type or energy level changes.
   */
  private createCurrentStateNode(): Node {
    return new Node({
      tagName: "p",
      innerContent: new DerivedProperty(
        [
          this.model.potentialTypeProperty,
          this.model.selectedEnergyLevelIndexProperty,
        ],
        (potentialType, levelIndex) => {
          const energyLevels = this.model.getEnergyLevels();
          if (energyLevels.length === 0) {
            return `Currently exploring a ${potentialType.name} potential well. No bound states found.`;
          }

          const energy = energyLevels[levelIndex];
          const levelNumber = levelIndex + 1; // Convert to 1-indexed for display
          const totalLevels = energyLevels.length;

          return (
            `Currently exploring a ${potentialType.name} potential well. ` +
            `Selected energy level ${levelNumber} of ${totalLevels} ` +
            `with energy ${energy.toFixed(3)} electron volts.`
          );
        }
      ),
    });
  }

  /**
   * Creates a node describing particle and potential parameters.
   * This updates when mass, width, or depth properties change.
   */
  private createParametersNode(): Node {
    // Check if model has depth property (not all potentials do)
    const hasDepth = "depthProperty" in this.model;

    const properties = hasDepth
      ? [
          this.model.massProperty,
          this.model.widthProperty,
          (this.model as OneWellModel | TwoWellsModel).depthProperty,
        ]
      : [this.model.massProperty, this.model.widthProperty];

    return new Node({
      tagName: "p",
      innerContent: new DerivedProperty(properties, (...values) => {
        const mass = values[0];
        const width = values[1];
        let description = `Particle mass: ${mass.toFixed(2)} electron masses. Well width: ${width.toFixed(2)} nanometers.`;

        if (hasDepth && values[2] !== undefined) {
          const depth = values[2];
          description += ` Well depth: ${depth.toFixed(2)} electron volts.`;
        }

        return description;
      }),
    });
  }

  /**
   * Creates a node describing wavefunction statistical properties.
   * This updates when position statistics change.
   */
  private createStatisticsNode(): Node {
    return new Node({
      tagName: "p",
      innerContent: new DerivedProperty(
        [this.model.rmsPositionProperty, this.model.averagePositionProperty],
        (rmsPosition, avgPosition) => {
          return (
            `Position statistics: Average position ${avgPosition.toFixed(2)} nanometers, ` +
            `position uncertainty (RMS) ${rmsPosition.toFixed(2)} nanometers.`
          );
        }
      ),
    });
  }

  /**
   * Creates an accessible description of the energy spectrum.
   * This can be called to get detailed information about all energy levels.
   */
  private createEnergySpectrumDescription(): string {
    const energyLevels = this.model.getEnergyLevels();
    if (energyLevels.length === 0) {
      return "No bound states exist for the current potential configuration.";
    }

    const groundEnergy = energyLevels[0];
    let description = `Energy spectrum contains ${energyLevels.length} bound state${energyLevels.length !== 1 ? "s" : ""}. `;
    description += `Ground state energy: ${groundEnergy.toFixed(3)} eV. `;

    if (energyLevels.length > 1) {
      const firstExcited = energyLevels[1];
      const spacing = firstExcited - groundEnergy;
      description += `First excited state: ${firstExcited.toFixed(3)} eV. `;
      description += `Energy spacing: ${spacing.toFixed(3)} eV.`;
    }

    return description;
  }
}
