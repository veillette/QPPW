/**
 * ZerosVisualization displays markers at the zeros (nodes) of the wavefunction.
 * These are points where the wavefunction crosses zero.
 */

import { Node, Circle } from "scenerystack/scenery";
import { BooleanProperty, Property, DerivedProperty } from "scenerystack/axon";
import QuantumConstants from "../../model/QuantumConstants.js";
import QPPWColors from "../../../QPPWColors.js";
import {
  BaseVisualization,
  type BaseVisualizationOptions,
} from "./BaseVisualization.js";

export type ZerosVisualizationOptions = BaseVisualizationOptions;

export class ZerosVisualization extends BaseVisualization {
  public readonly showProperty: BooleanProperty;
  private readonly zerosPositionsProperty: Property<number[]>;
  private readonly accessibleDescription: Node;

  constructor(options: ZerosVisualizationOptions) {
    super(options);

    // Set PDOM attributes for accessibility
    this.tagName = "div";
    this.labelTagName = "h3";
    this.labelContent = "Wavefunction Zeros";

    this.showProperty = new BooleanProperty(false);
    this.zerosPositionsProperty = new Property<number[]>([]);

    // Create accessible description with live region
    this.accessibleDescription = new Node({
      tagName: "div",
      ariaRole: "status",
      pdomAttributes: [{ attribute: "aria-live", value: "polite" }],
      innerContent: new DerivedProperty(
        [this.showProperty, this.zerosPositionsProperty],
        (show, zeros) => {
          if (!show || zeros.length === 0) {
            return "";
          }
          const positions = zeros.map((z) => z.toFixed(2)).join(", ");
          return (
            `Wavefunction has ${zeros.length} node${zeros.length !== 1 ? "s" : ""} ` +
            `(zero crossing${zeros.length !== 1 ? "s" : ""}) ` +
            `at positions: ${positions} nanometers.`
          );
        },
      ),
    });
    this.addChild(this.accessibleDescription);

    // Link visibility to property
    this.showProperty.link((show: boolean) => {
      if (show) {
        this.show();
      } else {
        this.hide();
      }
    });
  }

  /**
   * Update the visualization with new wavefunction data
   */
  public update(xGrid: number[], wavefunction: number[]): void {
    // Clear existing zeros
    this.container.removeAllChildren();

    // Only show if enabled
    if (!this.showProperty.value) {
      this.hide();
      this.zerosPositionsProperty.value = [];
      return;
    }

    // Find zeros
    const zeros = this.findZeros(xGrid, wavefunction);

    // Update the accessible description with zero positions
    this.zerosPositionsProperty.value = zeros;

    // Create circles at each zero position
    zeros.forEach((zeroX) => {
      const x = this.dataToViewX(zeroX);
      const y = this.dataToViewY(0); // Zeros are at y=0

      const circle = new Circle(4, {
        fill: QPPWColors.energyLevelSelectedProperty,
        stroke: QPPWColors.backgroundColorProperty,
        lineWidth: 1.5,
        centerX: x,
        centerY: y,
      });

      this.container.addChild(circle);
    });

    if (zeros.length > 0) {
      this.show();
    } else {
      this.hide();
    }
  }

  /**
   * Finds zeros (sign changes) in the wavefunction.
   * Uses linear interpolation for accurate zero positions.
   */
  private findZeros(xGrid: number[], wavefunction: number[]): number[] {
    const zeros: number[] = [];

    for (let i = 0; i < wavefunction.length - 1; i++) {
      const y1 = wavefunction[i];
      const y2 = wavefunction[i + 1];

      // Check for sign change (zero crossing)
      if (y1 * y2 < 0) {
        // Linear interpolation to find more accurate zero position
        const x1 = xGrid[i];
        const x2 = xGrid[i + 1];
        const zeroX = x1 - (y1 * (x2 - x1)) / (y2 - y1);
        zeros.push(zeroX * QuantumConstants.M_TO_NM); // Convert to nm
      }
    }

    return zeros;
  }
}
