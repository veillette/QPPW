/**
 * ManyWellsModel represents the physics model for multiple quantum potential wells.
 * Similar to TwoWellsModel but generalized to N wells (1-10).
 * Supports multi-square wells and multi-Coulomb 1D potentials.
 */

import { NumberProperty, Property } from "scenerystack/axon";
import { Range } from "scenerystack/dot";
import { BaseModel } from "../../common/model/BaseModel.js";
import {
  WellParameters,
  NumericalMethod,
} from "../../common/model/Schrodinger1DSolver.js";
import { PotentialType } from "../../common/model/PotentialFunction.js";
import QuantumConstants from "../../common/model/QuantumConstants.js";
import { SuperpositionType } from "../../common/model/SuperpositionType.js";
import QPPWPreferences from "../../QPPWPreferences.js";

export type DisplayMode = "probabilityDensity" | "waveFunction" | "phaseColor";

export class ManyWellsModel extends BaseModel {
  // Number of wells (1-10)
  public readonly numberOfWellsProperty: NumberProperty;

  // Model-specific well parameters
  public readonly wellSeparationProperty: NumberProperty;

  // Electric field (eV/nm) - creates a linear tilt in the potential
  public readonly electricFieldProperty: NumberProperty;

  // Model-specific display settings
  public readonly displayModeProperty: Property<DisplayMode>; // Extends base with "phaseColor" mode
  public readonly energyLevelProperty: NumberProperty; // Deprecated, kept for compatibility

  public constructor() {
    super();

    // Override potential type to default to multi-square well
    this.potentialTypeProperty.value = PotentialType.MULTI_SQUARE_WELL;

    // Override well width range for multi-square well
    this.wellWidthProperty.setValueAndRange(1.0, new Range(0.1, 3.0));

    // Initialize number of wells (default: 3 wells)
    this.numberOfWellsProperty = new NumberProperty(3, {
      range: new Range(1, 10),
    });

    // Initialize model-specific well parameters
    this.wellSeparationProperty = new NumberProperty(0.2, {
      range: new Range(0.05, 0.7),
    }); // in nanometers (spacing between wells)

    // Initialize electric field (default: 0, range: ±0.625 eV/nm for ±5eV tilt over 8nm)
    this.electricFieldProperty = new NumberProperty(0.0, {
      range: new Range(-0.625, 0.625),
    }); // in eV/nm

    // Initialize model-specific display settings
    this.displayModeProperty = new Property<DisplayMode>("probabilityDensity");
    this.energyLevelProperty = new NumberProperty(0); // Deprecated

    // Override superposition config default
    this.superpositionConfigProperty.value = {
      type: SuperpositionType.PSI_I_PSI_J,
      amplitudes: [0.7, 0.7], // Default to equal superposition
      phases: [0, 0],
    };
  }

  /**
   * Setup cache invalidation listeners, including model-specific properties.
   * Extends the base implementation to add ManyWellsModel-specific invalidation.
   */
  protected override setupCacheInvalidation(): void {
    super.setupCacheInvalidation();

    const invalidateCache = () => {
      this.boundStateResult = null;
    };

    this.numberOfWellsProperty.link(invalidateCache);
    this.wellSeparationProperty.link(invalidateCache);
    this.electricFieldProperty.link(invalidateCache);
  }

  /**
   * Called when the solver method changes.
   * Invalidates the cached bound state results.
   * @param _method - The new numerical method (unused but required by interface)
   */
  protected override onSolverMethodChanged(_method: NumericalMethod): void {
    this.boundStateResult = null; // Invalidate cache
  }

  /**
   * Resets all properties to their initial state.
   * Override from BaseModel to reset model-specific properties.
   */
  public override resetAll(): void {
    super.resetAll();
    this.numberOfWellsProperty.reset();
    this.wellSeparationProperty.reset();
    this.electricFieldProperty.reset();
    this.energyLevelProperty.reset();
    this.displayModeProperty.reset();
  }

  /**
   * Steps the model forward in time.
   * @param dt - The time step in seconds (can be negative for backward stepping)
   * @param forced - If true, steps even when paused (for manual stepping buttons)
   */
  public override step(dt: number, forced = false): void {
    super.step(dt, forced);
    // Add dynamics here if needed
  }

  /**
   * Calculate bound states using the Schrödinger solver.
   * Results are cached until well parameters change.
   * Override from BaseModel.
   */
  protected override calculateBoundStates(): void {
    const wellWidth = this.wellWidthProperty.value * QuantumConstants.NM_TO_M;
    const mass =
      this.particleMassProperty.value * QuantumConstants.ELECTRON_MASS;

    // Calculate number of states based on potential type
    const numStates = 80; // Default for multi-well potentials

    // Grid configuration
    const CHART_DISPLAY_RANGE_NM = 4;
    const method = this.solver.getNumericalMethod();
    let numGridPoints = QPPWPreferences.gridPointsProperty.value;

    if (method === "fgh") {
      // FGH: round to nearest power of 2 for FFT efficiency
      numGridPoints = Math.pow(2, Math.round(Math.log2(numGridPoints)));
    }

    const gridConfig = {
      xMin: -CHART_DISPLAY_RANGE_NM * QuantumConstants.NM_TO_M,
      xMax: CHART_DISPLAY_RANGE_NM * QuantumConstants.NM_TO_M,
      numPoints: numGridPoints,
    };

    try {
      // Build potential parameters based on type
      const potentialParams: WellParameters = {
        type: this.potentialTypeProperty.value,
        numberOfWells: this.numberOfWellsProperty.value,
        wellWidth: wellWidth,
      };

      // Add type-specific parameters
      switch (this.potentialTypeProperty.value) {
        case PotentialType.MULTI_SQUARE_WELL: {
          potentialParams.wellDepth =
            this.wellDepthProperty.value * QuantumConstants.EV_TO_JOULES;
          potentialParams.wellSeparation =
            this.wellSeparationProperty.value * QuantumConstants.NM_TO_M;
          break;
        }
        case PotentialType.MULTI_COULOMB_1D: {
          // For Coulomb potentials, use coulombStrength parameter α = k*e²
          const coulombConstant = 8.9875517923e9; // Coulomb's constant in N·m²/C²
          potentialParams.coulombStrength =
            coulombConstant *
            QuantumConstants.ELEMENTARY_CHARGE *
            QuantumConstants.ELEMENTARY_CHARGE;
          potentialParams.wellSeparation =
            this.wellSeparationProperty.value * QuantumConstants.NM_TO_M;
          break;
        }
        default:
          break;
      }

      // Solve using the analytical (or numerical) solution
      this.boundStateResult = this.solver.solveAnalyticalIfPossible(
        potentialParams,
        mass,
        numStates,
        gridConfig,
      );

      // Ensure selected energy level index is within bounds
      if (this.boundStateResult) {
        const maxIndex = this.boundStateResult.energies.length - 1;
        if (this.selectedEnergyLevelIndexProperty.value > maxIndex) {
          this.selectedEnergyLevelIndexProperty.value = Math.max(0, maxIndex);
        }
      }
    } catch (error) {
      console.error("Error calculating bound states:", error);
      this.boundStateResult = null;
    }
  }

  /**
   * Calculate the classical probability density for a given energy level.
   * Override from BaseModel to provide potential-specific implementations.
   * The classical probability density is inversely proportional to the velocity:
   * P(x) ∝ 1/v(x) = 1/√[2(E - V(x))/m]
   *
   * @param energyIndex - Index of the energy level (0-indexed)
   * @returns Array of classical probability density values, or null if unavailable
   */
  public override getClassicalProbabilityDensity(energyIndex: number): number[] | null {
    if (!this.boundStateResult) {
      this.calculateBoundStates();
    }

    if (
      !this.boundStateResult ||
      energyIndex < 0 ||
      energyIndex >= this.boundStateResult.energies.length
    ) {
      return null;
    }

    const energy = this.boundStateResult.energies[energyIndex];
    const xGrid = this.boundStateResult.xGrid;
    const mass =
      this.particleMassProperty.value * QuantumConstants.ELECTRON_MASS;

    // Calculate potential at each grid point
    const potential = this.calculatePotentialEnergy(xGrid);

    // Use BaseModel's common method to calculate classical probability density
    return this.calculateClassicalProbabilityDensity(
      potential,
      energy,
      mass,
      xGrid,
    );
  }

  /**
   * Calculate the potential energy at given positions.
   * @param xGrid - Array of x positions in meters
   * @returns Array of potential energy values in Joules
   */
  private calculatePotentialEnergy(xGrid: number[]): number[] {
    const numberOfWells = this.numberOfWellsProperty.value;
    const wellWidth = this.wellWidthProperty.value * QuantumConstants.NM_TO_M;
    const wellDepth =
      this.wellDepthProperty.value * QuantumConstants.EV_TO_JOULES;
    const wellSeparation =
      this.wellSeparationProperty.value * QuantumConstants.NM_TO_M;

    const potential: number[] = [];

    for (let i = 0; i < xGrid.length; i++) {
      const x = xGrid[i];
      let V = 0;

      switch (this.potentialTypeProperty.value) {
        case PotentialType.MULTI_SQUARE_WELL: {
          // Multiple square wells arranged periodically
          const halfWellWidth = wellWidth / 2;
          const period = wellWidth + wellSeparation;

          // Calculate total extent of the well array
          const totalExtent =
            numberOfWells * wellWidth + (numberOfWells - 1) * wellSeparation;
          const arrayStart = -totalExtent / 2;

          let inWell = false;

          // Check if x is inside any of the wells
          for (let wellIndex = 0; wellIndex < numberOfWells; wellIndex++) {
            const wellCenter = arrayStart + wellIndex * period + wellWidth / 2;
            const wellStart = wellCenter - halfWellWidth;
            const wellEnd = wellCenter + halfWellWidth;

            if (x >= wellStart && x <= wellEnd) {
              inWell = true;
              break;
            }
          }

          V = inWell ? 0 : wellDepth;
          break;
        }

        case PotentialType.MULTI_COULOMB_1D: {
          // Multiple Coulomb centers arranged periodically
          const coulombConstant = 8.9875517923e9;
          const coulombStrength =
            coulombConstant *
            QuantumConstants.ELEMENTARY_CHARGE *
            QuantumConstants.ELEMENTARY_CHARGE;

          // Calculate total extent of the Coulomb centers
          const totalExtent = (numberOfWells - 1) * wellSeparation;
          const arrayStart = -totalExtent / 2;

          // Sum contributions from all Coulomb centers
          V = 0;
          for (let wellIndex = 0; wellIndex < numberOfWells; wellIndex++) {
            const centerPosition = arrayStart + wellIndex * wellSeparation;
            const r = Math.abs(x - centerPosition);

            if (r > 1e-12) {
              // Avoid singularity at center
              V += -coulombStrength / r;
            } else {
              V += -coulombStrength / 1e-12;
            }
          }
          break;
        }

        default:
          V = 0;
          break;
      }

      potential.push(V);
    }

    return potential;
  }
}
