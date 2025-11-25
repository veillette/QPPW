/**
 * ManyWellsModel represents the physics model for multiple quantum potential wells.
 * Similar to TwoWellsModel but generalized to N wells (1-10).
 * Supports multi-square wells and multi-Coulomb 1D potentials.
 */

import { NumberProperty, Property } from "scenerystack/axon";
import { Range } from "scenerystack/dot";
import { BaseModel } from "../../common/model/BaseModel.js";
import Schrodinger1DSolver, {
  WellParameters,
  NumericalMethod,
} from "../../common/model/Schrodinger1DSolver.js";
import {
  PotentialType,
  BoundStateResult,
} from "../../common/model/PotentialFunction.js";
import QuantumConstants from "../../common/model/QuantumConstants.js";
import {
  SuperpositionType,
  SuperpositionConfig,
} from "../../common/model/SuperpositionType.js";
import QPPWPreferences from "../../QPPWPreferences.js";

export type DisplayMode = "probabilityDensity" | "waveFunction" | "phaseColor";

export class ManyWellsModel extends BaseModel {
  // Potential type selection (limited to multi-well potentials)
  public readonly potentialTypeProperty: Property<PotentialType>;

  // Number of wells (1-10)
  public readonly numberOfWellsProperty: NumberProperty;

  // Well parameters (used for different potential types)
  public readonly wellWidthProperty: NumberProperty;
  public readonly wellDepthProperty: NumberProperty;
  public readonly wellSeparationProperty: NumberProperty;

  // Electric field (eV/nm) - creates a linear tilt in the potential
  public readonly electricFieldProperty: NumberProperty;

  // Particle properties
  public readonly particleMassProperty: NumberProperty; // In units of electron mass

  // Energy level selection
  public readonly selectedEnergyLevelIndexProperty: NumberProperty; // 0-indexed
  public readonly energyLevelProperty: NumberProperty; // Deprecated, kept for compatibility

  // Display settings
  public readonly displayModeProperty: Property<DisplayMode>;
  public readonly showRealPartProperty: Property<boolean>;
  public readonly showImaginaryPartProperty: Property<boolean>;
  public readonly showMagnitudeProperty: Property<boolean>;
  public readonly showPhaseProperty: Property<boolean>;
  public readonly showClassicalProbabilityProperty: Property<boolean>;

  // Chart visibility
  public readonly showTotalEnergyProperty: Property<boolean>;
  public readonly showPotentialEnergyProperty: Property<boolean>;

  // Superposition state
  public readonly superpositionTypeProperty: Property<SuperpositionType>;
  public readonly superpositionConfigProperty: Property<SuperpositionConfig>;

  // Cached bound state results
  protected boundStateResult: BoundStateResult | null = null;

  public constructor() {
    super();

    // Initialize potential type (multi-square well by default)
    this.potentialTypeProperty = new Property<PotentialType>(
      PotentialType.MULTI_SQUARE_WELL,
    );

    // Initialize number of wells (default: 3 wells)
    this.numberOfWellsProperty = new NumberProperty(3, {
      range: new Range(1, 10),
    });

    // Initialize well parameters with default values
    this.wellWidthProperty = new NumberProperty(1.0, {
      range: new Range(0.1, 3.0),
    }); // in nanometers
    this.wellDepthProperty = new NumberProperty(5.0, {
      range: new Range(0.1, 15.0),
    }); // in eV
    this.wellSeparationProperty = new NumberProperty(0.2, {
      range: new Range(0.05, 0.7),
    }); // in nanometers (spacing between wells)

    // Initialize electric field (default: 0, range: ±0.625 eV/nm for ±5eV tilt over 8nm)
    this.electricFieldProperty = new NumberProperty(0.0, {
      range: new Range(-0.625, 0.625),
    }); // in eV/nm

    // Initialize particle mass (1.0 = electron mass)
    this.particleMassProperty = new NumberProperty(1.0, {
      range: new Range(0.5, 1.1),
    }); // 0.5 to 1.1 times electron mass

    // Initialize energy level selection (ground state by default)
    this.selectedEnergyLevelIndexProperty = new NumberProperty(0, {
      range: new Range(0, 99),
    });
    this.energyLevelProperty = new NumberProperty(0); // Deprecated

    // Initialize display settings
    this.displayModeProperty = new Property<DisplayMode>("probabilityDensity");
    this.showRealPartProperty = new Property<boolean>(true);
    this.showImaginaryPartProperty = new Property<boolean>(false);
    this.showMagnitudeProperty = new Property<boolean>(false);
    this.showPhaseProperty = new Property<boolean>(false);
    this.showClassicalProbabilityProperty = new Property<boolean>(false);

    // Initialize chart visibility
    this.showTotalEnergyProperty = new Property<boolean>(true);
    this.showPotentialEnergyProperty = new Property<boolean>(true);

    // Initialize superposition state
    this.superpositionTypeProperty = new Property<SuperpositionType>(
      SuperpositionType.SINGLE,
    );
    this.superpositionConfigProperty = new Property<SuperpositionConfig>({
      type: SuperpositionType.PSI_I_PSI_J,
      amplitudes: [0.7, 0.7], // Default to equal superposition
      phases: [0, 0],
    });

    // Recalculate bound states when parameters change
    const invalidateCache = () => {
      this.boundStateResult = null;
    };

    this.potentialTypeProperty.link(invalidateCache);
    this.numberOfWellsProperty.link(invalidateCache);
    this.wellWidthProperty.link(invalidateCache);
    this.wellDepthProperty.link(invalidateCache);
    this.wellSeparationProperty.link(invalidateCache);
    this.electricFieldProperty.link(invalidateCache);
    this.particleMassProperty.link(invalidateCache);
  }

  /**
   * Called when the solver method changes.
   * Invalidates the cached bound state results.
   * @param _method - The new numerical method (unused but required by interface)
   */
  protected onSolverMethodChanged(_method: NumericalMethod): void {
    this.boundStateResult = null; // Invalidate cache
  }

  /**
   * Resets the model to its initial state.
   * This is the public API method that delegates to resetAll().
   */
  public reset(): void {
    this.resetAll();
  }

  /**
   * Resets all properties to their initial state.
   * Override from BaseModel.
   */
  public override resetAll(): void {
    super.resetAll();
    this.potentialTypeProperty.reset();
    this.numberOfWellsProperty.reset();
    this.wellWidthProperty.reset();
    this.wellDepthProperty.reset();
    this.wellSeparationProperty.reset();
    this.particleMassProperty.reset();
    this.selectedEnergyLevelIndexProperty.reset();
    this.energyLevelProperty.reset();
    this.displayModeProperty.reset();
    this.showRealPartProperty.reset();
    this.showImaginaryPartProperty.reset();
    this.showMagnitudeProperty.reset();
    this.showPhaseProperty.reset();
    this.showClassicalProbabilityProperty.reset();
    this.showTotalEnergyProperty.reset();
    this.showPotentialEnergyProperty.reset();
    this.superpositionTypeProperty.reset();
    this.superpositionConfigProperty.reset();
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
   * Calculates the energy eigenvalues for the current well parameters.
   * @param n - The quantum number (1, 2, 3, ...)
   * @returns The energy in eV
   */
  public getEnergyLevel(n: number): number {
    // Ensure bound states are calculated
    if (!this.boundStateResult) {
      this.calculateBoundStates();
    }

    // Return energy for quantum number n (1-indexed)
    if (
      this.boundStateResult &&
      n > 0 &&
      n <= this.boundStateResult.energies.length
    ) {
      const energyJoules = this.boundStateResult.energies[n - 1];
      return Schrodinger1DSolver.joulesToEV(energyJoules);
    }

    return 0;
  }

  /**
   * Get all bound state energies and wavefunctions for the current well.
   * @returns BoundStateResult with energies (in eV) and wavefunctions
   */
  public getBoundStates(): BoundStateResult | null {
    if (!this.boundStateResult) {
      this.calculateBoundStates();
    }
    return this.boundStateResult;
  }

  /**
   * Calculate bound states using the Schrödinger solver.
   * Results are cached until well parameters change.
   */
  private calculateBoundStates(): void {
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
   * Get the wavefunction for a specific quantum number.
   * @param n - The quantum number (1, 2, 3, ...)
   * @returns Array of wavefunction values at grid points, or null if unavailable
   */
  public getWavefunction(n: number): number[] | null {
    if (!this.boundStateResult) {
      this.calculateBoundStates();
    }

    if (
      this.boundStateResult &&
      n > 0 &&
      n <= this.boundStateResult.wavefunctions.length
    ) {
      return this.boundStateResult.wavefunctions[n - 1];
    }

    return null;
  }

  /**
   * Get the spatial grid points for wavefunctions.
   * @returns Array of x positions in nanometers
   */
  public getXGrid(): number[] | null {
    if (!this.boundStateResult) {
      this.calculateBoundStates();
    }

    if (this.boundStateResult) {
      // Convert from meters to nanometers
      return this.boundStateResult.xGrid.map(
        (x) => x * QuantumConstants.M_TO_NM,
      );
    }

    return null;
  }

  /**
   * Calculate the classical probability density for a given energy level.
   * The classical probability density is inversely proportional to the velocity:
   * P(x) ∝ 1/v(x) = 1/√[2(E - V(x))/m]
   *
   * @param energyIndex - Index of the energy level (0-indexed)
   * @returns Array of classical probability density values, or null if unavailable
   */
  public getClassicalProbabilityDensity(energyIndex: number): number[] | null {
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

    // Calculate classical probability density
    const classicalProbability: number[] = [];
    let integralSum = 0;

    // First pass: calculate unnormalized probability and find turning points
    for (let i = 0; i < xGrid.length; i++) {
      const kineticEnergy = energy - potential[i];

      // Classical turning points: where E = V(x)
      // In classical mechanics, the particle cannot exist beyond turning points
      if (kineticEnergy <= 0) {
        classicalProbability.push(0);
      } else {
        // P(x) ∝ 1/v(x) = 1/√[2(E - V(x))/m] = √[m/(2(E - V(x)))]
        const probability = 1 / Math.sqrt(2 * kineticEnergy / mass);
        classicalProbability.push(probability);

        // For normalization (using trapezoidal rule)
        if (i > 0) {
          const dx = xGrid[i] - xGrid[i - 1];
          integralSum += (probability + classicalProbability[i - 1]) * dx / 2;
        }
      }
    }

    // Normalize so that ∫P(x)dx = 1
    if (integralSum > 0) {
      for (let i = 0; i < classicalProbability.length; i++) {
        classicalProbability[i] /= integralSum;
      }
    }

    return classicalProbability;
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
          const totalExtent = numberOfWells * wellWidth + (numberOfWells - 1) * wellSeparation;
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
