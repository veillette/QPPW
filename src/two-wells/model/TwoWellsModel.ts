/**
 * TwoWellsModel represents the physics model for a double quantum potential well.
 * It handles quantum tunneling and double-well dynamics.
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

export class TwoWellsModel extends BaseModel {
  // Model-specific well parameters
  public readonly wellSeparationProperty: NumberProperty;

  // Barrier parameters
  public readonly barrierHeightProperty: NumberProperty;
  public readonly barrierWidthProperty: NumberProperty;

  // Tunneling visualization
  public readonly tunnelingProbabilityProperty: NumberProperty;

  // Model-specific display settings
  public readonly displayModeProperty: Property<DisplayMode>; // Extends base with "phaseColor" mode
  public readonly energyLevelProperty: NumberProperty; // Deprecated, kept for compatibility

  public constructor() {
    super();

    // Override potential type to default to double square well
    this.potentialTypeProperty.value = PotentialType.DOUBLE_SQUARE_WELL;

    // Override well width range for double square well
    this.wellWidthProperty.setValueAndRange(1.0, new Range(0.1, 3.0));

    // Initialize model-specific well parameters
    this.wellSeparationProperty = new NumberProperty(0.2, {
      range: new Range(0.05, 0.7),
    }); // in nanometers (0.05-0.7 nm for double square well)

    // Initialize barrier parameters
    this.barrierHeightProperty = new NumberProperty(3, {
      range: new Range(0.1, 15.0),
    }); // in eV
    this.barrierWidthProperty = new NumberProperty(2, {
      range: new Range(0.1, 5.0),
    }); // in nanometers

    // Initialize tunneling probability
    this.tunnelingProbabilityProperty = new NumberProperty(0);

    // Initialize model-specific display settings
    this.displayModeProperty = new Property<DisplayMode>("probabilityDensity");
    this.energyLevelProperty = new NumberProperty(0); // Deprecated

    // Override superposition config default
    this.superpositionConfigProperty.value = {
      type: SuperpositionType.PSI_I_PSI_J,
      amplitudes: [0.7, 0.7], // Default to equal superposition of first two states (normalized)
      phases: [0, 0],
    };
  }

  /**
   * Setup cache invalidation listeners, including model-specific properties.
   * Extends the base implementation to add TwoWellsModel-specific invalidation.
   */
  protected override setupCacheInvalidation(): void {
    super.setupCacheInvalidation();

    const invalidateCache = () => {
      this.boundStateResult = null;
    };

    this.wellSeparationProperty.link(invalidateCache);
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
    this.wellSeparationProperty.reset();
    this.barrierHeightProperty.reset();
    this.barrierWidthProperty.reset();
    this.tunnelingProbabilityProperty.reset();
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
    if (this.isPlayingProperty.value || forced) {
      this.updateTunnelingProbability();
      // Add quantum tunneling dynamics here
    }
  }

  /**
   * Updates the tunneling probability based on current parameters.
   * Uses WKB approximation for quantum tunneling.
   */
  private updateTunnelingProbability(): void {
    const energy = this.energyLevelProperty.value;
    const barrierHeight = this.barrierHeightProperty.value;
    const barrierWidth = this.barrierWidthProperty.value * 1e-9; // Convert to meters

    if (energy >= barrierHeight) {
      // Classical regime - particle goes over the barrier
      this.tunnelingProbabilityProperty.value = 1.0;
    } else {
      // Quantum tunneling regime
      const hbar = 1.054571817e-34; // Reduced Planck constant (J·s)
      const electronMass = 9.10938356e-31; // Electron mass (kg)
      const eV = 1.602176634e-19; // Electron volt in joules

      const V0 = barrierHeight * eV;
      const E = energy * eV;
      const kappa = Math.sqrt((2 * electronMass * (V0 - E)) / (hbar * hbar));

      this.tunnelingProbabilityProperty.value = Math.exp(
        -2 * kappa * barrierWidth,
      );
    }
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

    // Calculate number of states based on potential type and energy range
    let numStates = 10; // Default for most potentials

    // For infinite well, calculate states up to 15 eV
    if (this.potentialTypeProperty.value === PotentialType.INFINITE_WELL) {
      const maxEnergy = 15 * QuantumConstants.EV_TO_JOULES; // 15 eV
      // E_n = (ℏ²π²n²)/(2mL²), solve for n
      const maxN = Math.floor(
        Math.sqrt(
          (2 * mass * wellWidth * wellWidth * maxEnergy) /
            (QuantumConstants.HBAR * QuantumConstants.HBAR * Math.PI * Math.PI),
        ),
      );
      numStates = Math.max(1, Math.min(maxN, 100)); // Cap at 100 for safety
    } else if (this.potentialTypeProperty.value === PotentialType.COULOMB_1D) {
      numStates = 80; // Use more states for Coulomb potential
    } else if (
      this.potentialTypeProperty.value === PotentialType.DOUBLE_SQUARE_WELL
    ) {
      numStates = 80; // Use more states for double well to capture splitting
    }

    // Grid configuration
    let gridConfig;

    // For all potentials, span the full chart display range (-4 nm to +4 nm)
    const CHART_DISPLAY_RANGE_NM = 4;

    if (this.potentialTypeProperty.value === PotentialType.DOUBLE_SQUARE_WELL) {
      // For double square well, use analytical solution with fixed 2000 grid points
      // spanning the full chart range (-4 nm to +4 nm)
      gridConfig = {
        xMin: -CHART_DISPLAY_RANGE_NM * QuantumConstants.NM_TO_M,
        xMax: CHART_DISPLAY_RANGE_NM * QuantumConstants.NM_TO_M,
        numPoints: 2000,
      };
    } else {
      const method = this.solver.getNumericalMethod();
      let numGridPoints = QPPWPreferences.gridPointsProperty.value;

      if (method === "fgh") {
        // FGH: round to nearest power of 2 for FFT efficiency
        numGridPoints = Math.pow(2, Math.round(Math.log2(numGridPoints)));
      }

      gridConfig = {
        xMin: -CHART_DISPLAY_RANGE_NM * QuantumConstants.NM_TO_M,
        xMax: CHART_DISPLAY_RANGE_NM * QuantumConstants.NM_TO_M,
        numPoints: numGridPoints,
      };
    }

    try {
      // Build potential parameters based on type
      const potentialParams: WellParameters = {
        type: this.potentialTypeProperty.value,
        wellWidth: wellWidth,
      };

      // Add type-specific parameters
      switch (this.potentialTypeProperty.value) {
        case PotentialType.INFINITE_WELL:
          // No additional parameters needed
          break;
        case PotentialType.COULOMB_1D: {
          // For Coulomb potentials, use coulombStrength parameter α = k*e²
          // where k = 1/(4πε₀) ≈ 8.9875517923e9 N·m²/C²
          // α ≈ 2.307e-28 J·m for electron charge
          const coulombConstant = 8.9875517923e9; // Coulomb's constant in N·m²/C²
          potentialParams.coulombStrength =
            coulombConstant *
            QuantumConstants.ELEMENTARY_CHARGE *
            QuantumConstants.ELEMENTARY_CHARGE;
          break;
        }
        case PotentialType.DOUBLE_SQUARE_WELL: {
          // For double square well, we need width, depth, and separation
          potentialParams.wellDepth =
            this.wellDepthProperty.value * QuantumConstants.EV_TO_JOULES;
          potentialParams.wellSeparation =
            this.wellSeparationProperty.value * QuantumConstants.NM_TO_M;
          break;
        }
        default:
          // For other potential types, use numerical solution
          break;
      }

      // Attempt analytical solution first
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
        case PotentialType.INFINITE_WELL:
          // V = 0 inside [-L/2, L/2], infinity outside
          V = Math.abs(x) <= wellWidth / 2 ? 0 : Infinity;
          break;

        case PotentialType.DOUBLE_SQUARE_WELL: {
          // Double square well: two wells separated by a barrier
          const halfWellWidth = wellWidth / 2;
          const halfSeparation = wellSeparation / 2;

          // Left well: centered at -halfSeparation
          const leftWellStart = -halfSeparation - halfWellWidth;
          const leftWellEnd = -halfSeparation + halfWellWidth;

          // Right well: centered at +halfSeparation
          const rightWellStart = halfSeparation - halfWellWidth;
          const rightWellEnd = halfSeparation + halfWellWidth;

          if (
            (x >= leftWellStart && x <= leftWellEnd) ||
            (x >= rightWellStart && x <= rightWellEnd)
          ) {
            V = 0; // Inside wells
          } else {
            V = wellDepth; // Outside wells (barrier or exterior)
          }
          break;
        }

        case PotentialType.COULOMB_1D: {
          // V(x) = -α/|x| where α = ke²
          const coulombConstant = 8.9875517923e9;
          const coulombStrength =
            coulombConstant *
            QuantumConstants.ELEMENTARY_CHARGE *
            QuantumConstants.ELEMENTARY_CHARGE;
          const r = Math.abs(x);
          if (r > 1e-12) {
            // Avoid singularity at origin
            V = -coulombStrength / r;
          } else {
            V = -coulombStrength / 1e-12;
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
