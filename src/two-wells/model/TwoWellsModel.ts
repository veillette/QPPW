/**
 * TwoWellsModel represents the physics model for a double quantum potential well.
 * It handles quantum tunneling and double-well dynamics.
 */

import { NumberProperty, Property } from "scenerystack/axon";
import { Range } from "scenerystack/dot";
import { BaseModel } from "../../common/model/BaseModel.js";
import Schrodinger1DSolver, { WellParameters, NumericalMethod } from "../../common/model/Schrodinger1DSolver.js";
import { PotentialType, BoundStateResult } from "../../common/model/PotentialFunction.js";
import QuantumConstants from "../../common/model/QuantumConstants.js";
import { SuperpositionType, SuperpositionConfig } from "../../common/model/SuperpositionType.js";

export type DisplayMode = "probabilityDensity" | "waveFunction" | "phaseColor";

export class TwoWellsModel extends BaseModel {
  // Potential type selection
  public readonly potentialTypeProperty: Property<PotentialType>;

  // Well parameters (used for different potential types)
  public readonly wellWidthProperty: NumberProperty;
  public readonly wellDepthProperty: NumberProperty;
  public readonly wellOffsetProperty: NumberProperty; // For asymmetric wells
  public readonly wellSeparationProperty: NumberProperty;

  // Barrier parameters
  public readonly barrierHeightProperty: NumberProperty;
  public readonly barrierWidthProperty: NumberProperty;

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

  // Chart visibility
  public readonly showTotalEnergyProperty: Property<boolean>;
  public readonly showPotentialEnergyProperty: Property<boolean>;

  // Tunneling visualization
  public readonly tunnelingProbabilityProperty: NumberProperty;

  // Superposition state
  public readonly superpositionTypeProperty: Property<SuperpositionType>;
  public readonly superpositionConfigProperty: Property<SuperpositionConfig>;

  // Cached bound state results
  protected boundStateResult: BoundStateResult | null = null;

  public constructor() {
    super();
    // Initialize potential type (double square well by default)
    this.potentialTypeProperty = new Property<PotentialType>(PotentialType.DOUBLE_SQUARE_WELL);

    // Initialize well parameters with default values
    this.wellWidthProperty = new NumberProperty(1.0, { range: new Range(0.1, 3.0) }); // in nanometers (0.1-3 nm for double square well)
    this.wellDepthProperty = new NumberProperty(5.0, { range: new Range(0.1, 15.0) }); // in eV (0.1-15 eV for double square well)
    this.wellOffsetProperty = new NumberProperty(0.5, { range: new Range(0.0, 1.0) }); // normalized position
    this.wellSeparationProperty = new NumberProperty(0.2, { range: new Range(0.05, 0.7) }); // in nanometers (0.05-0.7 nm for double square well)

    // Initialize barrier parameters
    this.barrierHeightProperty = new NumberProperty(3); // in eV
    this.barrierWidthProperty = new NumberProperty(2); // in nanometers

    // Initialize particle mass (1.0 = electron mass)
    this.particleMassProperty = new NumberProperty(1.0, { range: new Range(0.5, 1.1) }); // 0.5 to 1.1 times electron mass

    // Initialize energy level selection (ground state by default)
    // Use a large range to accommodate potentials with many states (up to 100)
    this.selectedEnergyLevelIndexProperty = new NumberProperty(0, { range: new Range(0, 99) });
    this.energyLevelProperty = new NumberProperty(0); // Deprecated

    // Initialize display settings
    this.displayModeProperty = new Property<DisplayMode>("probabilityDensity");
    this.showRealPartProperty = new Property<boolean>(true);
    this.showImaginaryPartProperty = new Property<boolean>(false);
    this.showMagnitudeProperty = new Property<boolean>(false);
    this.showPhaseProperty = new Property<boolean>(false);

    // Initialize chart visibility
    this.showTotalEnergyProperty = new Property<boolean>(true);
    this.showPotentialEnergyProperty = new Property<boolean>(true);

    // Initialize tunneling probability
    this.tunnelingProbabilityProperty = new NumberProperty(0);

    // Initialize superposition state
    this.superpositionTypeProperty = new Property<SuperpositionType>(SuperpositionType.PSI_I_PSI_J);
    this.superpositionConfigProperty = new Property<SuperpositionConfig>({
      type: SuperpositionType.PSI_I_PSI_J,
      amplitudes: [0.7, 0.7], // Default to equal superposition of first two states (normalized)
      phases: [0, 0],
    });

    // Recalculate bound states when parameters change
    const invalidateCache = () => {
      this.boundStateResult = null;
    };

    this.potentialTypeProperty.link(invalidateCache);
    this.wellWidthProperty.link(invalidateCache);
    this.wellDepthProperty.link(invalidateCache);
    this.wellOffsetProperty.link(invalidateCache);
    this.wellSeparationProperty.link(invalidateCache);
    this.particleMassProperty.link(invalidateCache);
  }

  /**
   * Called when the solver method changes.
   * Invalidates the cached bound state results.
   * @param _method - The new numerical method (unused but required by interface)
   */
  // eslint-disable-next-line @typescript-eslint/no-unused-vars
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
    this.wellWidthProperty.reset();
    this.wellDepthProperty.reset();
    this.wellOffsetProperty.reset();
    this.wellSeparationProperty.reset();
    this.barrierHeightProperty.reset();
    this.barrierWidthProperty.reset();
    this.particleMassProperty.reset();
    this.selectedEnergyLevelIndexProperty.reset();
    this.energyLevelProperty.reset();
    this.displayModeProperty.reset();
    this.showRealPartProperty.reset();
    this.showImaginaryPartProperty.reset();
    this.showMagnitudeProperty.reset();
    this.showPhaseProperty.reset();
    this.showTotalEnergyProperty.reset();
    this.showPotentialEnergyProperty.reset();
    this.tunnelingProbabilityProperty.reset();
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
      const transmissionCoeff = Math.exp(-2 * kappa * barrierWidth);

      this.tunnelingProbabilityProperty.value = transmissionCoeff;
    }
  }

  /**
   * Calculates the energy eigenvalues for the current well parameters.
   * Uses the Schrödinger solver with analytical solution for infinite well.
   * @param n - The quantum number (1, 2, 3, ...)
   * @returns The energy in eV
   */
  public getEnergyLevel(n: number): number {
    // Ensure bound states are calculated
    if (!this.boundStateResult) {
      this.calculateBoundStates();
    }

    // Return energy for quantum number n (1-indexed)
    if (this.boundStateResult && n > 0 && n <= this.boundStateResult.energies.length) {
      const energyJoules = this.boundStateResult.energies[n - 1];
      return Schrodinger1DSolver.joulesToEV(energyJoules);
    }

    // Fallback to analytical formula if solver fails
    const L = this.wellWidthProperty.value * QuantumConstants.NM_TO_M;
    const energy =
      (n * n * Math.PI * Math.PI * QuantumConstants.HBAR * QuantumConstants.HBAR) /
      (2 * QuantumConstants.ELECTRON_MASS * L * L);
    return Schrodinger1DSolver.joulesToEV(energy);
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
    const mass = this.particleMassProperty.value * QuantumConstants.ELECTRON_MASS;

    // Calculate number of states based on potential type and energy range
    let numStates = 10; // Default for most potentials

    // For infinite well, calculate states up to 15 eV
    if (this.potentialTypeProperty.value === PotentialType.INFINITE_WELL) {
      const maxEnergy = 15 * QuantumConstants.EV_TO_JOULES; // 15 eV
      // E_n = (ℏ²π²n²)/(2mL²), solve for n
      const maxN = Math.floor(Math.sqrt((2 * mass * wellWidth * wellWidth * maxEnergy) / (QuantumConstants.HBAR * QuantumConstants.HBAR * Math.PI * Math.PI)));
      numStates = Math.max(1, Math.min(maxN, 100)); // Cap at 100 for safety
    }
    else if (this.potentialTypeProperty.value === PotentialType.COULOMB_1D) {
      numStates = 30; // Use more states for Coulomb potential
    }
    else if (this.potentialTypeProperty.value === PotentialType.DOUBLE_SQUARE_WELL) {
      numStates = 20; // Use more states for double well to capture splitting
    }

    // Grid configuration
    let gridConfig;

    if (this.potentialTypeProperty.value === PotentialType.DOUBLE_SQUARE_WELL) {
      // For double square well, use focused grid to keep kinetic energy manageable
      const method = this.solver.getNumericalMethod();
      const separation = this.wellSeparationProperty.value;
      const wellCenter = (separation / 2 + wellWidth / (2 * QuantumConstants.NM_TO_M));

      // For Numerov, use minimal grid: just the wells plus tiny margin
      const gridRange = method === 'numerov'
        ? wellCenter + 0.3 * (wellWidth / QuantumConstants.NM_TO_M)
        : wellCenter + 1.5 * (wellWidth / QuantumConstants.NM_TO_M);

      if (method === 'dvr') {
        // DVR: use minimal points to keep kinetic energy low
        gridConfig = {
          xMin: -gridRange * QuantumConstants.NM_TO_M,
          xMax: gridRange * QuantumConstants.NM_TO_M,
          numPoints: 20,
        };
      } else if (method === 'spectral') {
        // Spectral: use more points than DVR for smooth wavefunctions
        gridConfig = {
          xMin: -gridRange * QuantumConstants.NM_TO_M,
          xMax: gridRange * QuantumConstants.NM_TO_M,
          numPoints: 60,
        };
      } else if (method === 'fgh') {
        // FGH: use power-of-2 points
        gridConfig = {
          xMin: -gridRange * QuantumConstants.NM_TO_M,
          xMax: gridRange * QuantumConstants.NM_TO_M,
          numPoints: 64,
        };
      } else {
        // Numerov: use focused grid with fewer points to avoid numerical instability
        gridConfig = {
          xMin: -gridRange * QuantumConstants.NM_TO_M,
          xMax: gridRange * QuantumConstants.NM_TO_M,
          numPoints: 200,
        };
      }
    } else {
      // For other potentials, span the full chart display range (-4 nm to +4 nm)
      const CHART_DISPLAY_RANGE_NM = 4;
      const method = this.solver.getNumericalMethod();
      let numGridPoints = 100;

      if (method === 'numerov') {
        numGridPoints = 1000;
      } else if (method === 'spectral' || method === 'fgh') {
        numGridPoints = 100; // Keep moderate for spectral methods
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
          potentialParams.coulombStrength = coulombConstant * QuantumConstants.ELEMENTARY_CHARGE * QuantumConstants.ELEMENTARY_CHARGE;
          break;
        }
        case PotentialType.DOUBLE_SQUARE_WELL: {
          // For double square well, we need width, depth, and separation
          potentialParams.wellDepth = this.wellDepthProperty.value * QuantumConstants.EV_TO_JOULES;
          potentialParams.wellSeparation = this.wellSeparationProperty.value * QuantumConstants.NM_TO_M;
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
   * Get the wavefunction for a specific quantum number.
   * @param n - The quantum number (1, 2, 3, ...)
   * @returns Array of wavefunction values at grid points, or null if unavailable
   */
  public getWavefunction(n: number): number[] | null {
    if (!this.boundStateResult) {
      this.calculateBoundStates();
    }

    if (this.boundStateResult && n > 0 && n <= this.boundStateResult.wavefunctions.length) {
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
      return this.boundStateResult.xGrid.map((x) => x * QuantumConstants.M_TO_NM);
    }

    return null;
  }
}
