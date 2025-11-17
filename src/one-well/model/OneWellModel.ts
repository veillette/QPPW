/**
 * OneWellModel represents the physics model for a single quantum potential well.
 * It handles the quantum mechanical calculations for a particle in a single well.
 */

import { NumberProperty, Property } from "scenerystack/axon";
import { Range } from "scenerystack/dot";
import { BaseModel } from "../../common/model/BaseModel.js";
import Schrodinger1DSolver, { WellParameters, NumericalMethod } from "../../common/model/Schrodinger1DSolver.js";
import { PotentialType, BoundStateResult } from "../../common/model/PotentialFunction.js";
import QuantumConstants from "../../common/model/QuantumConstants.js";
import { SuperpositionType, SuperpositionConfig } from "../../common/model/SuperpositionType.js";
import { calculateCoherentStateCoefficients } from "../../common/model/analytical-solutions/harmonic-oscillator.js";
import QPPWPreferences from "../../QPPWPreferences.js";

export type DisplayMode = "probabilityDensity" | "waveFunction" | "phaseColor";

export class OneWellModel extends BaseModel {
  // Potential type selection
  public readonly potentialTypeProperty: Property<PotentialType>;

  // Well parameters (used for different potential types)
  public readonly wellWidthProperty: NumberProperty;
  public readonly wellDepthProperty: NumberProperty;
  public readonly wellOffsetProperty: NumberProperty; // For asymmetric wells
  public readonly barrierHeightProperty: NumberProperty; // For Rosen-Morse and Eckart potentials

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

  // Superposition state
  public readonly superpositionTypeProperty: Property<SuperpositionType>;
  public readonly superpositionConfigProperty: Property<SuperpositionConfig>;

  // Coherent state parameter
  public readonly coherentDisplacementProperty: NumberProperty; // Displacement in nm

  // Cached bound state results
  protected boundStateResult: BoundStateResult | null = null;

  public constructor() {
    super();
    // Initialize potential type (square/infinite well by default)
    this.potentialTypeProperty = new Property<PotentialType>(PotentialType.INFINITE_WELL);

    // Initialize well parameters with default values
    this.wellWidthProperty = new NumberProperty(1.0, { range: new Range(0.1, 6.0) }); // in nanometers (max 6 nm)
    this.wellDepthProperty = new NumberProperty(5.0, { range: new Range(0.1, 15.0) }); // in eV (within energy graph bounds)
    this.wellOffsetProperty = new NumberProperty(0.5, { range: new Range(0.0, 1.0) }); // normalized position
    this.barrierHeightProperty = new NumberProperty(2.0, { range: new Range(0.0, 10.0) }); // in eV (for Rosen-Morse and Eckart)

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

    // Initialize superposition state
    this.superpositionTypeProperty = new Property<SuperpositionType>(SuperpositionType.PSI_I_PSI_J);
    this.superpositionConfigProperty = new Property<SuperpositionConfig>({
      type: SuperpositionType.PSI_I_PSI_J,
      amplitudes: [0.7, 0.7], // Default to equal superposition of first two states (normalized)
      phases: [0, 0],
    });

    // Initialize coherent state displacement (0.5 nm default, range 0 to 2 nm)
    this.coherentDisplacementProperty = new NumberProperty(0.5, { range: new Range(0.0, 2.0) });

    // Recalculate bound states when parameters change
    const invalidateCache = () => {
      this.boundStateResult = null;
      // Also update superposition coefficients when bound states change
      this.updateSuperpositionCoefficients();
    };

    this.potentialTypeProperty.link(invalidateCache);
    this.wellWidthProperty.link(invalidateCache);
    this.wellDepthProperty.link(invalidateCache);
    this.wellOffsetProperty.link(invalidateCache);
    this.barrierHeightProperty.link(invalidateCache);
    this.particleMassProperty.link(invalidateCache);

    // Update superposition coefficients when superposition type or coherent displacement changes
    this.superpositionTypeProperty.link(() => this.updateSuperpositionCoefficients());
    this.coherentDisplacementProperty.link(() => {
      if (this.superpositionTypeProperty.value === SuperpositionType.COHERENT) {
        this.updateSuperpositionCoefficients();
      }
    });
  }

  /**
   * Called when the solver method changes.
   * For OneWellModel, we ignore solver method changes because all potentials
   * in the one-well screen have analytical solutions and don't depend on the
   * numerical method preference.
   * @param _method - The new numerical method (unused)
   */
  // eslint-disable-next-line @typescript-eslint/no-unused-vars
  protected onSolverMethodChanged(_method: NumericalMethod): void {
    // Do nothing - one well screen always uses analytical solutions
    // regardless of the solver preference setting
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
    this.barrierHeightProperty.reset();
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
    this.superpositionTypeProperty.reset();
    this.superpositionConfigProperty.reset();
    this.coherentDisplacementProperty.reset();
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
    const wellDepth = this.wellDepthProperty.value * QuantumConstants.EV_TO_JOULES;
    const mass = this.particleMassProperty.value * QuantumConstants.ELECTRON_MASS;

    // Calculate number of states based on potential type and energy range
    let numStates = 10; // Default for most potentials

    // For harmonic oscillator, calculate states up to 15 eV
    if (this.potentialTypeProperty.value === PotentialType.HARMONIC_OSCILLATOR) {
      const springConstant = (4 * wellDepth) / (wellWidth * wellWidth);
      const omega = Math.sqrt(springConstant / mass);
      const maxEnergy = 15 * QuantumConstants.EV_TO_JOULES; // 15 eV
      // E_n = ℏω(n + 1/2), solve for n: n = E/(ℏω) - 1/2
      const maxN = Math.floor(maxEnergy / (QuantumConstants.HBAR * omega) - 0.5);
      numStates = Math.max(1, Math.min(maxN + 1, 100)); // Cap at 100 for safety
    }
    // For infinite well, calculate states up to 15 eV
    else if (this.potentialTypeProperty.value === PotentialType.INFINITE_WELL) {
      const maxEnergy = 15 * QuantumConstants.EV_TO_JOULES; // 15 eV
      // E_n = (ℏ²π²n²)/(2mL²), solve for n
      const maxN = Math.floor(Math.sqrt((2 * mass * wellWidth * wellWidth * maxEnergy) / (QuantumConstants.HBAR * QuantumConstants.HBAR * Math.PI * Math.PI)));
      numStates = Math.max(1, Math.min(maxN, 100)); // Cap at 100 for safety
    }
    // For finite well, estimate maximum number of bound states
    else if (this.potentialTypeProperty.value === PotentialType.FINITE_WELL) {
      // Approximate number of bound states: n_max ≈ (1/π) * sqrt(2mV₀L²/ℏ²)
      // Use generous estimate to ensure we get all states
      const estimatedMax = Math.ceil((1 / Math.PI) * Math.sqrt((2 * mass * wellDepth * wellWidth * wellWidth) / (QuantumConstants.HBAR * QuantumConstants.HBAR)));
      // Request more states than estimated to ensure we capture all bound states
      numStates = Math.max(10, Math.min(estimatedMax * 2, 100)); // At least 10, cap at 100
    }
    // For asymmetric triangle, calculate states that fit in the energy range
    else if (this.potentialTypeProperty.value === PotentialType.ASYMMETRIC_TRIANGLE) {
      numStates = 80; // Asymmetric triangle may have many states, use larger number
    }
    else {
      numStates = 80; // Use more states for other potentials
    }

    // Grid configuration spans the full chart display range (-4 nm to +4 nm)
    // This ensures wavefunctions are calculated across the entire visible area
    const CHART_DISPLAY_RANGE_NM = 4; // Match the chart x-axis range
    const gridConfig = {
      xMin: -CHART_DISPLAY_RANGE_NM * QuantumConstants.NM_TO_M,
      xMax: CHART_DISPLAY_RANGE_NM * QuantumConstants.NM_TO_M,
      numPoints: QPPWPreferences.gridPointsProperty.value,
    };

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
        case PotentialType.FINITE_WELL:
          potentialParams.wellDepth = wellDepth;
          break;
        case PotentialType.HARMONIC_OSCILLATOR:
          // Convert well depth to spring constant: k = mω² = m(4V₀/mL²) = 4V₀/L²
          potentialParams.springConstant = (4 * wellDepth) / (wellWidth * wellWidth);
          break;
        case PotentialType.MORSE:
          // Morse potential uses wellDepth as dissociation energy, equilibrium position at 0
          potentialParams.dissociationEnergy = wellDepth;
          potentialParams.equilibriumPosition = 0; // Center of the well
          break;
        case PotentialType.POSCHL_TELLER:
          // Pöschl-Teller potential uses wellDepth as potential depth
          potentialParams.potentialDepth = wellDepth;
          break;
        case PotentialType.ROSEN_MORSE:
          // Rosen-Morse potential uses wellDepth as potential depth and barrierHeight
          potentialParams.potentialDepth = wellDepth;
          potentialParams.barrierHeight = this.barrierHeightProperty.value * QuantumConstants.EV_TO_JOULES;
          break;
        case PotentialType.ECKART:
          // Eckart potential uses wellDepth as potential depth and barrierHeight
          potentialParams.potentialDepth = wellDepth;
          potentialParams.barrierHeight = this.barrierHeightProperty.value * QuantumConstants.EV_TO_JOULES;
          break;
        case PotentialType.ASYMMETRIC_TRIANGLE:
          // Slope is the field strength
          potentialParams.slope = wellDepth / wellWidth;
          break;
        case PotentialType.COULOMB_1D: {
          // For 1D Coulomb potential: E_0 = -2mα²/ℏ² = -13.6 eV
          // Solving for α: α = √(13.6 eV * ℏ² / (2m))
          const groundStateEnergy = 13.6 * QuantumConstants.EV_TO_JOULES; // -13.6 eV in Joules
          potentialParams.coulombStrength = Math.sqrt(
            (groundStateEnergy * QuantumConstants.HBAR * QuantumConstants.HBAR) / (2 * mass)
          );
          break;
        }
        case PotentialType.COULOMB_3D: {
          // For 3D Coulomb potential (hydrogen): E_1 = -mα²/(2ℏ²) = -13.6 eV
          // Solving for α: α = √(2ℏ² * 13.6 eV / m)
          const groundStateEnergy = 13.6 * QuantumConstants.EV_TO_JOULES; // -13.6 eV in Joules
          potentialParams.coulombStrength = Math.sqrt(
            (2 * QuantumConstants.HBAR * QuantumConstants.HBAR * groundStateEnergy) / mass
          );
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

  /**
   * Get the superposition wavefunction based on current superposition config.
   * Returns the real-space wavefunction: ψ(x) = Σ c_n * e^(iφ_n) * ψ_n(x)
   * Since we're returning real values for display, this returns the real part at t=0.
   *
   * @returns Object with wavefunction array and energy (weighted average), or null
   */
  public getSuperpositionWavefunction(): { wavefunction: number[]; energy: number } | null {
    if (!this.boundStateResult) {
      this.calculateBoundStates();
    }

    if (!this.boundStateResult) {
      return null;
    }

    const config = this.superpositionConfigProperty.value;
    const numPoints = this.boundStateResult.xGrid.length;
    const wavefunction = new Array(numPoints).fill(0);
    let weightedEnergy = 0;

    // Compute superposition: ψ(x) = Σ c_n * e^(iφ_n) * ψ_n(x)
    // At t=0, this is: ψ(x) = Σ c_n * (cos(φ_n) + i*sin(φ_n)) * ψ_n(x)
    // We return the real part: Σ c_n * cos(φ_n) * ψ_n(x)
    for (let n = 0; n < config.amplitudes.length; n++) {
      const amplitude = config.amplitudes[n];
      const phase = config.phases[n];

      if (amplitude === 0) {
        continue;
      }

      if (n < this.boundStateResult.wavefunctions.length) {
        const eigenfunction = this.boundStateResult.wavefunctions[n];
        const coeff = amplitude * Math.cos(phase);

        for (let i = 0; i < numPoints; i++) {
          wavefunction[i] += coeff * eigenfunction[i];
        }

        // Add weighted energy contribution
        weightedEnergy += amplitude * amplitude * this.boundStateResult.energies[n];
      }
    }

    return {
      wavefunction,
      energy: weightedEnergy,
    };
  }

  /**
   * Updates the superposition coefficients based on the selected superposition type.
   * This method computes the amplitudes and phases for predefined superposition types.
   */
  public updateSuperpositionCoefficients(): void {
    const type = this.superpositionTypeProperty.value;

    // Skip CUSTOM type - coefficients are set manually by the user
    if (type === SuperpositionType.CUSTOM) {
      return;
    }

    // Ensure we have bound states
    if (!this.boundStateResult) {
      this.calculateBoundStates();
    }

    if (!this.boundStateResult) {
      return;
    }

    const numStates = this.boundStateResult.energies.length;
    let amplitudes: number[];
    let phases: number[];

    switch (type) {
      case SuperpositionType.PSI_I_PSI_J:
        // Equal superposition of first two states: (|0⟩ + |1⟩)/√2
        amplitudes = new Array(numStates).fill(0);
        phases = new Array(numStates).fill(0);
        if (numStates >= 2) {
          amplitudes[0] = 1 / Math.sqrt(2);
          amplitudes[1] = 1 / Math.sqrt(2);
        }
        break;

      case SuperpositionType.PSI_K:
        // Single eigenstate (ground state)
        amplitudes = new Array(numStates).fill(0);
        phases = new Array(numStates).fill(0);
        amplitudes[0] = 1;
        break;

      case SuperpositionType.LOCALIZED_NARROW: {
        // Narrow localized state: superposition of first few states
        amplitudes = new Array(numStates).fill(0);
        phases = new Array(numStates).fill(0);
        const numStatesNarrow = Math.min(5, numStates);
        for (let i = 0; i < numStatesNarrow; i++) {
          amplitudes[i] = 1;
        }
        // Normalize
        const normNarrow = Math.sqrt(amplitudes.reduce((sum, a) => sum + a * a, 0));
        amplitudes = amplitudes.map((a) => a / normNarrow);
        break;
      }

      case SuperpositionType.LOCALIZED_WIDE: {
        // Wide localized state: superposition of more states
        amplitudes = new Array(numStates).fill(0);
        phases = new Array(numStates).fill(0);
        const numStatesWide = Math.min(10, numStates);
        for (let i = 0; i < numStatesWide; i++) {
          amplitudes[i] = 1;
        }
        // Normalize
        const normWide = Math.sqrt(amplitudes.reduce((sum, a) => sum + a * a, 0));
        amplitudes = amplitudes.map((a) => a / normWide);
        break;
      }

      case SuperpositionType.COHERENT:
        // Coherent state - only valid for harmonic oscillator
        if (this.potentialTypeProperty.value === PotentialType.HARMONIC_OSCILLATOR) {
          const wellWidth = this.wellWidthProperty.value * QuantumConstants.NM_TO_M;
          const wellDepth = this.wellDepthProperty.value * QuantumConstants.EV_TO_JOULES;
          const mass = this.particleMassProperty.value * QuantumConstants.ELECTRON_MASS;
          const springConstant = (4 * wellDepth) / (wellWidth * wellWidth);
          const displacement = this.coherentDisplacementProperty.value * QuantumConstants.NM_TO_M;

          const result = calculateCoherentStateCoefficients(
            displacement,
            springConstant,
            mass,
            numStates
          );
          amplitudes = result.amplitudes;
          phases = result.phases;
        } else {
          // Fall back to ground state for non-harmonic potentials
          amplitudes = new Array(numStates).fill(0);
          phases = new Array(numStates).fill(0);
          amplitudes[0] = 1;
        }
        break;

      default:
        // Default to ground state
        amplitudes = new Array(numStates).fill(0);
        phases = new Array(numStates).fill(0);
        amplitudes[0] = 1;
        break;
    }

    // Update the superposition config
    this.superpositionConfigProperty.value = {
      type,
      amplitudes,
      phases,
      displacement: type === SuperpositionType.COHERENT ? this.coherentDisplacementProperty.value : undefined,
    };
  }
}
