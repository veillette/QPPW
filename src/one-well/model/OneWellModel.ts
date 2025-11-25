/**
 * OneWellModel represents the physics model for a single quantum potential well.
 * It handles the quantum mechanical calculations for a particle in a single well.
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
import {
  calculateCoherentStateCoefficients,
  createHarmonicOscillatorPotential,
  calculateHarmonicOscillatorClassicalProbability,
} from "../../common/model/analytical-solutions/harmonic-oscillator.js";
import {
  createInfiniteWellPotential,
  calculateInfiniteWellClassicalProbability,
} from "../../common/model/analytical-solutions/infinite-square-well.js";
import {
  createFiniteWellPotential,
  calculateFiniteWellClassicalProbability,
} from "../../common/model/analytical-solutions/finite-square-well.js";

export type DisplayMode = "probabilityDensity" | "waveFunction" | "phaseColor";

export class OneWellModel extends BaseModel {
  // Potential type selection
  public readonly potentialTypeProperty: Property<PotentialType>;

  // Well parameters (used for different potential types)
  public readonly wellWidthProperty: NumberProperty;
  public readonly wellDepthProperty: NumberProperty;
  public readonly wellOffsetProperty: NumberProperty; // For asymmetric wells
  public readonly barrierHeightProperty: NumberProperty; // For Rosen-Morse and Eckart potentials
  public readonly potentialOffsetProperty: NumberProperty; // For triangular potential

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

  // Coherent state parameter
  public readonly coherentDisplacementProperty: NumberProperty; // Displacement in nm

  // Cached bound state results
  protected boundStateResult: BoundStateResult | null = null;

  public constructor() {
    super();
    // Initialize potential type (square/infinite well by default)
    this.potentialTypeProperty = new Property<PotentialType>(
      PotentialType.INFINITE_WELL,
    );

    // Initialize well parameters with default values
    this.wellWidthProperty = new NumberProperty(1.0, {
      range: new Range(0.1, 6.0),
    }); // in nanometers (max 6 nm)
    this.wellDepthProperty = new NumberProperty(5.0, {
      range: new Range(0.1, 15.0),
    }); // in eV (within energy graph bounds)
    this.wellOffsetProperty = new NumberProperty(0.5, {
      range: new Range(0.0, 1.0),
    }); // normalized position
    this.barrierHeightProperty = new NumberProperty(0.5, {
      range: new Range(0.0, 10.0),
    }); // in eV (for Rosen-Morse and Eckart)
    this.potentialOffsetProperty = new NumberProperty(0.0, {
      range: new Range(-5.0, 15.0),
    }); // in eV (for triangular potential)

    // Initialize particle mass (1.0 = electron mass)
    this.particleMassProperty = new NumberProperty(1.0, {
      range: new Range(0.5, 1.1),
    }); // 0.5 to 1.1 times electron mass

    // Initialize energy level selection (ground state by default)
    // Use a large range to accommodate potentials with many states (up to 100)
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
    // Default to PSI_K (single eigenstate) so selecting energy levels works properly
    this.superpositionTypeProperty = new Property<SuperpositionType>(
      SuperpositionType.SINGLE,
    );
    this.superpositionConfigProperty = new Property<SuperpositionConfig>({
      type: SuperpositionType.SINGLE,
      amplitudes: [1.0], // Default to ground state
      phases: [0],
    });

    // Initialize coherent state displacement (0.5 nm default, range 0 to 2 nm)
    this.coherentDisplacementProperty = new NumberProperty(0.5, {
      range: new Range(0.0, 2.0),
    });

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
    this.potentialOffsetProperty.link(invalidateCache);
    this.particleMassProperty.link(invalidateCache);

    // Update superposition coefficients when superposition type or coherent displacement changes
    this.superpositionTypeProperty.link(() =>
      this.updateSuperpositionCoefficients(),
    );
    this.coherentDisplacementProperty.link(() => {
      if (this.superpositionTypeProperty.value === SuperpositionType.COHERENT) {
        this.updateSuperpositionCoefficients();
      }
    });
  }

  /**
   * Called when the solver method or grid points changes.
   * For OneWellModel, solver method changes don't matter since we always use analytical solutions,
   * but grid points changes do matter because they affect the wavefunction grid.
   * @param _method - The new numerical method (unused)
   */

  protected onSolverMethodChanged(_method: NumericalMethod): void {
    // Invalidate the cached bound state result so wavefunctions are recalculated
    // with the new grid configuration
    this.boundStateResult = null;
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
    this.potentialOffsetProperty.reset();
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
    this.coherentDisplacementProperty.reset();
  }

  /**
   * Calculates the classical turning points for a given energy level.
   * Returns the x positions (in nm) where the potential energy equals the eigenstate energy.
   * Returns null if turning points cannot be calculated.
   */
  public getClassicalTurningPoints(energyLevel: number): {
    left: number;
    right: number;
  } | null {
    if (!this.boundStateResult) {
      this.calculateBoundStates();
    }

    if (
      !this.boundStateResult ||
      energyLevel < 0 ||
      energyLevel >= this.boundStateResult.energies.length
    ) {
      return null;
    }

    const energy = this.boundStateResult.energies[energyLevel]; // in Joules
    const energyEV = energy * QuantumConstants.JOULES_TO_EV; // in eV
    const xGrid = this.boundStateResult.xGrid; // in meters

    // Find the turning points by finding where V(x) ≈ E
    let leftTurningPoint: number | null = null;
    let rightTurningPoint: number | null = null;

    for (let i = 0; i < xGrid.length - 1; i++) {
      const x = xGrid[i] * QuantumConstants.M_TO_NM; // Convert to nm
      const V = this.getPotentialAtPosition(x); // in eV

      const xNext = xGrid[i + 1] * QuantumConstants.M_TO_NM;
      const VNext = this.getPotentialAtPosition(xNext);

      // Check if energy crosses potential between these two points
      if (
        (V <= energyEV && VNext >= energyEV) ||
        (V >= energyEV && VNext <= energyEV)
      ) {
        // Linear interpolation to find crossing point
        const t = (energyEV - V) / (VNext - V);
        const turningPoint = x + t * (xNext - x);

        if (leftTurningPoint === null) {
          leftTurningPoint = turningPoint;
        } else {
          rightTurningPoint = turningPoint;
        }
      }
    }

    if (leftTurningPoint !== null && rightTurningPoint !== null) {
      return { left: leftTurningPoint, right: rightTurningPoint };
    }

    return null;
  }

  /**
   * Helper method to get potential energy at a given position (in nm).
   * Returns potential in eV.
   * Uses the analytical solver's potential function when available.
   */
  private getPotentialAtPosition(xNm: number): number {
    const potentialType = this.potentialTypeProperty.value;
    const wellWidth = this.wellWidthProperty.value * QuantumConstants.NM_TO_M;
    const wellDepth = this.wellDepthProperty.value;
    const xMeters = xNm * QuantumConstants.NM_TO_M;

    // Try to use analytical solver potential functions
    let potentialFunction: ((x: number) => number) | null = null;

    switch (potentialType) {
      case PotentialType.INFINITE_WELL:
        potentialFunction = createInfiniteWellPotential(wellWidth);
        break;

      case PotentialType.FINITE_WELL:
        potentialFunction = createFiniteWellPotential(
          wellWidth,
          wellDepth * QuantumConstants.EV_TO_JOULES,
        );
        break;

      case PotentialType.HARMONIC_OSCILLATOR: {
        const wellDepthJoules = wellDepth * QuantumConstants.EV_TO_JOULES;
        const springConstant = (8 * wellDepthJoules) / (wellWidth * wellWidth);
        potentialFunction = createHarmonicOscillatorPotential(springConstant);
        break;
      }

      default:
        // Fall through to legacy implementation for other potentials
        break;
    }

    // If we have an analytical potential function, use it
    if (potentialFunction) {
      const VJoules = potentialFunction(xMeters);
      return VJoules * QuantumConstants.JOULES_TO_EV;
    }

    // Legacy implementation for potentials without analytical solver support
    if (potentialType === PotentialType.ASYMMETRIC_TRIANGLE) {
      if (xNm < 0) {
        return 1000; // Infinite wall
      } else {
        const F_eV_per_nm = wellDepth / (wellWidth * QuantumConstants.M_TO_NM);
        return F_eV_per_nm * xNm;
      }
    } else if (potentialType === PotentialType.TRIANGULAR) {
      const offset = this.potentialOffsetProperty.value;
      const height = wellDepth;
      const widthNm = wellWidth * QuantumConstants.M_TO_NM;
      if (xNm < 0) {
        return height + offset;
      } else if (xNm < widthNm) {
        return offset + (height / widthNm) * xNm;
      } else {
        return height + offset;
      }
    } else if (
      potentialType === PotentialType.COULOMB_1D ||
      potentialType === PotentialType.COULOMB_3D
    ) {
      const widthNm = wellWidth * QuantumConstants.M_TO_NM;
      const k = wellDepth * (widthNm / 2);
      const distance = Math.max(Math.abs(xNm), 0.01);
      return -k / distance;
    } else if (potentialType === PotentialType.MORSE) {
      const widthNm = wellWidth * QuantumConstants.M_TO_NM;
      const exponent = Math.exp(-xNm / widthNm);
      return wellDepth * Math.pow(1 - exponent, 2) - wellDepth;
    } else if (potentialType === PotentialType.POSCHL_TELLER) {
      const widthNm = wellWidth * QuantumConstants.M_TO_NM;
      const coshVal = Math.cosh(xNm / widthNm);
      return -wellDepth / (coshVal * coshVal);
    } else if (potentialType === PotentialType.ROSEN_MORSE) {
      const widthNm = wellWidth * QuantumConstants.M_TO_NM;
      const barrierHeight = this.barrierHeightProperty.value;
      const coshVal = Math.cosh(xNm / widthNm);
      const tanhVal = Math.tanh(xNm / widthNm);
      return -wellDepth / (coshVal * coshVal) + barrierHeight * tanhVal;
    } else if (potentialType === PotentialType.ECKART) {
      const widthNm = wellWidth * QuantumConstants.M_TO_NM;
      const barrierHeight = this.barrierHeightProperty.value;
      const expVal = Math.exp(xNm / widthNm);
      const denom = 1 + expVal;
      return wellDepth / (denom * denom) - barrierHeight / denom;
    }

    // Default case
    return 0;
  }

  /**
   * Calculates the probability of finding the particle in the classically forbidden region.
   * Returns a percentage (0-100).
   */
  public getClassicallyForbiddenProbability(energyLevel: number): number {
    if (!this.boundStateResult) {
      this.calculateBoundStates();
    }

    if (
      !this.boundStateResult ||
      energyLevel < 0 ||
      energyLevel >= this.boundStateResult.energies.length
    ) {
      return 0;
    }

    const turningPoints = this.getClassicalTurningPoints(energyLevel);
    if (!turningPoints) {
      return 0;
    }

    const wavefunction = this.boundStateResult.wavefunctions[energyLevel];
    const xGrid = this.boundStateResult.xGrid;

    let forbiddenProbability = 0;
    let totalProbability = 0;

    for (let i = 0; i < xGrid.length; i++) {
      const x = xGrid[i] * QuantumConstants.M_TO_NM; // Convert to nm
      const psi = wavefunction[i];
      const probabilityDensity = psi * psi;

      // Calculate dx for integration (using trapezoidal rule)
      let dx = 0;
      if (i === 0) {
        dx =
          ((xGrid[1] - xGrid[0]) / 2) * QuantumConstants.M_TO_NM;
      } else if (i === xGrid.length - 1) {
        dx =
          ((xGrid[i] - xGrid[i - 1]) / 2) * QuantumConstants.M_TO_NM;
      } else {
        dx =
          ((xGrid[i + 1] - xGrid[i - 1]) / 2) * QuantumConstants.M_TO_NM;
      }

      totalProbability += probabilityDensity * dx;

      // Check if we're in the classically forbidden region
      if (x < turningPoints.left || x > turningPoints.right) {
        forbiddenProbability += probabilityDensity * dx;
      }
    }

    // Return as percentage
    return totalProbability > 0
      ? (forbiddenProbability / totalProbability) * 100
      : 0;
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
    if (
      this.boundStateResult &&
      n > 0 &&
      n <= this.boundStateResult.energies.length
    ) {
      const energyJoules = this.boundStateResult.energies[n - 1];
      return Schrodinger1DSolver.joulesToEV(energyJoules);
    }

    // Fallback to analytical formula if solver fails
    const L = this.wellWidthProperty.value * QuantumConstants.NM_TO_M;
    const energy =
      (n *
        n *
        Math.PI *
        Math.PI *
        QuantumConstants.HBAR *
        QuantumConstants.HBAR) /
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
    // All potentials now use wellWidth as a width parameter in nanometers (converted to meters)
    const wellWidth = this.wellWidthProperty.value * QuantumConstants.NM_TO_M; // width in meters
    const wellDepth =
      this.wellDepthProperty.value * QuantumConstants.EV_TO_JOULES;
    const mass =
      this.particleMassProperty.value * QuantumConstants.ELECTRON_MASS;

    // Calculate number of states based on potential type and energy range
    let numStates = 10; // Default for most potentials

    // For harmonic oscillator, calculate states up to 15 eV
    if (
      this.potentialTypeProperty.value === PotentialType.HARMONIC_OSCILLATOR
    ) {
      const springConstant = (8 * wellDepth) / (wellWidth * wellWidth);
      const omega = Math.sqrt(springConstant / mass);
      const maxEnergy = 15 * QuantumConstants.EV_TO_JOULES; // 15 eV
      // E_n = ℏω(n + 1/2), solve for n: n = E/(ℏω) - 1/2
      const maxN = Math.floor(
        maxEnergy / (QuantumConstants.HBAR * omega) - 0.5,
      );
      numStates = Math.max(1, Math.min(maxN + 1, 100)); // Cap at 100 for safety
    }
    // For infinite well, calculate states up to 15 eV
    else if (this.potentialTypeProperty.value === PotentialType.INFINITE_WELL) {
      const maxEnergy = 15 * QuantumConstants.EV_TO_JOULES; // 15 eV
      // E_n = (ℏ²π²n²)/(2mL²), solve for n
      const maxN = Math.floor(
        Math.sqrt(
          (2 * mass * wellWidth * wellWidth * maxEnergy) /
          (QuantumConstants.HBAR * QuantumConstants.HBAR * Math.PI * Math.PI),
        ),
      );
      numStates = Math.max(1, Math.min(maxN, 100)); // Cap at 100 for safety
    }
    // For finite well, estimate maximum number of bound states
    else if (this.potentialTypeProperty.value === PotentialType.FINITE_WELL) {
      // Approximate number of bound states: n_max ≈ (1/π) * sqrt(2mV₀L²/ℏ²)
      // Use generous estimate to ensure we get all states
      const estimatedMax = Math.ceil(
        (1 / Math.PI) *
        Math.sqrt(
          (2 * mass * wellDepth * wellWidth * wellWidth) /
          (QuantumConstants.HBAR * QuantumConstants.HBAR),
        ),
      );
      // Request more states than estimated to ensure we capture all bound states
      numStates = Math.max(10, Math.min(estimatedMax * 2, 100)); // At least 10, cap at 100
    }
    // For asymmetric triangle, calculate states that fit in the energy range
    else if (
      this.potentialTypeProperty.value === PotentialType.ASYMMETRIC_TRIANGLE
    ) {
      numStates = 80; // Asymmetric triangle may have many states, use larger number
    }
    // For triangular potential, calculate states based on well depth
    else if (this.potentialTypeProperty.value === PotentialType.TRIANGULAR) {
      numStates = 50; // Triangular well typically has fewer states than asymmetric
    } else {
      numStates = 80; // Use more states for other potentials
    }

    // Grid configuration spans the full chart display range (-4 nm to +4 nm)
    // This ensures wavefunctions are calculated across the entire visible area
    const CHART_DISPLAY_RANGE_NM = 4; // Match the chart x-axis range

    // For analytical solutions, use a fixed high-resolution grid (1000 points)
    // since we can evaluate them exactly at any point
    const ANALYTICAL_GRID_POINTS = 1000;

    const gridConfig = {
      xMin: -CHART_DISPLAY_RANGE_NM * QuantumConstants.NM_TO_M,
      xMax: CHART_DISPLAY_RANGE_NM * QuantumConstants.NM_TO_M,
      numPoints: ANALYTICAL_GRID_POINTS, // Use fixed high-resolution grid for analytical solutions
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
          potentialParams.springConstant =
            (8 * wellDepth) / (wellWidth * wellWidth);
          break;
        case PotentialType.MORSE:
          // Morse potential: V(x) = D_e * (1 - exp(-(x - x_e)/a))^2
          // wellWidth is the width parameter 'a' in meters
          potentialParams.dissociationEnergy = wellDepth;
          potentialParams.equilibriumPosition = 0; // Center of the well
          potentialParams.wellWidth = wellWidth; // width in meters
          break;
        case PotentialType.POSCHL_TELLER:
          // Pöschl-Teller potential: V(x) = -V_0 / cosh²(x/a)
          // wellWidth is the width parameter 'a' in meters
          potentialParams.potentialDepth = wellDepth;
          potentialParams.wellWidth = wellWidth; // width in meters
          break;
        case PotentialType.ROSEN_MORSE:
          // Rosen-Morse potential: V(x) = -V_0 / cosh²(x/a) + V_1 * tanh(x/a)
          // wellWidth is the width parameter 'a' in meters
          potentialParams.potentialDepth = wellDepth;
          potentialParams.barrierHeight =
            this.barrierHeightProperty.value * QuantumConstants.EV_TO_JOULES;
          potentialParams.wellWidth = wellWidth; // width in meters
          break;
        case PotentialType.ECKART:
          // Eckart potential: V(x) = V_0 / (1 + exp(x/a))² - V_1 / (1 + exp(x/a))
          // wellWidth is the width parameter 'a' in meters
          potentialParams.potentialDepth = wellDepth;
          potentialParams.barrierHeight =
            this.barrierHeightProperty.value * QuantumConstants.EV_TO_JOULES;
          potentialParams.wellWidth = wellWidth; // width in meters
          break;
        case PotentialType.ASYMMETRIC_TRIANGLE:
          // Slope is the field strength
          potentialParams.slope = wellDepth / wellWidth;
          break;
        case PotentialType.TRIANGULAR:
          // Triangular potential:
          // V(x) = height + offset for x < 0
          // V(x) = offset at x = 0
          // V(x) = offset + (height/width) * x for 0 < x < width
          // V(x) = height + offset for x > width
          potentialParams.potentialDepth = wellDepth; // height in Joules
          potentialParams.wellWidth = wellWidth; // width in meters
          potentialParams.energyOffset =
            this.potentialOffsetProperty.value * QuantumConstants.EV_TO_JOULES; // offset in Joules
          break;
        case PotentialType.COULOMB_1D:
        case PotentialType.COULOMB_3D: {
          // For Coulomb potentials, use coulombStrength parameter α = k*e²
          // where k = 1/(4πε₀) ≈ 8.9875517923e9 N·m²/C²
          // α ≈ 2.307e-28 J·m for electron charge
          // Energy then scales naturally with mass: E_n = -mα²/(2ℏ²n²)
          // With electron mass, this gives E_1 = -13.6 eV
          const coulombConstant = 8.9875517923e9; // Coulomb's constant in N·m²/C²
          potentialParams.coulombStrength =
            coulombConstant *
            QuantumConstants.ELEMENTARY_CHARGE *
            QuantumConstants.ELEMENTARY_CHARGE;
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
   * For many potentials, the analytical solver provides an exact normalized solution.
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
    const wellWidth = this.wellWidthProperty.value * QuantumConstants.NM_TO_M;
    const wellDepth =
      this.wellDepthProperty.value * QuantumConstants.EV_TO_JOULES;

    // Use analytical solver methods when available
    const potentialType = this.potentialTypeProperty.value;

    try {
      switch (potentialType) {
        case PotentialType.INFINITE_WELL:
          return calculateInfiniteWellClassicalProbability(
            wellWidth,
            energy,
            mass,
            xGrid,
          );

        case PotentialType.FINITE_WELL:
          return calculateFiniteWellClassicalProbability(
            wellWidth,
            wellDepth,
            energy,
            mass,
            xGrid,
          );

        case PotentialType.HARMONIC_OSCILLATOR: {
          const springConstant = (8 * wellDepth) / (wellWidth * wellWidth);
          return calculateHarmonicOscillatorClassicalProbability(
            springConstant,
            energy,
            mass,
            xGrid,
          );
        }

        default:
          // Fall back to numerical calculation for other potentials
          break;
      }
    } catch (error) {
      console.warn(
        `Failed to use analytical classical probability for ${potentialType}, falling back to numerical:`,
        error,
      );
    }

    // Fallback: numerical calculation using potential function
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
        // Add epsilon to avoid division by zero at turning points
        const epsilon = 1e-10 * QuantumConstants.EV_TO_JOULES;
        const probability =
          1 / Math.sqrt((2 * Math.max(kineticEnergy, epsilon)) / mass);
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
   * Uses the analytical solver's potential function when available.
   * @param xGrid - Array of x positions in meters
   * @returns Array of potential energy values in Joules
   */
  private calculatePotentialEnergy(xGrid: number[]): number[] {
    const wellWidth = this.wellWidthProperty.value * QuantumConstants.NM_TO_M;
    const wellDepth =
      this.wellDepthProperty.value * QuantumConstants.EV_TO_JOULES;

    // Try to use analytical solver potential functions
    let potentialFunction: ((x: number) => number) | null = null;

    switch (this.potentialTypeProperty.value) {
      case PotentialType.INFINITE_WELL:
        potentialFunction = createInfiniteWellPotential(wellWidth);
        break;

      case PotentialType.FINITE_WELL:
        potentialFunction = createFiniteWellPotential(wellWidth, wellDepth);
        break;

      case PotentialType.HARMONIC_OSCILLATOR: {
        const springConstant = (8 * wellDepth) / (wellWidth * wellWidth);
        potentialFunction = createHarmonicOscillatorPotential(springConstant);
        break;
      }

      default:
        // Fall through to legacy implementation for other potentials
        break;
    }

    // If we have an analytical potential function, use it
    if (potentialFunction) {
      return xGrid.map((x) => potentialFunction!(x));
    }

    // Legacy implementation for potentials without analytical solver support
    const potential: number[] = [];

    for (let i = 0; i < xGrid.length; i++) {
      const x = xGrid[i];
      let V = 0;

      switch (this.potentialTypeProperty.value) {
        case PotentialType.MORSE: {
          // V(x) = D_e * (1 - exp(-a(x - x_e)))^2
          const a = 1 / wellWidth;
          const exponential = Math.exp(-a * x);
          V = wellDepth * Math.pow(1 - exponential, 2);
          break;
        }

        case PotentialType.POSCHL_TELLER: {
          // V(x) = -V_0 / cosh^2(x/a)
          const coshPT = Math.cosh(x / wellWidth);
          V = -wellDepth / (coshPT * coshPT);
          break;
        }

        case PotentialType.ROSEN_MORSE: {
          // V(x) = -V_0 / cosh^2(x/a) + V_1 * tanh(x/a)
          const barrierHeight =
            this.barrierHeightProperty.value * QuantumConstants.EV_TO_JOULES;
          const coshRM = Math.cosh(x / wellWidth);
          const tanhRM = Math.tanh(x / wellWidth);
          V = -wellDepth / (coshRM * coshRM) + barrierHeight * tanhRM;
          break;
        }

        case PotentialType.ECKART: {
          // V(x) = V_0 / (1 + exp(x/a))^2 - V_1 / (1 + exp(x/a))
          const barrierHeightE =
            this.barrierHeightProperty.value * QuantumConstants.EV_TO_JOULES;
          const expE = Math.exp(x / wellWidth);
          const denomE = 1 + expE;
          V = wellDepth / (denomE * denomE) - barrierHeightE / denomE;
          break;
        }

        case PotentialType.ASYMMETRIC_TRIANGLE: {
          // V(x) = F * x (linear potential with electric field)
          const slope = wellDepth / wellWidth;
          V = slope * x;
          break;
        }

        case PotentialType.TRIANGULAR: {
          // Triangular potential
          const energyOffset =
            this.potentialOffsetProperty.value * QuantumConstants.EV_TO_JOULES;
          if (x < 0) {
            V = wellDepth + energyOffset;
          } else if (x < wellWidth) {
            V = energyOffset + (wellDepth / wellWidth) * x;
          } else {
            V = wellDepth + energyOffset;
          }
          break;
        }

        case PotentialType.COULOMB_1D:
        case PotentialType.COULOMB_3D: {
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

  /**
   * Get the superposition wavefunction based on current superposition config.
   * Returns the real-space wavefunction: ψ(x) = Σ c_n * e^(iφ_n) * ψ_n(x)
   * Since we're returning real values for display, this returns the real part at t=0.
   *
   * @returns Object with wavefunction array and energy (weighted average), or null
   */
  public getSuperpositionWavefunction(): {
    wavefunction: number[];
    energy: number;
  } | null {
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
        weightedEnergy +=
          amplitude * amplitude * this.boundStateResult.energies[n];
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

      case SuperpositionType.SINGLE:
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
        const normNarrow = Math.sqrt(
          amplitudes.reduce((sum, a) => sum + a * a, 0),
        );
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
        const normWide = Math.sqrt(
          amplitudes.reduce((sum, a) => sum + a * a, 0),
        );
        amplitudes = amplitudes.map((a) => a / normWide);
        break;
      }

      case SuperpositionType.COHERENT: {
        const wellWidth =
          this.wellWidthProperty.value * QuantumConstants.NM_TO_M;
        const displacement =
          this.coherentDisplacementProperty.value * QuantumConstants.NM_TO_M;

        if (
          this.potentialTypeProperty.value === PotentialType.HARMONIC_OSCILLATOR
        ) {
          // True coherent state for harmonic oscillator
          const wellDepth =
            this.wellDepthProperty.value * QuantumConstants.EV_TO_JOULES;
          const mass =
            this.particleMassProperty.value * QuantumConstants.ELECTRON_MASS;
          const springConstant = (4 * wellDepth) / (wellWidth * wellWidth);

          const result = calculateCoherentStateCoefficients(
            displacement,
            springConstant,
            mass,
            numStates,
          );
          amplitudes = result.amplitudes;
          phases = result.phases;
        } else {
          // Coherent-state-like wavepacket for other potentials
          // Create a localized superposition centered at the displacement position
          // using a Gaussian envelope in eigenstate space
          amplitudes = new Array(numStates).fill(0);
          phases = new Array(numStates).fill(0);

          // For each eigenstate, compute its overlap with a position-localized state
          // We use the eigenfunction values at the displacement position as weights
          const displacementIndex = Math.round(
            ((displacement + 4 * QuantumConstants.NM_TO_M) /
              (8 * QuantumConstants.NM_TO_M)) *
            (this.boundStateResult!.xGrid.length - 1),
          );
          const clampedIndex = Math.max(
            0,
            Math.min(
              displacementIndex,
              this.boundStateResult!.xGrid.length - 1,
            ),
          );

          // Weight each eigenstate by its wavefunction value at the displacement position
          // and apply a Gaussian envelope in eigenstate index to create localization
          const sigma = 3.0; // Width of Gaussian in eigenstate space
          const n0 = Math.min(5, Math.floor(numStates / 2)); // Center around a middle eigenstate

          for (let n = 0; n < numStates; n++) {
            if (n < this.boundStateResult!.wavefunctions.length) {
              const psi_n_at_x =
                this.boundStateResult!.wavefunctions[n][clampedIndex];

              // Gaussian envelope in eigenstate space
              const gaussianWeight = Math.exp(
                -Math.pow(n - n0, 2) / (2 * sigma * sigma),
              );

              // Combine position-based and Gaussian weights
              amplitudes[n] = psi_n_at_x * gaussianWeight;
            }
          }

          // Normalize the amplitudes
          const norm = Math.sqrt(amplitudes.reduce((sum, a) => sum + a * a, 0));
          if (norm > 0) {
            amplitudes = amplitudes.map((a) => a / norm);
          } else {
            // Fallback: if all amplitudes are zero, just use ground state
            amplitudes[0] = 1;
          }

          // Set phases based on momentum (for a wavepacket moving toward the center)
          // Phase increases linearly with eigenstate index to create momentum
          const momentumDirection = displacement > 0 ? -1 : 1; // Move toward center
          for (let n = 0; n < numStates; n++) {
            phases[n] = momentumDirection * n * 0.1; // Small phase gradient
          }
        }
        break;
      }

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
      displacement:
        type === SuperpositionType.COHERENT
          ? this.coherentDisplacementProperty.value
          : undefined,
    };
  }
}
