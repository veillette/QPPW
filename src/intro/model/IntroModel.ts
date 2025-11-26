/**
 * IntroModel represents the physics model for the intro screen.
 * It extends BaseModel and provides a simplified interface for introductory exploration.
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

export type DisplayMode = "probabilityDensity" | "waveFunction";

export class IntroModel extends BaseModel {
  // ==================== CONSTANTS ====================

  /**
   * Default barrier height in electron volts.
   * Used for Rosen-Morse and Eckart potentials.
   */
  private static readonly DEFAULT_BARRIER_HEIGHT = 0.5;

  /**
   * Minimum barrier height in electron volts.
   */
  private static readonly BARRIER_HEIGHT_MIN = 0.0;

  /**
   * Maximum barrier height in electron volts.
   */
  private static readonly BARRIER_HEIGHT_MAX = 10.0;

  /**
   * Default potential offset in electron volts.
   * Used for triangular potential configuration.
   */
  private static readonly DEFAULT_POTENTIAL_OFFSET = 0.0;

  /**
   * Minimum potential offset in electron volts.
   */
  private static readonly POTENTIAL_OFFSET_MIN = -5.0;

  /**
   * Maximum potential offset in electron volts.
   */
  private static readonly POTENTIAL_OFFSET_MAX = 15.0;

  /**
   * Number of bound states to calculate.
   * Sufficient for most single-well potentials in the intro screen.
   */
  private static readonly NUM_BOUND_STATES = 10;

  /**
   * Chart display range in nanometers (extends from -RANGE to +RANGE).
   * Defines the spatial extent of the visualization.
   */
  private static readonly CHART_DISPLAY_RANGE_NM = 4;

  /**
   * Number of grid points for analytical solution evaluation.
   * High resolution ensures smooth wavefunction plots.
   */
  private static readonly ANALYTICAL_GRID_POINTS = 1000;

  /**
   * Spring constant multiplier for harmonic oscillator.
   * Factor used to convert well depth and width to spring constant: k = 8*V₀/L².
   */
  private static readonly SPRING_CONSTANT_MULTIPLIER = 8;

  /**
   * Coulomb's constant in N·m²/C².
   * Used for Coulomb potential calculations: k = 1/(4πε₀).
   */
  private static readonly COULOMB_CONSTANT = 8.9875517923e9;

  /**
   * Minimum distance for Coulomb potential calculations in meters.
   * Prevents singularity at the origin (r = 0).
   */
  private static readonly COULOMB_MIN_DISTANCE = 1e-12;

  /**
   * Half divisor for trapezoidal integration.
   * Used to calculate midpoint: (x[i+1] - x[i]) / 2.
   */
  private static readonly HALF_DIVISOR = 2;

  // ==================== PROPERTIES ====================

  // Model-specific well parameters
  public readonly barrierHeightProperty: NumberProperty;
  public readonly potentialOffsetProperty: NumberProperty;

  // Model-specific display settings (simplified, no "phaseColor" mode)
  public readonly displayModeProperty: Property<DisplayMode>;

  public constructor() {
    super();

    // Initialize model-specific well parameters
    this.barrierHeightProperty = new NumberProperty(
      IntroModel.DEFAULT_BARRIER_HEIGHT,
      {
        range: new Range(
          IntroModel.BARRIER_HEIGHT_MIN,
          IntroModel.BARRIER_HEIGHT_MAX,
        ),
      },
    );
    this.potentialOffsetProperty = new NumberProperty(
      IntroModel.DEFAULT_POTENTIAL_OFFSET,
      {
        range: new Range(
          IntroModel.POTENTIAL_OFFSET_MIN,
          IntroModel.POTENTIAL_OFFSET_MAX,
        ),
      },
    );

    // Initialize model-specific display settings (simplified)
    this.displayModeProperty = new Property<DisplayMode>("probabilityDensity");
  }

  /**
   * Setup cache invalidation listeners, including model-specific properties.
   * Extends the base implementation to add IntroModel-specific invalidation.
   */
  protected override setupCacheInvalidation(): void {
    super.setupCacheInvalidation();

    const invalidateCache = () => {
      this.boundStateResult = null;
    };

    this.barrierHeightProperty.link(invalidateCache);
    this.potentialOffsetProperty.link(invalidateCache);
  }

  /**
   * Called when the solver method or grid points changes.
   */
  protected override onSolverMethodChanged(_method: NumericalMethod): void {
    this.boundStateResult = null;
  }

  /**
   * Resets all properties to their initial state.
   * Override from BaseModel to reset model-specific properties.
   */
  public override resetAll(): void {
    super.resetAll();
    this.barrierHeightProperty.reset();
    this.potentialOffsetProperty.reset();
    this.displayModeProperty.reset();
  }

  /**
   * Calculate bound states using the Schrödinger solver.
   * Override from BaseModel.
   */
  protected override calculateBoundStates(): void {
    const wellWidth = this.wellWidthProperty.value * QuantumConstants.NM_TO_M;
    const wellDepth =
      this.wellDepthProperty.value * QuantumConstants.EV_TO_JOULES;
    const mass =
      this.particleMassProperty.value * QuantumConstants.ELECTRON_MASS;

    const numStates = IntroModel.NUM_BOUND_STATES;

    const gridConfig = {
      xMin: -IntroModel.CHART_DISPLAY_RANGE_NM * QuantumConstants.NM_TO_M,
      xMax: IntroModel.CHART_DISPLAY_RANGE_NM * QuantumConstants.NM_TO_M,
      numPoints: IntroModel.ANALYTICAL_GRID_POINTS,
    };

    try {
      const potentialParams: WellParameters = {
        type: this.potentialTypeProperty.value,
        wellWidth: wellWidth,
      };

      // Add type-specific parameters
      switch (this.potentialTypeProperty.value) {
        case PotentialType.FINITE_WELL:
          potentialParams.wellDepth = wellDepth;
          break;
        case PotentialType.HARMONIC_OSCILLATOR:
          potentialParams.springConstant =
            (IntroModel.SPRING_CONSTANT_MULTIPLIER * wellDepth) /
            (wellWidth * wellWidth);
          break;
        case PotentialType.MORSE:
          potentialParams.dissociationEnergy = wellDepth;
          potentialParams.equilibriumPosition = 0;
          potentialParams.wellWidth = wellWidth;
          break;
        case PotentialType.POSCHL_TELLER:
          potentialParams.potentialDepth = wellDepth;
          potentialParams.wellWidth = wellWidth;
          break;
        case PotentialType.ROSEN_MORSE:
          potentialParams.potentialDepth = wellDepth;
          potentialParams.barrierHeight =
            this.barrierHeightProperty.value * QuantumConstants.EV_TO_JOULES;
          potentialParams.wellWidth = wellWidth;
          break;
        case PotentialType.ECKART:
          potentialParams.potentialDepth = wellDepth;
          potentialParams.barrierHeight =
            this.barrierHeightProperty.value * QuantumConstants.EV_TO_JOULES;
          potentialParams.wellWidth = wellWidth;
          break;
        case PotentialType.ASYMMETRIC_TRIANGLE:
          potentialParams.slope = wellDepth / wellWidth;
          break;
        case PotentialType.TRIANGULAR:
          potentialParams.potentialDepth = wellDepth;
          potentialParams.wellWidth = wellWidth;
          potentialParams.energyOffset =
            this.potentialOffsetProperty.value * QuantumConstants.EV_TO_JOULES;
          break;
        case PotentialType.COULOMB_1D:
        case PotentialType.COULOMB_3D: {
          potentialParams.coulombStrength =
            IntroModel.COULOMB_CONSTANT *
            QuantumConstants.ELEMENTARY_CHARGE *
            QuantumConstants.ELEMENTARY_CHARGE;
          break;
        }
      }

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
   */
  public override getClassicalProbabilityDensity(
    energyIndex: number,
  ): number[] | null {
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

    // Fallback: numerical calculation using potential function
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
   * Calculates the classical turning points for a given energy level.
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

    const energy = this.boundStateResult.energies[energyLevel];
    const energyEV = energy * QuantumConstants.JOULES_TO_EV;
    const xGrid = this.boundStateResult.xGrid;

    let leftTurningPoint: number | null = null;
    let rightTurningPoint: number | null = null;

    // Handle infinite well directly (classical turning points at well edges)
    if (this.potentialTypeProperty.value === PotentialType.INFINITE_WELL) {
      const halfWidth = this.wellWidthProperty.value / 2;
      return { left: -halfWidth, right: halfWidth };
    }

    // Existing loop to find turning points for other potentials
    for (let i = 0; i < xGrid.length - 1; i++) {
      const x = xGrid[i] * QuantumConstants.M_TO_NM;
      const V = this.getPotentialAtPosition(x);

      const xNext = xGrid[i + 1] * QuantumConstants.M_TO_NM;
      const VNext = this.getPotentialAtPosition(xNext);

      if (
        (V <= energyEV && VNext >= energyEV) ||
        (V >= energyEV && VNext <= energyEV)
      ) {
        // Avoid division by zero when V and VNext are equal
        if (VNext === V) {
          continue; // skip this segment
        }
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
      // Clamp turning points to chart's X-axis range
      const minX = -IntroModel.CHART_DISPLAY_RANGE_NM;
      const maxX = IntroModel.CHART_DISPLAY_RANGE_NM;
      leftTurningPoint = Math.max(minX, Math.min(maxX, leftTurningPoint));
      rightTurningPoint = Math.max(minX, Math.min(maxX, rightTurningPoint));
      return { left: leftTurningPoint, right: rightTurningPoint };
    }

    return null;
  }

  /**
   * Calculates the probability of finding the particle in the classically forbidden region.
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
      const x = xGrid[i] * QuantumConstants.M_TO_NM;
      const psi = wavefunction[i];
      const probabilityDensity = psi * psi;

      let dx = 0;
      if (i === 0) {
        dx =
          ((xGrid[1] - xGrid[0]) / IntroModel.HALF_DIVISOR) *
          QuantumConstants.M_TO_NM;
      } else if (i === xGrid.length - 1) {
        dx =
          ((xGrid[i] - xGrid[i - 1]) / IntroModel.HALF_DIVISOR) *
          QuantumConstants.M_TO_NM;
      } else {
        dx =
          ((xGrid[i + 1] - xGrid[i - 1]) / IntroModel.HALF_DIVISOR) *
          QuantumConstants.M_TO_NM;
      }

      totalProbability += probabilityDensity * dx;

      if (x < turningPoints.left || x > turningPoints.right) {
        forbiddenProbability += probabilityDensity * dx;
      }
    }

    return totalProbability > 0
      ? (forbiddenProbability / totalProbability) * 100
      : 0;
  }

  /**
   * Calculate the potential energy at given positions.
   */
  private calculatePotentialEnergy(xGrid: number[]): number[] {
    const wellWidth = this.wellWidthProperty.value * QuantumConstants.NM_TO_M;
    const wellDepth =
      this.wellDepthProperty.value * QuantumConstants.EV_TO_JOULES;

    const potential: number[] = [];

    for (let i = 0; i < xGrid.length; i++) {
      const x = xGrid[i];
      let V = 0;

      switch (this.potentialTypeProperty.value) {
        case PotentialType.INFINITE_WELL:
          V = Math.abs(x) <= wellWidth / 2 ? 0 : Infinity;
          break;

        case PotentialType.FINITE_WELL:
          V = Math.abs(x) <= wellWidth / 2 ? 0 : wellDepth;
          break;

        case PotentialType.HARMONIC_OSCILLATOR: {
          const springConstant =
            (IntroModel.SPRING_CONSTANT_MULTIPLIER * wellDepth) /
            (wellWidth * wellWidth);
          V = 0.5 * springConstant * x * x;
          break;
        }

        case PotentialType.MORSE: {
          const a = 1 / wellWidth;
          const exponential = Math.exp(-a * x);
          V = wellDepth * Math.pow(1 - exponential, 2);
          break;
        }

        case PotentialType.POSCHL_TELLER: {
          const coshPT = Math.cosh(x / wellWidth);
          V = -wellDepth / (coshPT * coshPT);
          break;
        }

        case PotentialType.ROSEN_MORSE: {
          const barrierHeight =
            this.barrierHeightProperty.value * QuantumConstants.EV_TO_JOULES;
          const coshRM = Math.cosh(x / wellWidth);
          const tanhRM = Math.tanh(x / wellWidth);
          V = -wellDepth / (coshRM * coshRM) + barrierHeight * tanhRM;
          break;
        }

        case PotentialType.ECKART: {
          const barrierHeightE =
            this.barrierHeightProperty.value * QuantumConstants.EV_TO_JOULES;
          const expE = Math.exp(x / wellWidth);
          const denomE = 1 + expE;
          V = wellDepth / (denomE * denomE) - barrierHeightE / denomE;
          break;
        }

        case PotentialType.ASYMMETRIC_TRIANGLE: {
          const slope = wellDepth / wellWidth;
          V = slope * x;
          break;
        }

        case PotentialType.TRIANGULAR: {
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
          const coulombStrength =
            IntroModel.COULOMB_CONSTANT *
            QuantumConstants.ELEMENTARY_CHARGE *
            QuantumConstants.ELEMENTARY_CHARGE;
          const r = Math.abs(x);
          if (r > IntroModel.COULOMB_MIN_DISTANCE) {
            V = -coulombStrength / r;
          } else {
            V = -coulombStrength / IntroModel.COULOMB_MIN_DISTANCE;
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
   * Get the potential energy at a specific position (in nm).
   */
  private getPotentialAtPosition(xNm: number): number {
    const x = xNm * QuantumConstants.NM_TO_M;
    const potential = this.calculatePotentialEnergy([x]);
    return potential[0] * QuantumConstants.JOULES_TO_EV;
  }
}
