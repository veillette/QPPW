/**
 * Main solver for the 1D time-independent Schrödinger equation.
 * Provides a unified interface for analytical and numerical solutions.
 *
 * Usage:
 *   const solver = new Schrodinger1DSolver();
 *   const result = solver.solve(potential, mass, numStates, gridConfig);
 */

import {
  BoundStateResult,
  GridConfig,
  PotentialFunction,
  PotentialType,
} from "./PotentialFunction.js";
import {
  AnalyticalSolution,
  InfiniteSquareWellSolution,
  FiniteSquareWellSolution,
  HarmonicOscillatorSolution,
  MorsePotentialSolution,
  PoschlTellerPotentialSolution,
  RosenMorsePotentialSolution,
  EckartPotentialSolution,
  AsymmetricTrianglePotentialSolution,
  Coulomb1DPotentialSolution,
  Coulomb3DPotentialSolution,
  TriangularPotentialSolution,
  solveDoubleSquareWellAnalytical,
} from "./analytical-solutions";
import { solveMultiSquareWell } from "./analytical-solutions/multi-square-well.js";
import { solveMultiCoulomb1D } from "./analytical-solutions/multi-coulomb-1d.js";
import { solveNumerov } from "./NumerovSolver.js";
import { solveMatrixNumerov } from "./MatrixNumerovSolver.js";
import { solveDVR } from "./DVRSolver.js";
import { solveFGH } from "./FGHSolver.js";
import { solveSpectral } from "./SpectralSolver.js";
import { solveQuantumBound } from "./QuantumBoundStateSolver.js";
import { computeWavefunctionsNumerov } from "./WavefunctionNumerovSolver.js";
import QuantumConstants from "./QuantumConstants.js";
import QPPWPreferences from "../../QPPWPreferences.js";
import qppw from "../../QPPWNamespace.js";
import { NumericalMethod } from "./NumericalMethod.js";
import {
  BasePotential,
  AnalyticalPotential,
  InfiniteSquareWellPotential,
  FiniteSquareWellPotential,
  HarmonicOscillatorPotential,
  MorsePotential,
  PoschlTellerPotential,
  RosenMorsePotential,
  EckartPotential,
  AsymmetricTrianglePotential,
  Coulomb1DPotential,
  Coulomb3DPotential,
  TriangularPotential,
} from "./potentials/index.js";

// Re-export potential classes for external use
export {
  BasePotential,
  AnalyticalPotential,
  NumericalPotential,
} from "./potentials/index.js";

// Re-export NumericalMethod for backward compatibility
export { NumericalMethod };

/**
 * Configuration for potential well parameters (for analytical solutions)
 */
export interface WellParameters {
  /** Type of potential well */
  type: PotentialType;
  /** Well width for infinite/finite well (meters) */
  wellWidth?: number;
  /** Well depth for finite square well (Joules, positive value) */
  wellDepth?: number;
  /** Spring constant for harmonic oscillator (N/m) */
  springConstant?: number;
  /** Dissociation energy for Morse potential (Joules) */
  dissociationEnergy?: number;
  /** Equilibrium position for Morse potential (meters) */
  equilibriumPosition?: number;
  /** Potential depth for Pöschl-Teller, Rosen-Morse, and Eckart potentials (Joules) */
  potentialDepth?: number;
  /** Barrier height for Rosen-Morse and Eckart potentials (Joules) */
  barrierHeight?: number;
  /** Slope parameter for asymmetric triangle potential (Joules/meter) */
  slope?: number;
  /** Coulomb strength parameter α for Coulomb potentials (J·m) */
  coulombStrength?: number;
  /** Well separation for double square well (meters) */
  wellSeparation?: number;
  /** Energy offset for triangular potential (Joules) */
  energyOffset?: number;
  /** Number of wells for multi-square well and multi-Coulomb 1D (1-10) */
  numberOfWells?: number;
}

/**
 * Main class for solving the 1D time-independent Schrödinger equation.
 */
export class Schrodinger1DSolver {
  private numericalMethod: NumericalMethod;
  private analyticalSolution: AnalyticalSolution | null = null;

  /**
   * Create a new solver instance.
   * @param method - Numerical method to use (default: DVR)
   */
  constructor(method: NumericalMethod = NumericalMethod.DVR) {
    this.numericalMethod = method;
  }

  /**
   * Set the numerical method to use for calculations.
   * @param method - Numerov or DVR method
   */
  public setNumericalMethod(method: NumericalMethod): void {
    this.numericalMethod = method;
  }

  /**
   * Get the current numerical method setting.
   */
  public getNumericalMethod(): NumericalMethod {
    return this.numericalMethod;
  }

  /**
   * Get the current analytical solution instance (if one has been created).
   * This provides access to additional methods like calculateTurningPoints,
   * calculateWavefunctionZeros, etc.
   */
  public getAnalyticalSolution(): AnalyticalSolution | null {
    return this.analyticalSolution;
  }

  /**
   * Solve the Schrödinger equation using analytical solution if available,
   * otherwise use numerical method.
   *
   * @param wellParams - Parameters defining the potential well (for analytical solutions)
   * @param mass - Particle mass in kg
   * @param numStates - Number of bound states to calculate
   * @param gridConfig - Grid configuration for spatial discretization
   * @returns Bound state results with energies and wavefunctions
   */
  public solveAnalyticalIfPossible(
    wellParams: WellParameters,
    mass: number,
    numStates: number,
    gridConfig: GridConfig,
  ): BoundStateResult {
    // Try to create analytical solution instance and solve
    this.analyticalSolution = this.createAnalyticalSolution(wellParams, mass);

    if (this.analyticalSolution) {
      return this.analyticalSolution.solve(numStates, gridConfig);
    }

    // Handle special cases that don't have classes yet (multi-well potentials)
    const specialCaseResult = this.handleSpecialCases(
      wellParams,
      mass,
      numStates,
      gridConfig,
    );
    if (specialCaseResult) {
      return specialCaseResult;
    }

    // If no analytical solution available, throw error
    throw new Error(
      `No analytical solution available for potential type: ${wellParams.type}. Use solve() method with a custom potential function instead.`,
    );
  }

  /**
   * Create an analytical solution instance based on well parameters.
   * Returns null if the potential type doesn't have an analytical solution
   * or if required parameters are missing.
   *
   * @param wellParams - Parameters defining the potential well
   * @param mass - Particle mass in kg
   * @returns AnalyticalSolution instance or null
   */
  private createAnalyticalSolution(
    wellParams: WellParameters,
    mass: number,
  ): AnalyticalSolution | null {
    switch (wellParams.type) {
      case PotentialType.INFINITE_WELL:
        if (wellParams.wellWidth !== undefined) {
          return new InfiniteSquareWellSolution(wellParams.wellWidth, mass);
        }
        break;

      case PotentialType.FINITE_WELL:
        if (
          wellParams.wellWidth !== undefined &&
          wellParams.wellDepth !== undefined
        ) {
          return new FiniteSquareWellSolution(
            wellParams.wellWidth,
            wellParams.wellDepth,
            mass,
          );
        }
        break;

      case PotentialType.HARMONIC_OSCILLATOR:
        if (wellParams.springConstant !== undefined) {
          return new HarmonicOscillatorSolution(
            wellParams.springConstant,
            mass,
          );
        }
        break;

      case PotentialType.MORSE:
        if (
          wellParams.dissociationEnergy !== undefined &&
          wellParams.wellWidth !== undefined &&
          wellParams.equilibriumPosition !== undefined
        ) {
          return new MorsePotentialSolution(
            wellParams.dissociationEnergy,
            wellParams.wellWidth,
            wellParams.equilibriumPosition,
            mass,
          );
        }
        break;

      case PotentialType.POSCHL_TELLER:
        if (
          wellParams.potentialDepth !== undefined &&
          wellParams.wellWidth !== undefined
        ) {
          return new PoschlTellerPotentialSolution(
            wellParams.potentialDepth,
            wellParams.wellWidth,
            mass,
          );
        }
        break;

      case PotentialType.ROSEN_MORSE:
        if (
          wellParams.potentialDepth !== undefined &&
          wellParams.barrierHeight !== undefined &&
          wellParams.wellWidth !== undefined
        ) {
          return new RosenMorsePotentialSolution(
            wellParams.potentialDepth,
            wellParams.barrierHeight,
            wellParams.wellWidth,
            mass,
          );
        }
        break;

      case PotentialType.ECKART:
        if (
          wellParams.potentialDepth !== undefined &&
          wellParams.barrierHeight !== undefined &&
          wellParams.wellWidth !== undefined
        ) {
          return new EckartPotentialSolution(
            wellParams.potentialDepth,
            wellParams.barrierHeight,
            wellParams.wellWidth,
            mass,
          );
        }
        break;

      case PotentialType.ASYMMETRIC_TRIANGLE:
        if (
          wellParams.slope !== undefined &&
          wellParams.wellWidth !== undefined
        ) {
          return new AsymmetricTrianglePotentialSolution(
            wellParams.slope,
            wellParams.wellWidth,
            mass,
          );
        }
        break;

      case PotentialType.COULOMB_1D:
        if (wellParams.coulombStrength !== undefined) {
          return new Coulomb1DPotentialSolution(
            wellParams.coulombStrength,
            mass,
          );
        }
        break;

      case PotentialType.COULOMB_3D:
        if (wellParams.coulombStrength !== undefined) {
          return new Coulomb3DPotentialSolution(
            wellParams.coulombStrength,
            mass,
          );
        }
        break;

      case PotentialType.TRIANGULAR:
        if (
          wellParams.wellWidth !== undefined &&
          wellParams.wellDepth !== undefined &&
          wellParams.energyOffset !== undefined
        ) {
          return new TriangularPotentialSolution(
            wellParams.wellDepth,
            wellParams.wellWidth,
            wellParams.energyOffset,
            mass,
          );
        }
        break;

      case PotentialType.DOUBLE_SQUARE_WELL:
        // Double well doesn't have a class yet - still uses function
        // Keep existing behavior for backward compatibility
        if (
          wellParams.wellWidth !== undefined &&
          wellParams.wellDepth !== undefined &&
          wellParams.wellSeparation !== undefined
        ) {
          // We'll handle this in the calling method
          return null;
        }
        break;

      case PotentialType.MULTI_SQUARE_WELL:
        // Multi-square well doesn't have a class yet - still uses function
        return null;

      case PotentialType.MULTI_COULOMB_1D:
        // Multi-Coulomb 1D doesn't have a class yet - still uses function
        return null;

      default:
        return null;
    }

    return null;
  }

  /**
   * Handle special cases that don't have analytical solution classes yet.
   * These are typically multi-well potentials that still use standalone functions.
   *
   * @returns BoundStateResult if handled, null otherwise
   */
  private handleSpecialCases(
    wellParams: WellParameters,
    mass: number,
    numStates: number,
    gridConfig: GridConfig,
  ): BoundStateResult | null {
    switch (wellParams.type) {
      case PotentialType.DOUBLE_SQUARE_WELL:
        if (
          wellParams.wellWidth !== undefined &&
          wellParams.wellDepth !== undefined &&
          wellParams.wellSeparation !== undefined
        ) {
          // Always use analytical solution for double square well
          // Solves the transcendental equations from boundary conditions exactly
          return solveDoubleSquareWellAnalytical(
            wellParams.wellWidth,
            wellParams.wellDepth,
            wellParams.wellSeparation,
            mass,
            numStates,
            gridConfig,
          );
        }
        break;

      case PotentialType.MULTI_SQUARE_WELL:
        if (
          wellParams.numberOfWells !== undefined &&
          wellParams.wellWidth !== undefined &&
          wellParams.wellDepth !== undefined &&
          wellParams.wellSeparation !== undefined
        ) {
          return solveMultiSquareWell(
            wellParams.numberOfWells,
            wellParams.wellWidth,
            wellParams.wellDepth,
            wellParams.wellSeparation,
            mass,
            numStates,
            gridConfig,
            this, // Pass solver instance for numerical methods
          );
        }
        break;

      case PotentialType.MULTI_COULOMB_1D:
        if (
          wellParams.numberOfWells !== undefined &&
          wellParams.wellSeparation !== undefined &&
          wellParams.coulombStrength !== undefined
        ) {
          return solveMultiCoulomb1D(
            wellParams.numberOfWells,
            wellParams.wellSeparation,
            wellParams.coulombStrength,
            mass,
            numStates,
            gridConfig,
            this, // Pass solver instance for numerical methods
          );
        }
        break;
    }

    return null;
  }

  /**
   * Solve the Schrödinger equation numerically for an arbitrary potential.
   *
   * Uses a two-step approach:
   * 1. Find energies with selected solver on coarse grid (from preferences)
   * 2. Compute wavefunctions on finer grid (1000 points) using Numerov method
   *
   * @param potential - Function V(x) returning potential energy in Joules
   * @param mass - Particle mass in kg
   * @param numStates - Number of bound states to calculate
   * @param gridConfig - Grid configuration for spatial discretization
   * @param energyRange - Optional energy range for Numerov method [min, max] in Joules
   * @returns Bound state results with energies and wavefunctions
   */
  public solveNumerical(
    potential: PotentialFunction,
    mass: number,
    numStates: number,
    gridConfig: GridConfig,
    energyRange?: [number, number],
  ): BoundStateResult {
    // Configuration constants
    const ENERGIES_ONLY = false;
    const FINE_GRID_POINTS = 1000;

    const { xMin, xMax, numPoints } = gridConfig;

    // Step 1: Find energies (and optionally wavefunctions) using selected solver
    // Use the requested grid size from gridConfig, or fall back to preferences
    const coarseNumPoints =
      numPoints || QPPWPreferences.gridPointsProperty.value;
    const coarseGridConfig: GridConfig = {
      xMin,
      xMax,
      numPoints: coarseNumPoints,
    };

    let result: BoundStateResult;

    if (this.numericalMethod === NumericalMethod.NUMEROV) {
      // Numerov method requires energy range for shooting method
      if (!energyRange) {
        // Estimate energy range based on potential at grid points
        const dx = (xMax - xMin) / (coarseNumPoints - 1);
        let Vmin = Infinity;
        let Vmax = -Infinity;

        for (let i = 0; i < coarseNumPoints; i++) {
          const x = xMin + i * dx;
          const V = potential(x);
          if (V < Vmin) Vmin = V;
          if (V < 1e100 && V > Vmax) Vmax = V; // Ignore infinite barriers
        }

        // Search for bound states between Vmin and Vmax
        energyRange = [Vmin, Vmax];
      }

      // Numerov always returns full result with wavefunctions
      result = solveNumerov(
        potential,
        mass,
        numStates,
        coarseGridConfig,
        energyRange[0],
        energyRange[1],
      );
    } else if (this.numericalMethod === NumericalMethod.MATRIX_NUMEROV) {
      // Matrix Numerov method
      result = solveMatrixNumerov(
        potential,
        mass,
        numStates,
        coarseGridConfig,
        ENERGIES_ONLY,
      );
    } else if (this.numericalMethod === NumericalMethod.DVR) {
      // DVR method
      result = solveDVR(
        potential,
        mass,
        numStates,
        coarseGridConfig,
        ENERGIES_ONLY,
      );
    } else if (this.numericalMethod === NumericalMethod.FGH) {
      // Fourier Grid Hamiltonian method
      result = solveFGH(
        potential,
        mass,
        numStates,
        coarseGridConfig,
        ENERGIES_ONLY,
      );
    } else if (this.numericalMethod === NumericalMethod.QUANTUM_BOUND) {
      // Advanced shooting method with adaptive bracketing
      result = solveQuantumBound(potential, mass, numStates, coarseGridConfig);
    } else {
      // Spectral (Chebyshev) method
      result = solveSpectral(
        potential,
        mass,
        numStates,
        coarseGridConfig,
        ENERGIES_ONLY,
      );
    }

    // Step 2: If ENERGIES_ONLY was true, compute wavefunctions on finer grid using Numerov
    // Otherwise, return the result directly from the solver
    if (ENERGIES_ONLY) {
      const fineGridConfig: GridConfig = {
        xMin,
        xMax,
        numPoints: FINE_GRID_POINTS,
      };

      const wavefunctionResult = computeWavefunctionsNumerov(
        result.energies,
        potential,
        mass,
        fineGridConfig,
      );

      return {
        energies: result.energies,
        wavefunctions: wavefunctionResult.wavefunctions,
        xGrid: wavefunctionResult.xGrid,
        method: result.method,
      };
    } else {
      // Use wavefunctions directly from the solver
      return result;
    }
  }

  /**
   * Create a potential class instance from well parameters.
   * This is the new class-based approach that separates analytical and numerical potentials.
   *
   * @param wellParams - Parameters defining the potential well
   * @param mass - Particle mass in kg
   * @returns BasePotential instance (AnalyticalPotential or NumericalPotential)
   */
  public createPotential(
    wellParams: WellParameters,
    mass: number,
  ): BasePotential | null {
    switch (wellParams.type) {
      case PotentialType.INFINITE_WELL:
        if (wellParams.wellWidth !== undefined) {
          return new InfiniteSquareWellPotential(wellParams.wellWidth, mass);
        }
        break;

      case PotentialType.FINITE_WELL:
        if (
          wellParams.wellWidth !== undefined &&
          wellParams.wellDepth !== undefined
        ) {
          return new FiniteSquareWellPotential(
            wellParams.wellWidth,
            wellParams.wellDepth,
            mass,
          );
        }
        break;

      case PotentialType.HARMONIC_OSCILLATOR:
        if (wellParams.springConstant !== undefined) {
          return new HarmonicOscillatorPotential(
            wellParams.springConstant,
            mass,
          );
        }
        break;

      case PotentialType.MORSE:
        if (
          wellParams.dissociationEnergy !== undefined &&
          wellParams.wellWidth !== undefined &&
          wellParams.equilibriumPosition !== undefined
        ) {
          return new MorsePotential(
            wellParams.dissociationEnergy,
            wellParams.wellWidth,
            wellParams.equilibriumPosition,
            mass,
          );
        }
        break;

      case PotentialType.POSCHL_TELLER:
        if (
          wellParams.potentialDepth !== undefined &&
          wellParams.wellWidth !== undefined
        ) {
          return new PoschlTellerPotential(
            wellParams.potentialDepth,
            wellParams.wellWidth,
            mass,
          );
        }
        break;

      case PotentialType.ROSEN_MORSE:
        if (
          wellParams.potentialDepth !== undefined &&
          wellParams.barrierHeight !== undefined &&
          wellParams.wellWidth !== undefined
        ) {
          return new RosenMorsePotential(
            wellParams.potentialDepth,
            wellParams.barrierHeight,
            wellParams.wellWidth,
            mass,
          );
        }
        break;

      case PotentialType.ECKART:
        if (
          wellParams.potentialDepth !== undefined &&
          wellParams.barrierHeight !== undefined &&
          wellParams.wellWidth !== undefined
        ) {
          return new EckartPotential(
            wellParams.potentialDepth,
            wellParams.barrierHeight,
            wellParams.wellWidth,
            mass,
          );
        }
        break;

      case PotentialType.ASYMMETRIC_TRIANGLE:
        if (
          wellParams.slope !== undefined &&
          wellParams.wellWidth !== undefined
        ) {
          return new AsymmetricTrianglePotential(
            wellParams.slope,
            wellParams.wellWidth,
            mass,
          );
        }
        break;

      case PotentialType.COULOMB_1D:
        if (wellParams.coulombStrength !== undefined) {
          return new Coulomb1DPotential(wellParams.coulombStrength, mass);
        }
        break;

      case PotentialType.COULOMB_3D:
        if (wellParams.coulombStrength !== undefined) {
          return new Coulomb3DPotential(wellParams.coulombStrength, mass);
        }
        break;

      case PotentialType.TRIANGULAR:
        if (
          wellParams.wellWidth !== undefined &&
          wellParams.wellDepth !== undefined &&
          wellParams.energyOffset !== undefined
        ) {
          return new TriangularPotential(
            wellParams.wellDepth,
            wellParams.wellWidth,
            wellParams.energyOffset,
            mass,
          );
        }
        break;

      default:
        // For types not yet converted to classes (multi-well types, custom)
        return null;
    }

    return null;
  }

  /**
   * Solve the Schrödinger equation using a potential class instance.
   * This method automatically uses analytical or numerical methods based on
   * the potential type.
   *
   * @param potential - Potential class instance (AnalyticalPotential or NumericalPotential)
   * @param numStates - Number of bound states to calculate
   * @param gridConfig - Grid configuration for spatial discretization
   * @returns Bound state results with energies and wavefunctions
   */
  public solvePotential(
    potential: BasePotential,
    numStates: number,
    gridConfig: GridConfig,
  ): BoundStateResult {
    if (potential.hasAnalyticalSolution()) {
      // Use analytical solution
      const analyticalPotential = potential as AnalyticalPotential;
      return analyticalPotential.solve(numStates, gridConfig);
    } else {
      // Use numerical solution
      const potentialFunction = potential.createPotential();
      return this.solveNumerical(
        potentialFunction,
        potential.getMass(),
        numStates,
        gridConfig,
      );
    }
  }

  /**
   * Create a potential function for an infinite square well.
   * Centered at x=0, extending from -wellWidth/2 to +wellWidth/2.
   * @param wellWidth - Width of the well in meters
   * @param wellDepth - Depth of the well in Joules (0 inside, depth outside)
   * @returns Potential function V(x)
   */
  public static createInfiniteWellPotential(
    wellWidth: number,
    wellDepth = 1e100,
  ): PotentialFunction {
    const halfWidth = wellWidth / 2;
    return (x: number) => {
      if (x >= -halfWidth && x <= halfWidth) {
        return 0;
      } else {
        return wellDepth; // Very large value to approximate infinity
      }
    };
  }

  /**
   * Create a potential function for a finite square well.
   * @param wellWidth - Width of the well in meters
   * @param wellDepth - Depth of the well in Joules (V=0 outside, V=-depth inside)
   * @param center - Center position of well in meters (default 0)
   * @returns Potential function V(x)
   */
  public static createFiniteWellPotential(
    wellWidth: number,
    wellDepth: number,
    center = 0,
  ): PotentialFunction {
    const halfWidth = wellWidth / 2;
    return (x: number) => {
      const xShifted = x - center;
      if (Math.abs(xShifted) <= halfWidth) {
        return -wellDepth;
      } else {
        return 0;
      }
    };
  }

  /**
   * Convert energy from eV to Joules.
   */
  public static eVToJoules(eV: number): number {
    return eV * QuantumConstants.EV_TO_JOULES;
  }

  /**
   * Convert energy from Joules to eV.
   */
  public static joulesToEV(joules: number): number {
    return joules * QuantumConstants.JOULES_TO_EV;
  }
}

qppw.register("Schrodinger1DSolver", Schrodinger1DSolver);

export default Schrodinger1DSolver;
