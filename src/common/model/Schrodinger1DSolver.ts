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
  solveFiniteSquareWell,
  solveInfiniteWell,
  solveHarmonicOscillator,
  solveMorsePotential,
  solvePoschlTellerPotential,
  solveRosenMorsePotential,
  solveEckartPotential,
  solveAsymmetricTrianglePotential,
  solveCoulomb1DPotential,
  solveCoulomb3DPotential,
  solveTriangularPotential
} from "./analytical-solutions";
import { solveNumerov } from "./NumerovSolver.js";
import { solveMatrixNumerov } from "./MatrixNumerovSolver.js";
import { solveDVR } from "./DVRSolver.js";
import { solveFGH } from "./FGHSolver.js";
import { solveSpectral } from "./SpectralSolver.js";
import { solveDoubleWellNumerov } from "./DoubleWellNumerovSolver.js";
import QuantumConstants from "./QuantumConstants.js";
import qppw from "../../QPPWNamespace.js";

/**
 * Numerical method selection for solving the Schrödinger equation
 */
export enum NumericalMethod {
  NUMEROV = "numerov",
  MATRIX_NUMEROV = "matrix_numerov",
  DVR = "dvr",
  FGH = "fgh",
  SPECTRAL = "spectral",
}

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
}

/**
 * Main class for solving the 1D time-independent Schrödinger equation.
 */
export class Schrodinger1DSolver {
  private numericalMethod: NumericalMethod;

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
    // Try analytical solution first
    switch (wellParams.type) {
      case PotentialType.INFINITE_WELL:
        if (wellParams.wellWidth !== undefined) {
          return solveInfiniteWell(
            wellParams.wellWidth,
            mass,
            numStates,
            gridConfig,
          );
        }
        break;

      case PotentialType.FINITE_WELL:
        if (wellParams.wellWidth !== undefined && wellParams.wellDepth !== undefined) {
          return solveFiniteSquareWell(
            wellParams.wellWidth,
            wellParams.wellDepth,
            mass,
            numStates,
            gridConfig,
          );
        }
        break;

      case PotentialType.HARMONIC_OSCILLATOR:
        if (wellParams.springConstant !== undefined) {
          return solveHarmonicOscillator(
            wellParams.springConstant,
            mass,
            numStates,
            gridConfig,
          );
        }
        break;

      case PotentialType.MORSE:
        if (
          wellParams.dissociationEnergy !== undefined &&
          wellParams.wellWidth !== undefined &&
          wellParams.equilibriumPosition !== undefined
        ) {
          return solveMorsePotential(
            wellParams.dissociationEnergy,
            wellParams.wellWidth,
            wellParams.equilibriumPosition,
            mass,
            numStates,
            gridConfig,
          );
        }
        break;

      case PotentialType.POSCHL_TELLER:
        if (
          wellParams.potentialDepth !== undefined &&
          wellParams.wellWidth !== undefined
        ) {
          return solvePoschlTellerPotential(
            wellParams.potentialDepth,
            wellParams.wellWidth,
            mass,
            numStates,
            gridConfig,
          );
        }
        break;

      case PotentialType.ROSEN_MORSE:
        if (
          wellParams.potentialDepth !== undefined &&
          wellParams.barrierHeight !== undefined &&
          wellParams.wellWidth !== undefined
        ) {
          return solveRosenMorsePotential(
            wellParams.potentialDepth,
            wellParams.barrierHeight,
            wellParams.wellWidth,
            mass,
            numStates,
            gridConfig,
          );
        }
        break;

      case PotentialType.ECKART:
        if (
          wellParams.potentialDepth !== undefined &&
          wellParams.barrierHeight !== undefined &&
          wellParams.wellWidth !== undefined
        ) {
          return solveEckartPotential(
            wellParams.potentialDepth,
            wellParams.barrierHeight,
            wellParams.wellWidth,
            mass,
            numStates,
            gridConfig,
          );
        }
        break;

      case PotentialType.ASYMMETRIC_TRIANGLE:
        if (
          wellParams.slope !== undefined &&
          wellParams.wellWidth !== undefined
        ) {
          return solveAsymmetricTrianglePotential(
            wellParams.slope,
            wellParams.wellWidth,
            mass,
            numStates,
            gridConfig,
          );
        }
        break;

      case PotentialType.COULOMB_1D:
        if (wellParams.coulombStrength !== undefined) {
          return solveCoulomb1DPotential(
            wellParams.coulombStrength,
            mass,
            numStates,
            gridConfig,
          );
        }
        break;

      case PotentialType.COULOMB_3D:
        if (wellParams.coulombStrength !== undefined) {
          return solveCoulomb3DPotential(
            wellParams.coulombStrength,
            mass,
            numStates,
            gridConfig,
          );
        }
        break;

      case PotentialType.TRIANGULAR:
        if (
          wellParams.potentialDepth !== undefined &&
          wellParams.wellWidth !== undefined &&
          wellParams.energyOffset !== undefined
        ) {
          return solveTriangularPotential(
            wellParams.potentialDepth,
            wellParams.wellWidth,
            wellParams.energyOffset,
            mass,
            numStates,
            gridConfig,
          );
        }
        break;

      case PotentialType.DOUBLE_SQUARE_WELL:
        if (
          wellParams.wellWidth !== undefined &&
          wellParams.wellDepth !== undefined &&
          wellParams.wellSeparation !== undefined
        ) {
          // Use specialized double well solver for Numerov method
          // This provides much better energy search by using single well approximation
          if (this.numericalMethod === NumericalMethod.NUMEROV) {
            return solveDoubleWellNumerov(
              wellParams.wellWidth,
              wellParams.wellDepth,
              wellParams.wellSeparation,
              mass,
              numStates,
              gridConfig,
            );
          }

          // For other methods (DVR, FGH, Spectral), use generic potential
          const wellWidth = wellParams.wellWidth;
          const wellDepth = wellParams.wellDepth;
          const separation = wellParams.wellSeparation;

          // Two square wells centered at ±(separation/2 + wellWidth/2)
          // Wells are symmetric about x=0
          const leftWellCenter = -(separation / 2 + wellWidth / 2);
          const rightWellCenter = separation / 2 + wellWidth / 2;

          const potential: PotentialFunction = (x: number) => {
            // Check if x is inside either well
            const inLeftWell = Math.abs(x - leftWellCenter) <= wellWidth / 2;
            const inRightWell = Math.abs(x - rightWellCenter) <= wellWidth / 2;

            if (inLeftWell || inRightWell) {
              return 0; // Ground level is 0 inside the wells
            } else {
              return wellDepth; // Potential height outside the wells
            }
          };

          const result = this.solveNumerical(potential, mass, numStates, gridConfig);
          return result;
        }
        break;

      case PotentialType.CUSTOM:
        // No analytical solution, fall through to numerical
        break;
    }

    // If we reach here, we need to use numerical method
    // But we need a potential function, which should be provided separately
    throw new Error(
      "Analytical solution not available. Use solveNumerical() with a custom potential function.",
    );
  }

  /**
   * Solve the Schrödinger equation numerically for an arbitrary potential.
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
    if (this.numericalMethod === NumericalMethod.NUMEROV) {
      // Numerov method requires energy range for shooting method
      if (!energyRange) {
        // Estimate energy range based on potential at grid points
        const { xMin, xMax, numPoints } = gridConfig;
        const dx = (xMax - xMin) / (numPoints - 1);
        let Vmin = Infinity;
        let Vmax = -Infinity;

        for (let i = 0; i < numPoints; i++) {
          const x = xMin + i * dx;
          const V = potential(x);
          if (V < Vmin) Vmin = V;
          if (V < 1e100 && V > Vmax) Vmax = V; // Ignore infinite barriers
        }

        // Search for bound states between Vmin and Vmax
        energyRange = [Vmin, Vmax];
      }

      return solveNumerov(
        potential,
        mass,
        numStates,
        gridConfig,
        energyRange[0],
        energyRange[1],
      );
    } else if (this.numericalMethod === NumericalMethod.MATRIX_NUMEROV) {
      // Matrix Numerov method
      return solveMatrixNumerov(potential, mass, numStates, gridConfig);
    } else if (this.numericalMethod === NumericalMethod.DVR) {
      // DVR method
      return solveDVR(potential, mass, numStates, gridConfig);
    } else if (this.numericalMethod === NumericalMethod.FGH) {
      // Fourier Grid Hamiltonian method
      return solveFGH(potential, mass, numStates, gridConfig);
    } else {
      // Spectral (Chebyshev) method
      return solveSpectral(potential, mass, numStates, gridConfig);
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
