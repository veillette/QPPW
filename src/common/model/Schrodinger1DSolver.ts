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
  solveCoulomb3DPotential
} from "./AnalyticalSolutions.js";
import { solveNumerov } from "./NumerovSolver.js";
import { solveDVR } from "./DVRSolver.js";
import { solveFGH } from "./FGHSolver.js";
import { solveSpectral } from "./SpectralSolver.js";
import QuantumConstants from "./QuantumConstants.js";
import qppw from "../../QPPWNamespace.js";

/**
 * Numerical method selection for solving the Schrödinger equation
 */
export enum NumericalMethod {
  NUMEROV = "numerov",
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
   * @param wellWidth - Width of the well in meters
   * @param wellDepth - Depth of the well in Joules (0 inside, depth outside)
   * @returns Potential function V(x)
   */
  public static createInfiniteWellPotential(
    wellWidth: number,
    wellDepth = 1e100,
  ): PotentialFunction {
    return (x: number) => {
      if (x >= 0 && x <= wellWidth) {
        return 0;
      } else {
        return wellDepth; // Very large value to approximate infinity
      }
    };
  }

  /**
   * Create a potential function for a harmonic oscillator.
   * @param springConstant - Spring constant k in N/m
   * @param center - Center position of oscillator in meters (default 0)
   * @returns Potential function V(x) = (1/2) * k * (x - center)²
   */
  public static createHarmonicOscillatorPotential(
    springConstant: number,
    center = 0,
  ): PotentialFunction {
    return (x: number) => {
      const displacement = x - center;
      return 0.5 * springConstant * displacement * displacement;
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
   * Create a potential function for a Morse potential.
   * V(x) = D_e * (1 - exp(-a(x - x_e)))^2
   *
   * @param dissociationEnergy - Dissociation energy D_e in Joules
   * @param wellWidth - Width parameter a (inverse meters)
   * @param equilibriumPosition - Equilibrium position x_e in meters
   * @returns Potential function V(x)
   */
  public static createMorsePotential(
    dissociationEnergy: number,
    wellWidth: number,
    equilibriumPosition: number,
  ): PotentialFunction {
    return (x: number) => {
      const exponent = Math.exp(-wellWidth * (x - equilibriumPosition));
      return dissociationEnergy * Math.pow(1 - exponent, 2);
    };
  }

  /**
   * Create a potential function for a Pöschl-Teller potential.
   * V(x) = -V_0 / cosh²(ax)
   *
   * @param potentialDepth - Potential depth V_0 in Joules (positive value)
   * @param wellWidth - Width parameter a (inverse meters)
   * @returns Potential function V(x)
   */
  public static createPoschlTellerPotential(
    potentialDepth: number,
    wellWidth: number,
  ): PotentialFunction {
    return (x: number) => {
      const coshVal = Math.cosh(wellWidth * x);
      return -potentialDepth / (coshVal * coshVal);
    };
  }

  /**
   * Create a potential function for a Rosen-Morse potential.
   * V(x) = -V_0 / cosh²(ax) + V_1 * tanh(ax)
   *
   * @param potentialDepth - Potential depth V_0 in Joules (positive value)
   * @param barrierHeight - Barrier height V_1 in Joules
   * @param wellWidth - Width parameter a (inverse meters)
   * @returns Potential function V(x)
   */
  public static createRosenMorsePotential(
    potentialDepth: number,
    barrierHeight: number,
    wellWidth: number,
  ): PotentialFunction {
    return (x: number) => {
      const coshVal = Math.cosh(wellWidth * x);
      const tanhVal = Math.tanh(wellWidth * x);
      return -potentialDepth / (coshVal * coshVal) + barrierHeight * tanhVal;
    };
  }

  /**
   * Create a potential function for an Eckart potential.
   * V(x) = V_0 / (1 + exp(ax))² - V_1 / (1 + exp(ax))
   *
   * @param potentialDepth - Potential depth V_0 in Joules
   * @param barrierHeight - Barrier height V_1 in Joules
   * @param wellWidth - Width parameter a (inverse meters)
   * @returns Potential function V(x)
   */
  public static createEckartPotential(
    potentialDepth: number,
    barrierHeight: number,
    wellWidth: number,
  ): PotentialFunction {
    return (x: number) => {
      const expVal = Math.exp(wellWidth * x);
      const denom = 1 + expVal;
      return potentialDepth / (denom * denom) - barrierHeight / denom;
    };
  }

  /**
   * Create a potential function for an asymmetric triangle potential.
   * V(x) = 0 for x < 0
   * V(x) = -b(a-x) for 0 < x < a
   * V(x) = 0 for x > a
   *
   * @param slope - Slope parameter b in Joules/meter (positive value)
   * @param wellWidth - Width parameter a in meters
   * @returns Potential function V(x)
   */
  public static createAsymmetricTrianglePotential(
    slope: number,
    wellWidth: number,
  ): PotentialFunction {
    return (x: number) => {
      if (x < 0) {
        return 0;
      } else if (x <= wellWidth) {
        return -slope * (wellWidth - x);
      } else {
        return 0;
      }
    };
  }

  /**
   * Create a potential function for a 1D Coulomb potential.
   * V(x) = -α/|x|
   *
   * This potential has a singularity at x=0 and describes a 1D hydrogen-like atom.
   *
   * @param coulombStrength - Coulomb strength parameter α in J·m
   * @returns Potential function V(x)
   */
  public static createCoulomb1DPotential(
    coulombStrength: number,
  ): PotentialFunction {
    return (x: number) => {
      if (x === 0) {
        return -Infinity; // Singularity at x=0
      }
      return -coulombStrength / Math.abs(x);
    };
  }

  /**
   * Create a potential function for a 3D Coulomb potential (radial).
   * V(r) = -α/r
   *
   * This is the radial potential for the hydrogen atom.
   * Note: r should always be positive for the radial coordinate.
   *
   * @param coulombStrength - Coulomb strength parameter α in J·m
   * @returns Potential function V(r)
   */
  public static createCoulomb3DPotential(
    coulombStrength: number,
  ): PotentialFunction {
    return (r: number) => {
      if (r <= 0) {
        return -Infinity; // Singularity at r=0
      }
      return -coulombStrength / r;
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
