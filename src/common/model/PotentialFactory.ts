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
} from "./analytical-solutions/index.js";

import {
  BasePotential,
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

import { PotentialType, WellParameters } from "./PotentialFunction.js";

/**
 * Configuration for creating analytical solutions and potential classes.
 */
interface PotentialConfig {
  /** Constructor for analytical solution class */
  solutionClass: new (...args: any[]) => AnalyticalSolution;
  /** Constructor for potential class */
  potentialClass: new (...args: any[]) => BasePotential;
  /** Extract constructor arguments from well parameters */
  extractArgs: (params: WellParameters, mass: number) => any[];
  /** Required parameter names */
  requiredParams: (keyof WellParameters)[];
}

/**
 * Registry of all potential type configurations.
 * Single source of truth for creating analytical solutions and potentials.
 */
const POTENTIAL_CONFIGS: Partial<Record<PotentialType, PotentialConfig>> = {
  [PotentialType.INFINITE_WELL]: {
    solutionClass: InfiniteSquareWellSolution,
    potentialClass: InfiniteSquareWellPotential,
    extractArgs: (p, m) => [p.wellWidth!, m],
    requiredParams: ["wellWidth"],
  },

  [PotentialType.FINITE_WELL]: {
    solutionClass: FiniteSquareWellSolution,
    potentialClass: FiniteSquareWellPotential,
    extractArgs: (p, m) => [p.wellWidth!, p.wellDepth!, m],
    requiredParams: ["wellWidth", "wellDepth"],
  },

  [PotentialType.HARMONIC_OSCILLATOR]: {
    solutionClass: HarmonicOscillatorSolution,
    potentialClass: HarmonicOscillatorPotential,
    extractArgs: (p, m) => [p.springConstant!, m],
    requiredParams: ["springConstant"],
  },

  [PotentialType.MORSE]: {
    solutionClass: MorsePotentialSolution,
    potentialClass: MorsePotential,
    extractArgs: (p, m) => [
      p.dissociationEnergy!,
      p.wellWidth!,
      p.equilibriumPosition!,
      m,
    ],
    requiredParams: ["dissociationEnergy", "wellWidth", "equilibriumPosition"],
  },

  [PotentialType.POSCHL_TELLER]: {
    solutionClass: PoschlTellerPotentialSolution,
    potentialClass: PoschlTellerPotential,
    extractArgs: (p, m) => [p.potentialDepth!, p.wellWidth!, m],
    requiredParams: ["potentialDepth", "wellWidth"],
  },

  [PotentialType.ROSEN_MORSE]: {
    solutionClass: RosenMorsePotentialSolution,
    potentialClass: RosenMorsePotential,
    extractArgs: (p, m) => [
      p.potentialDepth!,
      p.barrierHeight!,
      p.wellWidth!,
      m,
    ],
    requiredParams: ["potentialDepth", "barrierHeight", "wellWidth"],
  },

  [PotentialType.ECKART]: {
    solutionClass: EckartPotentialSolution,
    potentialClass: EckartPotential,
    extractArgs: (p, m) => [
      p.potentialDepth!,
      p.barrierHeight!,
      p.wellWidth!,
      m,
    ],
    requiredParams: ["potentialDepth", "barrierHeight", "wellWidth"],
  },

  [PotentialType.ASYMMETRIC_TRIANGLE]: {
    solutionClass: AsymmetricTrianglePotentialSolution,
    potentialClass: AsymmetricTrianglePotential,
    extractArgs: (p, m) => [p.slope!, p.wellWidth!, m],
    requiredParams: ["slope", "wellWidth"],
  },

  [PotentialType.COULOMB_1D]: {
    solutionClass: Coulomb1DPotentialSolution,
    potentialClass: Coulomb1DPotential,
    extractArgs: (p, m) => [p.coulombStrength!, m],
    requiredParams: ["coulombStrength"],
  },

  [PotentialType.COULOMB_3D]: {
    solutionClass: Coulomb3DPotentialSolution,
    potentialClass: Coulomb3DPotential,
    extractArgs: (p, m) => [p.coulombStrength!, m],
    requiredParams: ["coulombStrength"],
  },

  [PotentialType.TRIANGULAR]: {
    solutionClass: TriangularPotentialSolution,
    potentialClass: TriangularPotential,
    extractArgs: (p, m) => [p.wellDepth!, p.wellWidth!, p.energyOffset!, m],
    requiredParams: ["wellDepth", "wellWidth", "energyOffset"],
  },
};

/**
 * Factory for creating analytical solutions and potential instances.
 * Eliminates code duplication in Schrodinger1DSolver.
 */
export class PotentialFactory {
  /**
   * Check if all required parameters are present.
   */
  private static hasRequiredParams(
    params: WellParameters,
    required: (keyof WellParameters)[],
  ): boolean {
    return required.every((key) => params[key] !== undefined);
  }

  /**
   * Create an analytical solution instance.
   * Returns null if the potential type doesn't support analytical solutions
   * or required parameters are missing.
   */
  static createAnalyticalSolution(
    wellParams: WellParameters,
    mass: number,
  ): AnalyticalSolution | null {
    const config = POTENTIAL_CONFIGS[wellParams.type];

    if (!config) {
      return null; // No analytical solution for this type
    }

    if (!this.hasRequiredParams(wellParams, config.requiredParams)) {
      return null; // Missing required parameters
    }

    const args = config.extractArgs(wellParams, mass);
    return new config.solutionClass(...args);
  }

  /**
   * Create a potential class instance.
   * Returns null if the potential type isn't supported
   * or required parameters are missing.
   */
  static createPotential(
    wellParams: WellParameters,
    mass: number,
  ): BasePotential | null {
    const config = POTENTIAL_CONFIGS[wellParams.type];

    if (!config) {
      return null; // No potential class for this type
    }

    if (!this.hasRequiredParams(wellParams, config.requiredParams)) {
      return null; // Missing required parameters
    }

    const args = config.extractArgs(wellParams, mass);
    return new config.potentialClass(...args);
  }
}
