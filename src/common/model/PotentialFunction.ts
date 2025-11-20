/**
 * Type definitions and interfaces for potential energy functions.
 */

import qppw from "../../QPPWNamespace.js";

/**
 * Type of potential well for analytical solution selection
 */
export enum PotentialType {
  INFINITE_WELL = "infiniteWell",
  FINITE_WELL = "finiteWell",
  HARMONIC_OSCILLATOR = "harmonicOscillator",
  MORSE = "morse",
  POSCHL_TELLER = "poschlTeller",
  ROSEN_MORSE = "rosenMorse",
  ECKART = "eckart",
  ASYMMETRIC_TRIANGLE = "asymmetricTriangle",
  COULOMB_1D = "coulomb1D",
  COULOMB_3D = "coulomb3D",
  DOUBLE_SQUARE_WELL = "doubleSquareWell",
  TRIANGULAR = "triangular",
  CUSTOM = "custom",
}

/**
 * A function that returns the potential energy at a given position.
 * @param x - Position in meters
 * @returns Potential energy in Joules
 */
export type PotentialFunction = (x: number) => number;

/**
 * Configuration for the spatial grid used in numerical calculations
 */
export interface GridConfig {
  /** Minimum x value in meters */
  xMin: number;
  /** Maximum x value in meters */
  xMax: number;
  /** Number of grid points */
  numPoints: number;
}

/**
 * Valid solver method identifiers
 */
export type SolverMethod = "analytical" | "numerov" | "matrix_numerov" | "dvr" | "fgh" | "spectral" | "quantum_bound";

/**
 * Result from solving the Schrödinger equation
 */
export interface BoundStateResult {
  /** Energy eigenvalues in Joules */
  energies: number[];
  /** Wavefunctions (each row is one eigenstate) */
  wavefunctions: number[][];
  /** Grid x-positions in meters */
  xGrid: number[];
  /** Whether analytical or numerical solution was used */
  method: SolverMethod;
}

/**
 * Result from solving only for energies (no wavefunctions)
 */
export interface EnergyOnlyResult {
  /** Energy eigenvalues in Joules */
  energies: number[];
  /** Method used to find energies */
  method: SolverMethod;
}

qppw.register("PotentialFunction", { PotentialType });

export default PotentialFunction;
