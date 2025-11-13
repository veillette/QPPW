/**
 * Type definitions and interfaces for potential energy functions.
 */

import qppw from "../../QPPWNamespace.js";

/**
 * Type of potential well for analytical solution selection
 */
export enum PotentialType {
  INFINITE_WELL = "infiniteWell",
  HARMONIC_OSCILLATOR = "harmonicOscillator",
  MORSE = "morse",
  POSCHL_TELLER = "poschlTeller",
  ROSEN_MORSE = "rosenMorse",
  ECKART = "eckart",
  ASYMMETRIC_TRIANGLE = "asymmetricTriangle",
  COULOMB_1D = "coulomb1D",
  COULOMB_3D = "coulomb3D",
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
  method: "analytical" | "numerov" | "dvr";
}

qppw.register("PotentialFunction", { PotentialType });

export default PotentialFunction;
