/**
 * Analytical solutions for well-known quantum potentials.
 * These provide exact solutions without numerical approximation.
 *
 * This module exports all analytical solution functions for various quantum potentials.
 */

import qppw from "../../../QPPWNamespace.js";

// Export all analytical solution functions
export { solveInfiniteWell } from "./infinite-square-well.js";
export { solveFiniteSquareWell } from "./finite-square-well.js";
export { solveHarmonicOscillator } from "./harmonic-oscillator.js";
export { solveMorsePotential } from "./morse-potential.js";
export { solvePoschlTellerPotential } from "./poschl-teller-potential.js";
export { solveRosenMorsePotential } from "./rosen-morse-potential.js";
export { solveEckartPotential } from "./eckart-potential.js";
export { solveAsymmetricTrianglePotential } from "./asymmetric-triangle-potential.js";
export { solveCoulomb1DPotential } from "./coulomb-1d-potential.js";
export { solveCoulomb3DPotential } from "./coulomb-3d-potential.js";

// Export utility functions that may be useful externally
export * from "./math-utilities.js";

// Import all functions for registration
import { solveInfiniteWell } from "./infinite-square-well.js";
import { solveFiniteSquareWell } from "./finite-square-well.js";
import { solveHarmonicOscillator } from "./harmonic-oscillator.js";
import { solveMorsePotential } from "./morse-potential.js";
import { solvePoschlTellerPotential } from "./poschl-teller-potential.js";
import { solveRosenMorsePotential } from "./rosen-morse-potential.js";
import { solveEckartPotential } from "./eckart-potential.js";
import { solveAsymmetricTrianglePotential } from "./asymmetric-triangle-potential.js";
import { solveCoulomb1DPotential } from "./coulomb-1d-potential.js";
import { solveCoulomb3DPotential } from "./coulomb-3d-potential.js";

// Register all solutions with the QPPW namespace
qppw.register("AnalyticalSolutions", {
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
});
