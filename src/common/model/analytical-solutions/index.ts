/**
 * Analytical solutions for well-known quantum potentials.
 * These provide exact solutions without numerical approximation.
 *
 * This module exports all analytical solution functions for various quantum potentials.
 */

import qppw from "../../../QPPWNamespace.js";

// Export all analytical solution functions
export {
  solveInfiniteWell,
  createInfiniteWellPotential,
  calculateInfiniteWellClassicalProbability,
  calculateInfiniteWellWavefunctionZeros,
  calculateInfiniteWellTurningPoints,
  calculateInfiniteWellWavefunctionDerivatives,
} from "./infinite-square-well.js";
export {
  solveFiniteSquareWell,
  createFiniteWellPotential,
  calculateFiniteWellClassicalProbability,
  calculateFiniteWellWavefunctionZeros,
  calculateFiniteWellTurningPoints,
  calculateFiniteWellWavefunctionDerivatives,
} from "./finite-square-well.js";
export {
  solveHarmonicOscillator,
  createHarmonicOscillatorPotential,
  calculateHarmonicOscillatorClassicalProbability,
  calculateHarmonicOscillatorWavefunctionZeros,
  calculateHarmonicOscillatorTurningPoints,
  calculateHarmonicOscillatorWavefunctionDerivatives,
} from "./harmonic-oscillator.js";
export { solveMorsePotential } from "./morse-potential.js";
export { solvePoschlTellerPotential } from "./poschl-teller-potential.js";
export { solveRosenMorsePotential } from "./rosen-morse-potential.js";
export { solveEckartPotential } from "./eckart-potential.js";
export { solveAsymmetricTrianglePotential } from "./asymmetric-triangle-potential.js";
export { solveCoulomb1DPotential } from "./coulomb-1d-potential.js";
export { solveCoulomb1DNumerical } from "./coulomb-1d-numerical-wrapper.js";
export { solveCoulomb3DPotential } from "./coulomb-3d-potential.js";
export { solveTriangularPotential } from "./triangular-potential.js";
export { solveDoubleSquareWellAnalytical } from "./double-square-well.js";
export { solveMultiSquareWell } from "./multi-square-well.js";
export { solveMultiCoulomb1D } from "./multi-coulomb-1d.js";

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
import { solveTriangularPotential } from "./triangular-potential.js";
import { solveDoubleSquareWellAnalytical } from "./double-square-well.js";
import { solveMultiSquareWell } from "./multi-square-well.js";
import { solveMultiCoulomb1D } from "./multi-coulomb-1d.js";

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
  solveCoulomb3DPotential,
  solveTriangularPotential,
  solveDoubleSquareWellAnalytical,
  solveMultiSquareWell,
  solveMultiCoulomb1D,
});
