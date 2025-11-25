/**
 * Analytical solutions for well-known quantum potentials.
 * These provide exact solutions without numerical approximation.
 *
 * This module exports all analytical solution functions and classes for various quantum potentials.
 */

import qppw from "../../../QPPWNamespace.js";

// Export the abstract base class
export { AnalyticalSolution } from "./AnalyticalSolution.js";

// Export all analytical solution functions and classes
export {
  solveInfiniteWell,
  createInfiniteWellPotential,
  calculateInfiniteWellClassicalProbability,
  calculateInfiniteWellWavefunctionZeros,
  calculateInfiniteWellTurningPoints,
  calculateInfiniteWellWavefunctionSecondDerivative,
  InfiniteSquareWellSolution,
} from "./infinite-square-well.js";
export {
  solveFiniteSquareWell,
  createFiniteWellPotential,
  calculateFiniteWellClassicalProbability,
  calculateFiniteWellWavefunctionZeros,
  calculateFiniteWellTurningPoints,
  calculateFiniteWellWavefunctionSecondDerivative,
  FiniteSquareWellSolution,
} from "./finite-square-well.js";
export {
  solveHarmonicOscillator,
  createHarmonicOscillatorPotential,
  calculateHarmonicOscillatorClassicalProbability,
  calculateHarmonicOscillatorWavefunctionZeros,
  calculateHarmonicOscillatorTurningPoints,
  calculateHarmonicOscillatorWavefunctionSecondDerivative,
  HarmonicOscillatorSolution,
} from "./harmonic-oscillator.js";
export {
  solveMorsePotential,
  createMorsePotential,
  calculateMorsePotentialClassicalProbability,
  calculateMorsePotentialWavefunctionZeros,
  calculateMorsePotentialTurningPoints,
  calculateMorsePotentialWavefunctionSecondDerivative,
  MorsePotentialSolution,
} from "./morse-potential.js";
export {
  solvePoschlTellerPotential,
  createPoschlTellerPotential,
  calculatePoschlTellerClassicalProbability,
  calculatePoschlTellerWavefunctionZeros,
  calculatePoschlTellerTurningPoints,
  calculatePoschlTellerWavefunctionSecondDerivative,
  PoschlTellerPotentialSolution,
} from "./poschl-teller-potential.js";
export {
  solveRosenMorsePotential,
  createRosenMorsePotential,
  calculateRosenMorsePotentialClassicalProbability,
  calculateRosenMorsePotentialWavefunctionZeros,
  calculateRosenMorsePotentialTurningPoints,
  calculateRosenMorsePotentialWavefunctionSecondDerivative,
  RosenMorsePotentialSolution,
} from "./rosen-morse-potential.js";
export {
  solveEckartPotential,
  createEckartPotential,
  calculateEckartPotentialClassicalProbability,
  calculateEckartPotentialWavefunctionZeros,
  calculateEckartPotentialTurningPoints,
  calculateEckartPotentialWavefunctionSecondDerivative,
  EckartPotentialSolution,
} from "./eckart-potential.js";
export {
  solveAsymmetricTrianglePotential,
  createAsymmetricTrianglePotential,
  calculateAsymmetricTriangleClassicalProbability,
  calculateAsymmetricTriangleWavefunctionZeros,
  calculateAsymmetricTriangleTurningPoints,
  calculateAsymmetricTriangleWavefunctionSecondDerivative,
  AsymmetricTrianglePotentialSolution,
} from "./asymmetric-triangle-potential.js";
export {
  solveCoulomb1DPotential,
  createCoulomb1DPotential,
  calculateCoulomb1DClassicalProbability,
  calculateCoulomb1DWavefunctionZeros,
  calculateCoulomb1DTurningPoints,
  calculateCoulomb1DWavefunctionSecondDerivative,
  Coulomb1DPotentialSolution,
} from "./coulomb-1d-potential.js";
export { solveCoulomb1DNumerical } from "./coulomb-1d-numerical-wrapper.js";
export {
  solveCoulomb3DPotential,
  createCoulomb3DPotential,
  calculateCoulomb3DClassicalProbability,
  calculateCoulomb3DWavefunctionZeros,
  calculateCoulomb3DTurningPoints,
  calculateCoulomb3DWavefunctionSecondDerivative,
  Coulomb3DPotentialSolution,
} from "./coulomb-3d-potential.js";
export {
  solveTriangularPotential,
  createTriangularPotential,
  calculateTriangularPotentialClassicalProbability,
  calculateTriangularPotentialWavefunctionZeros,
  calculateTriangularPotentialTurningPoints,
  calculateTriangularPotentialWavefunctionSecondDerivative,
  TriangularPotentialSolution,
} from "./triangular-potential.js";
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
