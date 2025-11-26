/**
 * Exports for quantum potential classes.
 *
 * This module provides the base classes and implementations for both
 * analytical and numerical quantum mechanical potentials.
 */

export { BasePotential } from "./BasePotential.js";
export { AnalyticalPotential } from "./AnalyticalPotential.js";
export { NumericalPotential } from "./NumericalPotential.js";

// Specific analytical potential implementations will be added here
export { InfiniteSquareWellPotential } from "./InfiniteSquareWellPotential.js";
export { FiniteSquareWellPotential } from "./FiniteSquareWellPotential.js";
export { HarmonicOscillatorPotential } from "./HarmonicOscillatorPotential.js";
export { MorsePotential } from "./MorsePotential.js";
export { PoschlTellerPotential } from "./PoschlTellerPotential.js";
export { RosenMorsePotential } from "./RosenMorsePotential.js";
export { EckartPotential } from "./EckartPotential.js";
export { AsymmetricTrianglePotential } from "./AsymmetricTrianglePotential.js";
export { Coulomb1DPotential } from "./Coulomb1DPotential.js";
export { Coulomb3DPotential } from "./Coulomb3DPotential.js";
export { TriangularPotential } from "./TriangularPotential.js";
