/**
 * Simple test to verify the new potential class hierarchy.
 * This demonstrates the separation of analytical and numerical potentials.
 */

import { Schrodinger1DSolver, WellParameters } from "./src/common/model/Schrodinger1DSolver.js";
import { PotentialType } from "./src/common/model/PotentialFunction.js";
import { NumericalPotential } from "./src/common/model/potentials/index.js";
import QuantumConstants from "./src/common/model/QuantumConstants.js";

// Test configuration
const mass = QuantumConstants.ELECTRON_MASS;
const numStates = 5;
const gridConfig = {
  xMin: -1e-9,
  xMax: 1e-9,
  numPoints: 500,
};

console.log("Testing Potential Class Hierarchy");
console.log("=".repeat(50));

// Test 1: Create analytical potentials
console.log("\n1. Testing Analytical Potentials:");
console.log("-".repeat(50));

const analyticalTypes = [
  PotentialType.INFINITE_WELL,
  PotentialType.FINITE_WELL,
  PotentialType.HARMONIC_OSCILLATOR,
  PotentialType.MORSE,
];

const solver = new Schrodinger1DSolver();

for (const type of analyticalTypes) {
  let wellParams: WellParameters;

  switch (type) {
    case PotentialType.INFINITE_WELL:
      wellParams = { type, wellWidth: 1e-9 };
      break;
    case PotentialType.FINITE_WELL:
      wellParams = { type, wellWidth: 1e-9, wellDepth: 1e-18 };
      break;
    case PotentialType.HARMONIC_OSCILLATOR:
      wellParams = { type, springConstant: 1.0 };
      break;
    case PotentialType.MORSE:
      wellParams = {
        type,
        dissociationEnergy: 1e-18,
        wellWidth: 1e-10,
        equilibriumPosition: 0,
      };
      break;
    default:
      continue;
  }

  const potential = solver.createPotential(wellParams, mass);

  if (potential) {
    console.log(`\n${type}:`);
    console.log(`  - Has analytical solution: ${potential.hasAnalyticalSolution()}`);
    console.log(`  - Type: ${potential.getType()}`);
    console.log(`  - Class: ${potential.constructor.name}`);

    try {
      const result = solver.solvePotential(potential, numStates, gridConfig);
      console.log(`  - Method used: ${result.method}`);
      console.log(`  - Number of energies found: ${result.energies.length}`);
      console.log(`  - First 3 energies (eV): ${result.energies.slice(0, 3).map(e =>
        (e * QuantumConstants.JOULES_TO_EV).toExponential(3)
      ).join(", ")}`);
    } catch (error) {
      console.log(`  - Error: ${error}`);
    }
  }
}

// Test 2: Create a numerical potential
console.log("\n\n2. Testing Numerical Potential:");
console.log("-".repeat(50));

// Create a custom potential (parabolic well)
const customPotential = (x: number) => {
  const k = 1.0; // spring constant
  return 0.5 * k * x * x;
};

const numericalPot = new NumericalPotential(customPotential, mass);
console.log(`\nCustom Parabolic Potential:`);
console.log(`  - Has analytical solution: ${numericalPot.hasAnalyticalSolution()}`);
console.log(`  - Type: ${numericalPot.getType()}`);
console.log(`  - Class: ${numericalPot.constructor.name}`);

try {
  const result = solver.solvePotential(numericalPot, numStates, gridConfig);
  console.log(`  - Method used: ${result.method}`);
  console.log(`  - Number of energies found: ${result.energies.length}`);
  console.log(`  - First 3 energies (eV): ${result.energies.slice(0, 3).map(e =>
    (e * QuantumConstants.JOULES_TO_EV).toExponential(3)
  ).join(", ")}`);
} catch (error) {
  console.log(`  - Error: ${error}`);
}

// Test 3: Verify solver can distinguish between analytical and numerical
console.log("\n\n3. Verifying Analytical vs Numerical Distinction:");
console.log("-".repeat(50));

const infiniteWellParams: WellParameters = {
  type: PotentialType.INFINITE_WELL,
  wellWidth: 1e-9,
};

const analyticalPot = solver.createPotential(infiniteWellParams, mass);

if (analyticalPot) {
  console.log("\nAnalytical Potential (Infinite Well):");
  console.log(`  - hasAnalyticalSolution(): ${analyticalPot.hasAnalyticalSolution()}`);

  const result1 = solver.solvePotential(analyticalPot, numStates, gridConfig);
  console.log(`  - Solution method: ${result1.method}`);
}

console.log("\nNumerical Potential (Custom):");
console.log(`  - hasAnalyticalSolution(): ${numericalPot.hasAnalyticalSolution()}`);

const result2 = solver.solvePotential(numericalPot, numStates, gridConfig);
console.log(`  - Solution method: ${result2.method}`);

console.log("\n" + "=".repeat(50));
console.log("All tests completed!");
