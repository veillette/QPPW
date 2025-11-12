/**
 * Simple test script for Morse and Pöschl-Teller potentials
 * Run with: node test-morse-poschl-teller.js
 */

import Schrodinger1DSolver from './src/common/model/Schrodinger1DSolver.js';
import { PotentialType } from './src/common/model/PotentialFunction.js';
import QuantumConstants from './src/common/model/QuantumConstants.js';

console.log("Testing Morse and Pöschl-Teller Analytical Solvers\n");
console.log("=".repeat(60));

// Test 1: Morse Potential
try {
  console.log("\n1. Testing Morse Potential");
  console.log("-".repeat(60));

  const solver1 = new Schrodinger1DSolver();
  const mass = QuantumConstants.ELECTRON_MASS;

  const dissociationEnergy = Schrodinger1DSolver.eVToJoules(5); // 5 eV
  const wellWidth = 1e10; // 1/Å
  const equilibriumPosition = 0;

  const gridConfig1 = {
    xMin: -2e-10,
    xMax: 5e-10,
    numPoints: 200,
  };

  const result1 = solver1.solveAnalyticalIfPossible(
    {
      type: PotentialType.MORSE,
      dissociationEnergy: dissociationEnergy,
      wellWidth: wellWidth,
      equilibriumPosition: equilibriumPosition,
    },
    mass,
    10, // numStates
    gridConfig1,
  );

  console.log(`Method: ${result1.method}`);
  console.log(`Number of bound states: ${result1.energies.length}`);
  console.log(`Grid points: ${result1.xGrid.length}`);
  console.log(`Energy levels (first 5):`);
  for (let i = 0; i < Math.min(5, result1.energies.length); i++) {
    const E_eV = Schrodinger1DSolver.joulesToEV(result1.energies[i]);
    console.log(`  E_${i} = ${E_eV.toFixed(6)} eV`);
  }
  console.log("✓ Morse potential test PASSED");

} catch (error) {
  console.error("✗ Morse potential test FAILED:", error.message);
}

// Test 2: Pöschl-Teller Potential
try {
  console.log("\n2. Testing Pöschl-Teller Potential");
  console.log("-".repeat(60));

  const solver2 = new Schrodinger1DSolver();
  const mass = QuantumConstants.ELECTRON_MASS;

  const potentialDepth = Schrodinger1DSolver.eVToJoules(10); // 10 eV
  const wellWidth = 5e9; // 5/nm

  const gridConfig2 = {
    xMin: -1e-9,
    xMax: 1e-9,
    numPoints: 200,
  };

  const result2 = solver2.solveAnalyticalIfPossible(
    {
      type: PotentialType.POSCHL_TELLER,
      potentialDepth: potentialDepth,
      wellWidth: wellWidth,
    },
    mass,
    5, // numStates
    gridConfig2,
  );

  console.log(`Method: ${result2.method}`);
  console.log(`Number of bound states: ${result2.energies.length}`);
  console.log(`Grid points: ${result2.xGrid.length}`);
  console.log(`Energy levels:`);
  for (let i = 0; i < result2.energies.length; i++) {
    const E_eV = Schrodinger1DSolver.joulesToEV(result2.energies[i]);
    console.log(`  E_${i} = ${E_eV.toFixed(6)} eV`);
  }
  console.log("✓ Pöschl-Teller potential test PASSED");

} catch (error) {
  console.error("✗ Pöschl-Teller potential test FAILED:", error.message);
}

console.log("\n" + "=".repeat(60));
console.log("All tests completed!");
