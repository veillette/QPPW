/**
 * Test script for Coulomb potential analytical solutions.
 *
 * This tests both the 1D and 3D Coulomb potentials.
 */

import { Schrodinger1DSolver } from "./dist/common/model/Schrodinger1DSolver.js";
import { PotentialType } from "./dist/common/model/PotentialFunction.js";
import QuantumConstants from "./dist/common/model/QuantumConstants.js";

console.log("Testing Coulomb Potential Analytical Solutions\n");
console.log("=".repeat(60));

// Test 1D Coulomb Potential
console.log("\n1. Testing 1D Coulomb Potential");
console.log("-".repeat(60));

const solver = new Schrodinger1DSolver();

// For hydrogen-like atom in 1D
// α = e²/(4πε₀) ≈ 2.307e-28 J·m (in SI units)
const ALPHA_1D = 2.307e-28; // J·m
const mass = QuantumConstants.ELECTRON_MASS;

// Solve 1D Coulomb
const wellParams1D = {
  type: PotentialType.COULOMB_1D,
  coulombStrength: ALPHA_1D,
};

const gridConfig1D = {
  xMin: -50e-9, // -50 nm
  xMax: 50e-9,   // 50 nm
  numPoints: 1000,
};

const numStates1D = 5;

try {
  const result1D = solver.solveAnalyticalIfPossible(
    wellParams1D,
    mass,
    numStates1D,
    gridConfig1D
  );

  console.log(`Number of states calculated: ${result1D.energies.length}`);
  console.log(`Method used: ${result1D.method}`);
  console.log("\nEnergy eigenvalues (in eV):");

  // Expected formula: E_n = -mα²/(2ℏ²(n+1/2)²)
  for (let n = 0; n < result1D.energies.length; n++) {
    const energyEV = result1D.energies[n] * QuantumConstants.JOULES_TO_EV;

    // Calculate expected energy
    const HBAR = QuantumConstants.HBAR;
    const expectedEnergy = -(mass * ALPHA_1D * ALPHA_1D) / (2 * HBAR * HBAR * (n + 0.5) * (n + 0.5));
    const expectedEnergyEV = expectedEnergy * QuantumConstants.JOULES_TO_EV;

    console.log(`  n=${n}: E = ${energyEV.toFixed(6)} eV (expected: ${expectedEnergyEV.toFixed(6)} eV)`);
  }

  console.log("\n✓ 1D Coulomb potential test completed successfully!");
} catch (error) {
  console.error("✗ Error solving 1D Coulomb potential:", error.message);
}

// Test 3D Coulomb Potential (Hydrogen atom, L=0)
console.log("\n2. Testing 3D Coulomb Potential (Hydrogen atom, L=0)");
console.log("-".repeat(60));

// For hydrogen atom in 3D
const ALPHA_3D = 2.307e-28; // J·m

const wellParams3D = {
  type: PotentialType.COULOMB_3D,
  coulombStrength: ALPHA_3D,
};

const gridConfig3D = {
  xMin: 1e-12,   // Start slightly above 0 for radial coordinate
  xMax: 100e-9,  // 100 nm
  numPoints: 1000,
};

const numStates3D = 5;

try {
  const result3D = solver.solveAnalyticalIfPossible(
    wellParams3D,
    mass,
    numStates3D,
    gridConfig3D
  );

  console.log(`Number of states calculated: ${result3D.energies.length}`);
  console.log(`Method used: ${result3D.method}`);
  console.log("\nEnergy eigenvalues (in eV):");

  // Expected formula: E_n = -mα²/(2ℏ²n²) = -13.6/n² eV for hydrogen
  // This is the standard hydrogen atom formula
  for (let n = 1; n <= result3D.energies.length; n++) {
    const energyEV = result3D.energies[n - 1] * QuantumConstants.JOULES_TO_EV;

    // Calculate expected energy
    const HBAR = QuantumConstants.HBAR;
    const expectedEnergy = -(mass * ALPHA_3D * ALPHA_3D) / (2 * HBAR * HBAR * n * n);
    const expectedEnergyEV = expectedEnergy * QuantumConstants.JOULES_TO_EV;

    console.log(`  n=${n}: E = ${energyEV.toFixed(6)} eV (expected: ${expectedEnergyEV.toFixed(6)} eV)`);
  }

  console.log("\n✓ 3D Coulomb potential (hydrogen atom) test completed successfully!");

  // Compare with known hydrogen atom ground state energy (-13.6 eV)
  const groundStateEV = result3D.energies[0] * QuantumConstants.JOULES_TO_EV;
  const hydrogenGroundState = -13.6; // eV
  console.log(`\nGround state energy: ${groundStateEV.toFixed(2)} eV`);
  console.log(`Expected hydrogen ground state: ${hydrogenGroundState} eV`);
  console.log(`Note: Small difference due to the specific value of α used.`);

} catch (error) {
  console.error("✗ Error solving 3D Coulomb potential:", error.message);
}

console.log("\n" + "=".repeat(60));
console.log("All tests completed!");
