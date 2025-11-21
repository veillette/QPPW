/**
 * Integration test to verify analytical double well solver is used by TwoWellsModel.
 */

import { TwoWellsModel } from "../src/two-wells/model/TwoWellsModel.js";
import QuantumConstants from "../src/common/model/QuantumConstants.js";

const { EV_TO_JOULES } = QuantumConstants;

console.log("Testing Double Well Analytical Solver Integration");
console.log("=".repeat(60));

// Create model instance
const model = new TwoWellsModel();

// Set parameters for a typical double well
model.wellWidthProperty.value = 5;  // 5 nm
model.wellDepthProperty.value = 0.3;  // 0.3 eV
model.wellSeparationProperty.value = 4;  // 4 nm separation
model.particleMassProperty.value = 0.067;  // GaAs effective mass

console.log("\nModel parameters:");
console.log(`  Well width: ${model.wellWidthProperty.value} nm`);
console.log(`  Well depth: ${model.wellDepthProperty.value} eV`);
console.log(`  Well separation: ${model.wellSeparationProperty.value} nm`);
console.log(`  Particle mass: ${model.particleMassProperty.value} m_e`);

// Get bound states
const boundStates = model.getBoundStates();

if (!boundStates) {
  console.error("\n✗ FAILED: No bound states returned");
  process.exit(1);
}

console.log(`\n✓ Found ${boundStates.energies.length} bound states`);

// Check that analytical method was used
if (boundStates.method !== "analytical") {
  console.error(`\n✗ FAILED: Expected analytical method, got ${boundStates.method}`);
  process.exit(1);
}

console.log(`✓ Method: ${boundStates.method}`);

// Verify energies are reasonable
console.log("\nBound state energies:");
for (let i = 0; i < boundStates.energies.length; i++) {
  const energyEV = boundStates.energies[i] / EV_TO_JOULES;
  console.log(`  State ${i}: ${energyEV.toFixed(6)} eV`);

  // All energies should be negative (bound states)
  if (energyEV >= 0) {
    console.error(`\n✗ FAILED: Bound state ${i} has non-negative energy`);
    process.exit(1);
  }

  // All energies should be greater than -wellDepth
  if (energyEV < -model.wellDepthProperty.value) {
    console.error(`\n✗ FAILED: Bound state ${i} energy below well depth`);
    process.exit(1);
  }
}

console.log("\n✓ All energies are valid bound states");

// Verify wavefunctions are normalized
console.log("\nNormalization check:");
const dx = boundStates.xGrid[1] - boundStates.xGrid[0];
for (let i = 0; i < boundStates.wavefunctions.length; i++) {
  const wf = boundStates.wavefunctions[i];
  const norm = wf.reduce((sum, val) => sum + val * val, 0) * dx;

  console.log(`  State ${i}: ${norm.toFixed(6)}`);

  if (Math.abs(norm - 1.0) > 0.01) {
    console.error(`\n✗ FAILED: State ${i} not normalized (norm = ${norm})`);
    process.exit(1);
  }
}

console.log("\n✓ All wavefunctions are properly normalized");

// Verify energy level method works
console.log("\nEnergy level access:");
for (let n = 1; n <= Math.min(3, boundStates.energies.length); n++) {
  const energy = model.getEnergyLevel(n);
  const expectedEnergy = boundStates.energies[n - 1] / EV_TO_JOULES;
  const diff = Math.abs(energy - expectedEnergy);

  console.log(`  n=${n}: ${energy.toFixed(6)} eV`);

  if (diff > 1e-6) {
    console.error(`\n✗ FAILED: getEnergyLevel(${n}) mismatch`);
    process.exit(1);
  }
}

console.log("\n✓ Energy level access working correctly");

console.log("\n" + "=".repeat(60));
console.log("Integration test PASSED ✓");
console.log("\nThe analytical solver is successfully integrated and");
console.log("being used as the default method for double well potential!");
