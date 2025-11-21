/**
 * Simple test to verify the analytical solver can be imported and used.
 */

import { solveDoubleSquareWellAnalytical } from "../src/common/model/analytical-solutions/double-square-well.js";
import QuantumConstants from "../src/common/model/QuantumConstants.js";

const { HBAR, ELECTRON_MASS, EV_TO_JOULES, NM_TO_M } = QuantumConstants;

console.log("Testing Analytical Solver Import");
console.log("=".repeat(60));

// Test parameters
const wellWidth = 5 * NM_TO_M;
const wellDepth = 0.3 * EV_TO_JOULES;
const wellSeparation = 4 * NM_TO_M;
const mass = 0.067 * ELECTRON_MASS;

const gridConfig = {
  xMin: -15 * NM_TO_M,
  xMax: 15 * NM_TO_M,
  numPoints: 300,
};

console.log("Calling solveDoubleSquareWellAnalytical...");

const result = solveDoubleSquareWellAnalytical(
  wellWidth,
  wellDepth,
  wellSeparation,
  mass,
  5,
  gridConfig
);

console.log(`\n✓ Function imported and executed successfully`);
console.log(`✓ Found ${result.energies.length} states`);
console.log(`✓ Method: ${result.method}`);
console.log(`✓ Grid has ${result.xGrid.length} points`);
console.log(`✓ Wavefunctions: ${result.wavefunctions.length}`);

if (result.method !== "analytical") {
  console.error(`✗ FAILED: Expected method='analytical', got '${result.method}'`);
  process.exit(1);
}

console.log("\n" + "=".repeat(60));
console.log("Import test PASSED ✓");
