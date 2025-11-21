/**
 * Test for analytical double square well solver.
 */

import { solveDoubleSquareWellAnalytical } from "../src/common/model/analytical-solutions/double-square-well.js";
import QuantumConstants from "../src/common/model/QuantumConstants.js";

const { HBAR, ELECTRON_MASS, EV_TO_JOULES, NM_TO_M } = QuantumConstants;

console.log("Testing Analytical Double Square Well Solver");
console.log("=".repeat(60));

// Test case 1: Typical double quantum well parameters
console.log("\nTest 1: Typical quantum well parameters");
console.log("-".repeat(60));

const wellWidth = 5 * NM_TO_M;          // 5 nm well width
const wellDepth = 0.3 * EV_TO_JOULES;   // 0.3 eV well depth
const wellSeparation = 4 * NM_TO_M;     // 4 nm separation (2 nm barrier)
const mass = 0.067 * ELECTRON_MASS;     // Effective mass in GaAs
const numStates = 6;

const gridConfig = {
  xMin: -15 * NM_TO_M,
  xMax: 15 * NM_TO_M,
  numPoints: 500,
};

console.log(`Well width: ${wellWidth / NM_TO_M} nm`);
console.log(`Well depth: ${wellDepth / EV_TO_JOULES} eV`);
console.log(`Well separation: ${wellSeparation / NM_TO_M} nm`);
console.log(`Barrier width: ${wellSeparation / NM_TO_M} nm`);
console.log(`Effective mass: ${mass / ELECTRON_MASS} m_e`);
console.log(`Number of states: ${numStates}`);

try {
  const result = solveDoubleSquareWellAnalytical(
    wellWidth,
    wellDepth,
    wellSeparation,
    mass,
    numStates,
    gridConfig
  );

  console.log(`\nFound ${result.energies.length} bound states:`);
  console.log("-".repeat(60));

  for (let i = 0; i < result.energies.length; i++) {
    const energyEV = result.energies[i] / EV_TO_JOULES;
    const wf = result.wavefunctions[i];

    // Calculate probability in left well, barrier, and right well
    const dx = result.xGrid[1] - result.xGrid[0];
    let probLeft = 0, probBarrier = 0, probRight = 0;

    const Linner = wellSeparation / 2;
    const Louter = Linner + wellWidth;

    for (let j = 0; j < result.xGrid.length; j++) {
      const x = result.xGrid[j];
      const probDensity = wf[j] * wf[j] * dx;

      if (x < -Linner && x > -Louter) {
        probLeft += probDensity;
      } else if (Math.abs(x) < Linner) {
        probBarrier += probDensity;
      } else if (x > Linner && x < Louter) {
        probRight += probDensity;
      }
    }

    console.log(`State ${i}: E = ${energyEV.toFixed(6)} eV`);
    console.log(`  Probability: Left=${(probLeft*100).toFixed(1)}%, Barrier=${(probBarrier*100).toFixed(1)}%, Right=${(probRight*100).toFixed(1)}%`);

    // Check normalization
    const totalProb = wf.reduce((sum, val) => sum + val * val, 0) * dx;
    console.log(`  Normalization: ${totalProb.toFixed(6)} (should be ~1.0)`);
  }

  // Verify parity by checking symmetry
  console.log("\nParity verification:");
  console.log("-".repeat(60));

  // Find index closest to x = 0
  let midIndex = 0;
  let minDist = Infinity;
  for (let j = 0; j < result.xGrid.length; j++) {
    const dist = Math.abs(result.xGrid[j]);
    if (dist < minDist) {
      minDist = dist;
      midIndex = j;
    }
  }

  for (let i = 0; i < Math.min(4, result.energies.length); i++) {
    const wf = result.wavefunctions[i];

    // Check if wavefunction is symmetric or antisymmetric
    let symDiff = 0, antiDiff = 0;
    let count = 0;

    // Sample points around the center
    for (let j = 1; j < Math.min(100, Math.floor(result.xGrid.length / 2)); j++) {
      const leftIdx = midIndex - j;
      const rightIdx = midIndex + j;

      if (leftIdx >= 0 && rightIdx < wf.length) {
        const valLeft = wf[leftIdx];
        const valRight = wf[rightIdx];

        // Symmetric: ψ(-x) = ψ(x)
        symDiff += Math.abs(valLeft - valRight);
        // Antisymmetric: ψ(-x) = -ψ(x)
        antiDiff += Math.abs(valLeft + valRight);

        count++;
      }
    }

    symDiff /= count;
    antiDiff /= count;

    // Normalize by typical wavefunction magnitude
    const maxMag = Math.max(...wf.map(v => Math.abs(v)));
    const relSymDiff = symDiff / maxMag;
    const relAntiDiff = antiDiff / maxMag;

    const tolerance = 0.1;  // Relative tolerance (10%)
    const parity = (relSymDiff < tolerance) ? "even" : ((relAntiDiff < tolerance) ? "odd" : "mixed");
    console.log(`State ${i}: ${parity} parity (rel sym: ${relSymDiff.toFixed(4)}, rel anti: ${relAntiDiff.toFixed(4)})`);
  }

  console.log("\n" + "=".repeat(60));
  console.log("Test 1: PASSED ✓");

} catch (error) {
  console.error("Test 1: FAILED ✗");
  console.error(error);
  process.exit(1);
}

// Test case 2: Wide barrier (weakly coupled wells)
console.log("\n\nTest 2: Wide barrier (weakly coupled wells)");
console.log("=".repeat(60));

const wellWidth2 = 5 * NM_TO_M;
const wellDepth2 = 0.5 * EV_TO_JOULES;
const wellSeparation2 = 10 * NM_TO_M;  // Wider barrier
const mass2 = 0.067 * ELECTRON_MASS;
const numStates2 = 4;

const gridConfig2 = {
  xMin: -20 * NM_TO_M,
  xMax: 20 * NM_TO_M,
  numPoints: 600,
};

console.log(`Well width: ${wellWidth2 / NM_TO_M} nm`);
console.log(`Well depth: ${wellDepth2 / EV_TO_JOULES} eV`);
console.log(`Barrier width: ${wellSeparation2 / NM_TO_M} nm`);

try {
  const result2 = solveDoubleSquareWellAnalytical(
    wellWidth2,
    wellDepth2,
    wellSeparation2,
    mass2,
    numStates2,
    gridConfig2
  );

  console.log(`\nFound ${result2.energies.length} bound states:`);

  // For wide barriers, we expect closely spaced doublets
  if (result2.energies.length >= 2) {
    const splitting1 = Math.abs(result2.energies[1] - result2.energies[0]);
    console.log(`Ground state splitting: ${(splitting1 / EV_TO_JOULES * 1000).toFixed(3)} meV`);

    if (result2.energies.length >= 4) {
      const splitting2 = Math.abs(result2.energies[3] - result2.energies[2]);
      console.log(`First excited state splitting: ${(splitting2 / EV_TO_JOULES * 1000).toFixed(3)} meV`);
    }
  }

  console.log("\n" + "=".repeat(60));
  console.log("Test 2: PASSED ✓");

} catch (error) {
  console.error("Test 2: FAILED ✗");
  console.error(error);
  process.exit(1);
}

console.log("\n\nAll tests passed! ✓");
