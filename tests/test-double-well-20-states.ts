/**
 * Comprehensive test for double quantum well with 20 states.
 * Tests parity alternation and edge behavior.
 */

import { solveDoubleSquareWellAnalytical } from "../src/common/model/analytical-solutions/double-square-well.js";
import QuantumConstants from "../src/common/model/QuantumConstants.js";

const { ELECTRON_MASS, EV_TO_JOULES, NM_TO_M } = QuantumConstants;

/**
 * Detect parity by checking symmetry of wavefunction
 */
function detectParityFromWavefunction(wf: number[], xGrid: number[]): "even" | "odd" | "mixed" {
  // Find index closest to x = 0
  let midIndex = 0;
  let minDist = Infinity;
  for (let j = 0; j < xGrid.length; j++) {
    const dist = Math.abs(xGrid[j]);
    if (dist < minDist) {
      minDist = dist;
      midIndex = j;
    }
  }

  // Check if wavefunction is symmetric or antisymmetric
  let symDiff = 0, antiDiff = 0;
  let count = 0;

  // Sample points around the center
  for (let j = 1; j < Math.min(100, Math.floor(xGrid.length / 2)); j++) {
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

  if (relSymDiff < tolerance) {
    return "even";
  } else if (relAntiDiff < tolerance) {
    return "odd";
  } else {
    return "mixed";
  }
}

/**
 * Check edge behavior of wavefunction
 */
function checkEdgeBehavior(wf: number[], xGrid: number[], Louter: number): {
  hasNode: boolean;
  edgeValue: number;
  penetrationOk: boolean;
  edgeSlope: number;
} {
  // Find index closest to x = Louter (outer edge of right well)
  let edgeIdx = 0;
  let minDist = Infinity;
  for (let j = 0; j < xGrid.length; j++) {
    const dist = Math.abs(xGrid[j] - Louter);
    if (dist < minDist) {
      minDist = dist;
      edgeIdx = j;
    }
  }

  const edgeValue = Math.abs(wf[edgeIdx]);

  // Check a few points outside the well to see if there's exponential decay
  const outsidePoints = [];
  for (let j = edgeIdx; j < Math.min(edgeIdx + 10, xGrid.length); j++) {
    if (xGrid[j] > Louter) {
      outsidePoints.push(Math.abs(wf[j]));
    }
  }

  // Check if there's a node right at the edge (which would be wrong)
  const hasNode = edgeValue < 1e-10;

  // Check if wavefunction decays outside (proper penetration)
  let penetrationOk = true;
  if (outsidePoints.length >= 2) {
    // Values should be decreasing
    for (let i = 1; i < outsidePoints.length; i++) {
      if (outsidePoints[i] > outsidePoints[i-1] * 1.1) {
        penetrationOk = false;
        break;
      }
    }
  }

  // Calculate slope at edge
  let edgeSlope = 0;
  if (edgeIdx + 1 < xGrid.length) {
    const dx = xGrid[edgeIdx + 1] - xGrid[edgeIdx];
    edgeSlope = (wf[edgeIdx + 1] - wf[edgeIdx]) / dx;
  }

  return { hasNode, edgeValue, penetrationOk, edgeSlope };
}

console.log("Testing Double Quantum Well with 20 States");
console.log("=".repeat(70));

// Test with deep wells to support 20 states
const wellWidth = 10 * NM_TO_M;         // 10 nm well width (wider to support more states)
const wellDepth = 2.0 * EV_TO_JOULES;   // 2.0 eV well depth (deeper to support more states)
const wellSeparation = 3 * NM_TO_M;     // 3 nm separation
const mass = 0.067 * ELECTRON_MASS;     // Effective mass in GaAs
const numStates = 20;

const Linner = wellSeparation / 2;
const Louter = Linner + wellWidth;

const gridConfig = {
  xMin: -25 * NM_TO_M,
  xMax: 25 * NM_TO_M,
  numPoints: 1000,
};

console.log(`Well width: ${wellWidth / NM_TO_M} nm`);
console.log(`Well depth: ${wellDepth / EV_TO_JOULES} eV`);
console.log(`Well separation: ${wellSeparation / NM_TO_M} nm`);
console.log(`Effective mass: ${mass / ELECTRON_MASS} m_e`);
console.log(`Requesting ${numStates} states\n`);

try {
  const result = solveDoubleSquareWellAnalytical(
    wellWidth,
    wellDepth,
    wellSeparation,
    mass,
    numStates,
    gridConfig
  );

  console.log(`Found ${result.energies.length} bound states\n`);

  // Check parity alternation for all states
  const parities: ("even" | "odd" | "mixed")[] = [];
  let parityErrors = 0;

  console.log("State | Energy (eV) | Parity  | Expected | Edge Value | Node? | Penetration");
  console.log("-".repeat(70));

  for (let i = 0; i < result.energies.length; i++) {
    const energyEV = result.energies[i] / EV_TO_JOULES;
    const parity = detectParityFromWavefunction(result.wavefunctions[i], result.xGrid);
    const expectedParity = (i % 2 === 0) ? "even" : "odd";
    const edgeInfo = checkEdgeBehavior(result.wavefunctions[i], result.xGrid, Louter);

    parities.push(parity);

    const parityMatch = parity === expectedParity ? "✓" : "✗";
    const nodeWarning = edgeInfo.hasNode ? "YES!" : "no";
    const penetrationStatus = edgeInfo.penetrationOk ? "OK" : "BAD!";

    console.log(
      `${i.toString().padStart(5)} | ${energyEV.toFixed(6)} | ${parity.padEnd(7)} | ${expectedParity.padEnd(8)} ${parityMatch} | ` +
      `${edgeInfo.edgeValue.toExponential(2)} | ${nodeWarning.padEnd(5)} | ${penetrationStatus}`
    );

    if (parity !== expectedParity) {
      parityErrors++;
    }

    if (edgeInfo.hasNode) {
      console.log(`      WARNING: State ${i} has a node at the edge!`);
    }

    if (!edgeInfo.penetrationOk) {
      console.log(`      WARNING: State ${i} has improper penetration into barrier!`);
    }
  }

  console.log("\n" + "=".repeat(70));
  console.log(`Parity alternation: ${parityErrors === 0 ? "PASS ✓" : `FAIL ✗ (${parityErrors} errors)`}`);

  if (parityErrors > 0) {
    console.log("\nParity sequence:");
    console.log(`Found:    ${parities.join(", ")}`);
    const expected = Array.from({ length: parities.length }, (_, i) => i % 2 === 0 ? "even" : "odd");
    console.log(`Expected: ${expected.join(", ")}`);
  }

  process.exit(parityErrors > 0 ? 1 : 0);

} catch (error) {
  console.error("Test FAILED:");
  console.error(error);
  process.exit(1);
}
