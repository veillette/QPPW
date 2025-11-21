/**
 * Test to verify parity alternation in double quantum well.
 *
 * For a symmetric double well, the bound states MUST alternate in parity:
 * - Ground state (n=0): EVEN
 * - 1st excited (n=1): ODD
 * - 2nd excited (n=2): EVEN
 * - 3rd excited (n=3): ODD
 * - etc.
 */

import { solveDoubleSquareWellAnalytical } from "../src/common/model/analytical-solutions/double-square-well.js";
import QuantumConstants from "../src/common/model/QuantumConstants.js";

const { ELECTRON_MASS, EV_TO_JOULES, NM_TO_M } = QuantumConstants;

/**
 * Detect parity by checking symmetry of wavefunction
 */
function detectParityFromWavefunction(wf: number[], xGrid: number[]): "even" | "odd" {
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
    throw new Error(`Wavefunction has neither clear even nor odd parity! sym=${relSymDiff.toFixed(4)}, anti=${relAntiDiff.toFixed(4)}`);
  }
}

console.log("Testing Parity Alternation in Double Quantum Well");
console.log("=".repeat(70));

// Test multiple parameter sets
const testCases = [
  {
    name: "Narrow barrier (strongly coupled)",
    wellWidth: 5 * NM_TO_M,
    wellDepth: 0.3 * EV_TO_JOULES,
    wellSeparation: 2 * NM_TO_M,
    numStates: 4,
  },
  {
    name: "Medium barrier",
    wellWidth: 5 * NM_TO_M,
    wellDepth: 0.4 * EV_TO_JOULES,
    wellSeparation: 4 * NM_TO_M,
    numStates: 5,
  },
  {
    name: "Wide barrier (weakly coupled)",
    wellWidth: 4 * NM_TO_M,
    wellDepth: 0.5 * EV_TO_JOULES,
    wellSeparation: 8 * NM_TO_M,
    numStates: 4,
  },
  {
    name: "Deep wells (many states)",
    wellWidth: 6 * NM_TO_M,
    wellDepth: 0.6 * EV_TO_JOULES,
    wellSeparation: 3 * NM_TO_M,
    numStates: 8,
  },
];

const mass = 0.067 * ELECTRON_MASS;
let allPassed = true;

for (const testCase of testCases) {
  console.log(`\n${testCase.name}`);
  console.log("-".repeat(70));
  console.log(`  Well width: ${testCase.wellWidth / NM_TO_M} nm`);
  console.log(`  Well depth: ${testCase.wellDepth / EV_TO_JOULES} eV`);
  console.log(`  Barrier width: ${testCase.wellSeparation / NM_TO_M} nm`);
  console.log(`  Requested states: ${testCase.numStates}`);

  const gridConfig = {
    xMin: -20 * NM_TO_M,
    xMax: 20 * NM_TO_M,
    numPoints: 600,
  };

  try {
    const result = solveDoubleSquareWellAnalytical(
      testCase.wellWidth,
      testCase.wellDepth,
      testCase.wellSeparation,
      mass,
      testCase.numStates,
      gridConfig
    );

    console.log(`  Found ${result.energies.length} bound states`);

    // Verify parity alternation
    const parities: ("even" | "odd")[] = [];
    for (let i = 0; i < result.energies.length; i++) {
      const parity = detectParityFromWavefunction(result.wavefunctions[i], result.xGrid);
      parities.push(parity);
    }

    console.log(`  Parities: ${parities.join(", ")}`);

    // Check that parities alternate: even, odd, even, odd, ...
    let alternationCorrect = true;
    for (let i = 0; i < parities.length; i++) {
      const expectedParity = (i % 2 === 0) ? "even" : "odd";
      if (parities[i] !== expectedParity) {
        console.log(`  ❌ FAIL: State ${i} has ${parities[i]} parity, expected ${expectedParity}`);
        alternationCorrect = false;
        allPassed = false;
      }
    }

    if (alternationCorrect) {
      console.log(`  ✓ PASS: All ${parities.length} states have correct alternating parity`);
    }

  } catch (error) {
    console.error(`  ❌ FAIL: ${error}`);
    allPassed = false;
  }
}

console.log("\n" + "=".repeat(70));
if (allPassed) {
  console.log("✓ All parity alternation tests PASSED!");
  process.exit(0);
} else {
  console.log("❌ Some parity alternation tests FAILED");
  process.exit(1);
}
