/**
 * Comprehensive test for double quantum well with up to 20 states.
 * Tests parity alternation, edge behavior, and derivative continuity.
 */

import { solveDoubleSquareWellAnalytical } from "../src/common/model/analytical-solutions/double-square-well.js";
import QuantumConstants from "../src/common/model/QuantumConstants.js";

const { HBAR, ELECTRON_MASS, EV_TO_JOULES, NM_TO_M } = QuantumConstants;

/**
 * Detect parity by checking symmetry of wavefunction
 */
function detectParityFromWavefunction(wf: number[], xGrid: number[]): "even" | "odd" | "mixed" {
  let midIndex = 0;
  let minDist = Infinity;
  for (let j = 0; j < xGrid.length; j++) {
    const dist = Math.abs(xGrid[j]);
    if (dist < minDist) {
      minDist = dist;
      midIndex = j;
    }
  }

  let symDiff = 0, antiDiff = 0;
  let count = 0;

  for (let j = 1; j < Math.min(100, Math.floor(xGrid.length / 2)); j++) {
    const leftIdx = midIndex - j;
    const rightIdx = midIndex + j;

    if (leftIdx >= 0 && rightIdx < wf.length) {
      const valLeft = wf[leftIdx];
      const valRight = wf[rightIdx];

      symDiff += Math.abs(valLeft - valRight);
      antiDiff += Math.abs(valLeft + valRight);
      count++;
    }
  }

  symDiff /= count;
  antiDiff /= count;

  const maxMag = Math.max(...wf.map(v => Math.abs(v)));
  const relSymDiff = symDiff / maxMag;
  const relAntiDiff = antiDiff / maxMag;

  const tolerance = 0.1;

  if (relSymDiff < tolerance) {
    return "even";
  } else if (relAntiDiff < tolerance) {
    return "odd";
  } else {
    return "mixed";
  }
}

/**
 * Check derivative continuity at outer boundary
 */
function checkDerivativeContinuity(
  E: number,
  parity: "even" | "odd",
  Linner: number,
  L: number,
  V0: number,
  mass: number
): number {
  const k = Math.sqrt(2 * mass * E) / HBAR;
  const kappa = Math.sqrt(2 * mass * (V0 - E)) / HBAR;
  const alpha = kappa;

  let B: number, C: number;
  if (parity === "even") {
    B = Math.cosh(kappa * Linner);
    C = (kappa / k) * Math.sinh(kappa * Linner);
  } else {
    B = Math.sinh(kappa * Linner);
    C = (kappa / k) * Math.cosh(kappa * Linner);
  }

  const D = B * Math.cos(k * L) + C * Math.sin(k * L);
  const numerator = -k * B * Math.sin(k * L) + k * C * Math.cos(k * L);

  const derivInside = numerator / D;
  const derivOutside = -alpha;

  return Math.abs((derivInside - derivOutside) / derivOutside);
}

console.log("Comprehensive Test: Double Quantum Well with up to 20 States");
console.log("=".repeat(75));

// Test multiple parameter sets to find one that supports many states
const testCases = [
  {
    name: "Moderate wells (testing first 10 states)",
    wellWidth: 8 * NM_TO_M,
    wellDepth: 0.6 * EV_TO_JOULES,
    wellSeparation: 3 * NM_TO_M,
    numStates: 10,
  },
  {
    name: "Wider barrier weakly coupled (testing first 12 states)",
    wellWidth: 7 * NM_TO_M,
    wellDepth: 0.8 * EV_TO_JOULES,
    wellSeparation: 8 * NM_TO_M,
    numStates: 12,
  },
  {
    name: "Very wide and deep wells (testing first 20 states)",
    wellWidth: 20 * NM_TO_M,
    wellDepth: 0.5 * EV_TO_JOULES,
    wellSeparation: 4 * NM_TO_M,
    numStates: 20,
  },
];

const mass = 0.067 * ELECTRON_MASS;
let allPassed = true;

for (const testCase of testCases) {
  console.log(`\n${testCase.name}`);
  console.log("-".repeat(75));
  console.log(`  Well width: ${testCase.wellWidth / NM_TO_M} nm`);
  console.log(`  Well depth: ${testCase.wellDepth / EV_TO_JOULES} eV`);
  console.log(`  Barrier width: ${testCase.wellSeparation / NM_TO_M} nm`);
  console.log(`  Target states: ${testCase.numStates}`);

  const Linner = testCase.wellSeparation / 2;
  const _Louter = Linner + testCase.wellWidth;
  const L = testCase.wellWidth;

  const gridConfig = {
    xMin: -30 * NM_TO_M,
    xMax: 30 * NM_TO_M,
    numPoints: 1500,
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

    console.log(`  Found ${result.energies.length} bound states\n`);

    if (result.energies.length >= 10) {
      console.log("  State | Energy (eV) | Parity | Expected | Deriv Error");
      console.log("  " + "-".repeat(65));

      let parityErrors = 0;
      let derivErrors = 0;

      for (let i = 0; i < result.energies.length; i++) {
        const energyEV = result.energies[i] / EV_TO_JOULES;
        const parity = detectParityFromWavefunction(result.wavefunctions[i], result.xGrid);
        const expectedParity = (i % 2 === 0) ? "even" : "odd";
        const derivError = checkDerivativeContinuity(
          result.energies[i],
          parity as "even" | "odd",
          Linner,
          L,
          testCase.wellDepth,
          mass
        );

        const parityMatch = parity === expectedParity ? "✓" : "✗";
        const derivOk = derivError < 0.02 ? "✓" : "✗";

        if (i < 20) {  // Only print first 20
          console.log(
            `  ${i.toString().padStart(5)} | ${energyEV.toFixed(6)} | ${parity.padEnd(6)} | ${expectedParity.padEnd(8)} ${parityMatch} | ${(derivError * 100).toFixed(2)}% ${derivOk}`
          );
        }

        if (parity !== expectedParity) {
          parityErrors++;
        }
        if (derivError >= 0.02) {
          derivErrors++;
        }
      }

      console.log("\n  Summary:");
      console.log(`    Total states found: ${result.energies.length}`);
      console.log(`    Parity alternation: ${parityErrors === 0 ? "PASS ✓" : `FAIL ✗ (${parityErrors} errors)`}`);
      console.log(`    Derivative continuity: ${derivErrors === 0 ? "PASS ✓" : `FAIL ✗ (${derivErrors} errors)`}`);

      if (parityErrors > 0 || derivErrors > 0) {
        allPassed = false;
      }
    } else {
      console.log(`  ⚠️  Only ${result.energies.length} states found (less than 10)`);
      console.log(`      Parameter set doesn't support enough bound states`);
    }

  } catch (error) {
    console.error(`  ❌ FAIL: ${error}`);
    allPassed = false;
  }
}

console.log("\n" + "=".repeat(75));
if (allPassed) {
  console.log("✓ All comprehensive tests PASSED!");
  process.exit(0);
} else {
  console.log("❌ Some tests FAILED");
  process.exit(1);
}
