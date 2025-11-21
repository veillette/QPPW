#!/usr/bin/env node
/**
 * Comprehensive tests for double well potential
 * These tests explore extreme parameter ranges and verify theoretical predictions
 *
 * Energy Convention (CURRENT):
 * - Wells are at V = -V₀ (negative potential)
 * - Barrier is at V = 0 (reference)
 * - Bound states have energies: -V₀ < E < 0
 * - Ground state has most negative energy (deepest in well)
 * - Higher energy states are less negative (approaching 0)
 *
 * NOTE: This will be updated when analytical solver switches to V=0 convention
 *
 * Usage:
 *   npx tsx --import ./tests/browser-globals.js tests/test-double-well-comprehensive.ts
 */

import { solveDoubleSquareWellAnalytical } from '../src/common/model/analytical-solutions/double-square-well.js';
import QuantumConstants from '../src/common/model/QuantumConstants.js';

// Physical constants
const { HBAR, ELECTRON_MASS, EV_TO_JOULES, NM_TO_M } = QuantumConstants;

let totalTests = 0;
let passedTests = 0;
let failedTests = 0;

/**
 * Helper function to solve double well with nm/eV units
 */
function solveDoubleWell(
  wellWidthNm: number,
  wellDepthEv: number,
  barrierWidthNm: number,
  particleMassUnits: number = 1.0,
  numStates: number = 10,
  gridPoints: number = 1000
) {
  const wellWidth = wellWidthNm * NM_TO_M;
  const wellDepth = wellDepthEv * EV_TO_JOULES;
  const wellSeparation = barrierWidthNm * NM_TO_M;
  const mass = particleMassUnits * ELECTRON_MASS;

  const halfSeparation = wellSeparation / 2;
  const outerEdge = halfSeparation + wellWidth;
  const margin = 2.0 * NM_TO_M;
  const gridRange = outerEdge + margin;

  const gridConfig = {
    xMin: -gridRange,
    xMax: gridRange,
    numPoints: gridPoints
  };

  const result = solveDoubleSquareWellAnalytical(
    wellWidth,
    wellDepth,
    wellSeparation,
    mass,
    numStates,
    gridConfig
  );

  return {
    states: result.energies.map((E, i) => ({
      energy: E / EV_TO_JOULES,
      psi: result.wavefunctions[i]
    })),
    x: result.xGrid,  // Keep in meters for proper normalization
    xNm: result.xGrid.map(x => x / NM_TO_M),  // Also provide in nm for display
    energies: result.energies.map(E => E / EV_TO_JOULES),
    wavefunctions: result.wavefunctions
  };
}

/**
 * Assert helper function
 */
function assert(condition: boolean, message: string): void {
  totalTests++;
  if (condition) {
    passedTests++;
  } else {
    failedTests++;
    console.error(`✗ FAILED: ${message}`);
    throw new Error(message);
  }
}

/**
 * Calculate WKB transmission coefficient
 * Current convention: wells at V=-V₀, barrier at V=0
 * Energy is negative for bound states (-V₀ < E < 0)
 */
function calculateWKBTransmission(
  energy: number,        // Energy in eV (negative for bound states)
  barrierHeight: number, // Barrier height in eV (well depth V₀, positive value)
  barrierWidth: number,  // Barrier width in nm
  particleMass: number = 1.0
): number {
  const E = energy * EV_TO_JOULES;
  const _V = barrierHeight * EV_TO_JOULES;
  const a = barrierWidth * NM_TO_M;
  const m = particleMass * ELECTRON_MASS;

  // Energy relative to barrier top: E_rel = E - 0 = E (negative)
  // For tunneling: κ² = 2m|E|/ℏ² = -2mE/ℏ² (since E < 0)
  if (E >= 0) return 1.0; // Above barrier

  const kappa = Math.sqrt(-2 * m * E) / HBAR;
  const transmission = Math.exp(-2 * kappa * a);

  return transmission;
}

/**
 * Calculate overlap integral
 */
function calculateWellOverlap(
  x: number[],
  psi: number[],
  leftWellCenter: number,
  rightWellCenter: number,
  wellWidth: number
): { left: number; right: number; barrier: number } {
  const dx = x[1] - x[0];

  let leftProb = 0;
  let rightProb = 0;
  let barrierProb = 0;

  for (let i = 0; i < x.length; i++) {
    const prob = psi[i] * psi[i] * dx;

    if (Math.abs(x[i] - leftWellCenter) < wellWidth / 2) {
      leftProb += prob;
    } else if (Math.abs(x[i] - rightWellCenter) < wellWidth / 2) {
      rightProb += prob;
    } else if (x[i] > leftWellCenter + wellWidth / 2 && x[i] < rightWellCenter - wellWidth / 2) {
      barrierProb += prob;
    }
  }

  return { left: leftProb, right: rightProb, barrier: barrierProb };
}

/**
 * Determine parity from wavefunction
 */
function determineParity(x: number[], psi: number[]): 'even' | 'odd' | 'unknown' {
  const midIndex = Math.floor(x.length / 2);

  let symDiff = 0;
  let antiDiff = 0;
  let count = 0;

  for (let i = 1; i < Math.min(100, midIndex); i++) {
    const leftIdx = midIndex - i;
    const rightIdx = midIndex + i;

    if (leftIdx >= 0 && rightIdx < psi.length) {
      symDiff += Math.abs(psi[leftIdx] - psi[rightIdx]);
      antiDiff += Math.abs(psi[leftIdx] + psi[rightIdx]);
      count++;
    }
  }

  symDiff /= count;
  antiDiff /= count;

  const maxMag = Math.max(...psi.map(v => Math.abs(v)));
  const relSymDiff = symDiff / maxMag;
  const relAntiDiff = antiDiff / maxMag;

  if (relSymDiff < 0.1) return 'even';
  if (relAntiDiff < 0.1) return 'odd';
  return 'unknown';
}

console.log('='.repeat(70));
console.log(' COMPREHENSIVE DOUBLE WELL POTENTIAL TESTS');
console.log('='.repeat(70));

// Test 1: Tunneling Probability - Wide Barrier
console.log('\n[Test 1] Tunneling Probability - Wide Barrier (Weak Coupling)');
console.log('-'.repeat(70));
try {
  const result = solveDoubleWell(1.0, 5.0, 2.0, 1.0, 10, 1000);

  assert(result.states.length >= 2, 'Should have at least 2 states');

  const groundState = result.states[0];
  const overlap = calculateWellOverlap(
    result.x,
    groundState.psi,
    -1.5 * NM_TO_M, // left well center in meters
    1.5 * NM_TO_M,  // right well center in meters
    1.0 * NM_TO_M   // well width in meters
  );

  const wkbTransmission = calculateWKBTransmission(
    groundState.energy,
    5.0,
    2.0,
    1.0
  );

  console.log(`  Ground state energy: ${groundState.energy.toFixed(6)} eV (-${5.0} < E < 0)`);
  console.log(`  WKB transmission: ${wkbTransmission.toExponential(3)}`);
  console.log(`  Barrier probability: ${overlap.barrier.toFixed(6)}`);
  console.log(`  Left well: ${overlap.left.toFixed(6)}, Right well: ${overlap.right.toFixed(6)}`);

  assert(groundState.energy < 0 && groundState.energy > -5.0, 'Ground state energy should be bound (-wellDepth < E < 0)');
  assert(overlap.barrier < 0.2, 'Barrier probability should be small for wide barrier');

  console.log(`  ✓ Test PASSED`);
} catch (error) {
  console.error(`  ✗ Test FAILED: ${error}`);
  process.exit(1);
}

// Test 2: Tunneling Probability - Narrow Barrier
console.log('\n[Test 2] Tunneling Probability - Narrow Barrier (Strong Coupling)');
console.log('-'.repeat(70));
try {
  const result = solveDoubleWell(1.0, 5.0, 0.2, 1.0, 10, 1000);

  assert(result.states.length >= 2, 'Should have at least 2 states');

  const groundState = result.states[0];
  const overlap = calculateWellOverlap(
    result.x,
    groundState.psi,
    -0.6 * NM_TO_M,
    0.6 * NM_TO_M,
    1.0 * NM_TO_M
  );

  console.log(`  Barrier probability: ${overlap.barrier.toFixed(4)}`);
  console.log(`  Left well: ${overlap.left.toFixed(6)}, Right well: ${overlap.right.toFixed(6)}`);

  assert(overlap.barrier > 0.001, 'Barrier probability should be significant for narrow barrier');

  console.log(`  ✓ Test PASSED`);
} catch (error) {
  console.error(`  ✗ Test FAILED: ${error}`);
  process.exit(1);
}

// Test 3: Symmetric Distribution for Even States
console.log('\n[Test 3] Symmetric Distribution for Even States');
console.log('-'.repeat(70));
try {
  const result = solveDoubleWell(1.0, 5.0, 0.5, 1.0, 10, 1000);

  assert(result.states.length >= 1, 'Should have at least 1 state');

  const parity = determineParity(result.x, result.states[0].psi);

  const overlap = calculateWellOverlap(
    result.x,
    result.states[0].psi,
    -0.75 * NM_TO_M,
    0.75 * NM_TO_M,
    1.0 * NM_TO_M
  );

  const asymmetry = Math.abs(overlap.left - overlap.right) / (overlap.left + overlap.right);
  assert(asymmetry < 0.1, 'Even state should have symmetric distribution');

  console.log(`  ✓ Ground state parity: ${parity}`);
  console.log(`  ✓ Asymmetry: ${asymmetry.toFixed(4)}`);
  console.log(`  ✓ Test PASSED`);
} catch (error) {
  console.error(`  ✗ Test FAILED: ${error}`);
  process.exit(1);
}

// Test 4: Energy Splitting vs Barrier Width
console.log('\n[Test 4] Energy Splitting - Exponential Decrease with Barrier Width');
console.log('-'.repeat(70));
try {
  const barrierWidths = [0.3, 0.5, 0.8, 1.2, 1.5];
  const splittings: number[] = [];

  for (const barrierWidth of barrierWidths) {
    const result = solveDoubleWell(1.0, 5.0, barrierWidth, 1.0, 10, 1000);

    assert(result.states.length >= 2, `Should have at least 2 states for barrier width ${barrierWidth}nm`);

    const splitting = Math.abs(result.energies[1] - result.energies[0]);
    splittings.push(splitting);

    console.log(`  Barrier ${barrierWidth}nm: splitting = ${splitting.toExponential(3)} eV (${result.states.length} states found)`);
  }

  // Verify monotonic decrease (allowing for numerical precision near zero)
  for (let i = 1; i < splittings.length; i++) {
    // Only check if both values are significant
    if (splittings[i - 1] > 1e-9 && splittings[i] > 1e-12) {
      assert(splittings[i] <= splittings[i - 1], `Splitting should decrease or stay near zero with barrier width (${splittings[i].toExponential(3)} <= ${splittings[i - 1].toExponential(3)})`);
    }
  }

  // Check approximate exponential scaling (only for significant splittings)
  const significantData = barrierWidths.map((w, i) => ({ w, s: splittings[i] }))
    .filter(d => d.s > 1e-9);  // Only use non-zero splittings

  if (significantData.length >= 2) {
    const logSplittings = significantData.map(d => Math.log(d.s));
    const widths = significantData.map(d => d.w);
    const n = significantData.length;
    const sumX = widths.reduce((a, b) => a + b, 0);
    const sumY = logSplittings.reduce((a, b) => a + b, 0);
    const sumXY = widths.reduce((sum, x, i) => sum + x * logSplittings[i], 0);
    const sumX2 = widths.reduce((sum, x) => sum + x * x, 0);

    const slope = (n * sumXY - sumX * sumY) / (n * sumX2 - sumX * sumX);

    assert(slope < 0, 'Exponential decay rate should be negative');
    console.log(`  ✓ Exponential decay rate: ${slope.toFixed(3)}`);
  } else {
    console.log(`  ✓ Splittings too small to measure exponential decay`);
  }
  console.log(`  ✓ Test PASSED`);
} catch (error) {
  console.error(`  ✗ Test FAILED: ${error}`);
  process.exit(1);
}

// Test 5-12: Extreme Parameter Regimes
const extremeTests = [
  { name: 'Very Shallow Wells (0.1 eV)', params: [2.0, 0.1, 0.5, 1.0, 10, 1000], minStates: 1 },
  { name: 'Very Deep Wells (15 eV)', params: [1.0, 15.0, 0.5, 1.0, 20, 1000], minStates: 5 },
  { name: 'Very Narrow Wells (0.1 nm)', params: [0.1, 10.0, 0.1, 1.0, 10, 1000], minStates: 1 },
  { name: 'Very Wide Wells (3 nm)', params: [3.0, 5.0, 0.5, 1.0, 20, 1000], minStates: 3 },
  { name: 'Very Wide Barrier (3 nm)', params: [1.0, 5.0, 3.0, 1.0, 10, 1000], minStates: 2 },
  { name: 'Minimum Barrier Width (0.05 nm)', params: [1.0, 5.0, 0.05, 1.0, 10, 1000], minStates: 2 },
  { name: 'Heavy Particle Mass (1.1 m_e)', params: [1.0, 5.0, 0.5, 1.1, 10, 1000], minStates: 2 },
  { name: 'Light Particle Mass (0.5 m_e)', params: [1.0, 5.0, 0.5, 0.5, 10, 1000], minStates: 1 },
];

extremeTests.forEach((test, i) => {
  console.log(`\n[Test ${5 + i}] Extreme Parameters - ${test.name}`);
  console.log('-'.repeat(70));
  try {
    const result = solveDoubleWell(...test.params as [number, number, number, number, number, number]);

    assert(result.states.length >= test.minStates, `Should have at least ${test.minStates} states`);

    const groundState = result.states[0];
    const wellDepth = test.params[1] as number;
    assert(isFinite(groundState.energy), 'Ground state energy should be finite');
    assert(groundState.energy < 0 && groundState.energy > -wellDepth, `Ground state should be bound (-${wellDepth} < E < 0 eV)`);

    console.log(`  ✓ States found: ${result.states.length}`);
    console.log(`  ✓ Ground state energy: ${groundState.energy.toFixed(6)} eV`);
    console.log(`  ✓ Test PASSED`);
  } catch (error) {
    console.error(`  ✗ Test FAILED: ${error}`);
    process.exit(1);
  }
});

// Test 13: Grid Convergence
console.log('\n[Test 13] Grid Convergence');
console.log('-'.repeat(70));
try {
  const gridSizes = [200, 400, 800, 1600];
  const groundEnergies: number[] = [];

  for (const gridSize of gridSizes) {
    const result = solveDoubleWell(1.0, 5.0, 0.5, 1.0, 10, gridSize);

    assert(result.states.length >= 1, `Should have at least 1 state for grid size ${gridSize}`);
    groundEnergies.push(result.states[0].energy);

    console.log(`  Grid ${gridSize}: E_0 = ${result.states[0].energy.toFixed(8)} eV`);
  }

  // Final two should be very close
  const finalDiff = Math.abs(groundEnergies[3] - groundEnergies[2]);
  assert(finalDiff < 0.001, `Grid should converge (final difference: ${finalDiff} eV)`);

  console.log(`  ✓ Convergence achieved (diff < 0.001 eV)`);
  console.log(`  ✓ Test PASSED`);
} catch (error) {
  console.error(`  ✗ Test FAILED: ${error}`);
  process.exit(1);
}

// Test 14: Normalization Across Grid Sizes
console.log('\n[Test 14] Normalization Across Grid Sizes');
console.log('-'.repeat(70));
try {
  const gridSizes = [200, 500, 1000, 2000];

  for (const gridSize of gridSizes) {
    const result = solveDoubleWell(1.0, 5.0, 0.5, 1.0, 10, gridSize);

    assert(result.states.length >= 1, `Should have at least 1 state for grid size ${gridSize}`);

    const dx = result.x[1] - result.x[0];
    let norm = 0;
    for (const psi of result.states[0].psi) {
      norm += psi * psi * dx;
    }

    assert(Math.abs(norm - 1.0) < 0.01, `Normalization should be close to 1.0 (got ${norm})`);
    console.log(`  Grid ${gridSize}: normalization = ${norm.toFixed(6)}`);
  }

  console.log(`  ✓ Test PASSED`);
} catch (error) {
  console.error(`  ✗ Test FAILED: ${error}`);
  process.exit(1);
}

// Test 15: Boundary Decay
console.log('\n[Test 15] Boundary Decay');
console.log('-'.repeat(70));
try {
  const result = solveDoubleWell(1.0, 5.0, 0.5, 1.0, 10, 1000);

  assert(result.states.length >= 1, 'Should have at least 1 state');

  const groundState = result.states[0];
  const psi = groundState.psi;

  const leftBoundary = Math.abs(psi[0]);
  const rightBoundary = Math.abs(psi[psi.length - 1]);

  const maxPsi = Math.max(...psi.map(Math.abs));

  assert(leftBoundary / maxPsi < 0.01, 'Left boundary should be small');
  assert(rightBoundary / maxPsi < 0.01, 'Right boundary should be small');

  console.log(`  ✓ Left boundary decay: ${(leftBoundary / maxPsi * 100).toFixed(3)}%`);
  console.log(`  ✓ Right boundary decay: ${(rightBoundary / maxPsi * 100).toFixed(3)}%`);
  console.log(`  ✓ Test PASSED`);
} catch (error) {
  console.error(`  ✗ Test FAILED: ${error}`);
  process.exit(1);
}

// Test 16: Ground State Has No Nodes
console.log('\n[Test 16] Ground State - No Nodes');
console.log('-'.repeat(70));
try {
  const result = solveDoubleWell(1.0, 5.0, 0.5, 1.0, 10, 1000);

  assert(result.states.length >= 1, 'Should have at least 1 state');

  const groundState = result.states[0];
  const psi = groundState.psi;

  let signChanges = 0;
  for (let i = 1; i < psi.length; i++) {
    if (psi[i] * psi[i - 1] < 0) {
      signChanges++;
    }
  }

  assert(signChanges === 0, `Ground state should have 0 nodes (found ${signChanges})`);

  console.log(`  ✓ Sign changes (nodes): ${signChanges}`);
  console.log(`  ✓ Test PASSED`);
} catch (error) {
  console.error(`  ✗ Test FAILED: ${error}`);
  process.exit(1);
}

// Test 17: First Excited State Has One Node
console.log('\n[Test 17] First Excited State - One Node');
console.log('-'.repeat(70));
try {
  const result = solveDoubleWell(1.0, 5.0, 0.5, 1.0, 10, 1000);

  assert(result.states.length >= 2, 'Should have at least 2 states');

  const firstExcited = result.states[1];
  const psi = firstExcited.psi;

  let signChanges = 0;
  for (let i = 1; i < psi.length; i++) {
    if (psi[i] * psi[i - 1] < 0) {
      signChanges++;
    }
  }

  assert(signChanges === 1, `First excited state should have 1 node (found ${signChanges})`);

  console.log(`  ✓ Sign changes (nodes): ${signChanges}`);
  console.log(`  ✓ Test PASSED`);
} catch (error) {
  console.error(`  ✗ Test FAILED: ${error}`);
  process.exit(1);
}

// Test 18: Monotonically Increasing Energies
console.log('\n[Test 18] Monotonically Increasing Energies');
console.log('-'.repeat(70));
try {
  const result = solveDoubleWell(1.0, 5.0, 0.5, 1.0, 10, 1000);

  assert(result.states.length >= 2, 'Should have at least 2 states');

  for (let i = 1; i < result.states.length; i++) {
    // Allow for numerical precision (energies should be >= within tolerance)
    const tolerance = 1e-10; // 0.0001 meV
    assert(result.states[i].energy >= result.states[i - 1].energy - tolerance,
      `Energy ${i} should be >= energy ${i - 1} (${result.states[i].energy} >= ${result.states[i - 1].energy})`);
  }

  console.log(`  ✓ All ${result.states.length} energies are monotonically non-decreasing`);
  console.log(`  ✓ Test PASSED`);
} catch (error) {
  console.error(`  ✗ Test FAILED: ${error}`);
  process.exit(1);
}

// Test 19: Alternating Parity
console.log('\n[Test 19] Alternating Parity');
console.log('-'.repeat(70));
try {
  const result = solveDoubleWell(1.0, 5.0, 0.5, 1.0, 10, 1000);

  assert(result.states.length >= 4, 'Should have at least 4 states');

  const parities = result.states.slice(0, 4).map(s => determineParity(result.x, s.psi));

  assert(parities[0] === 'even', `State 0 should be even (got ${parities[0]})`);
  assert(parities[1] === 'odd', `State 1 should be odd (got ${parities[1]})`);
  assert(parities[2] === 'even', `State 2 should be even (got ${parities[2]})`);
  assert(parities[3] === 'odd', `State 3 should be odd (got ${parities[3]})`);

  console.log(`  ✓ Parities: ${parities.join(', ')}`);
  console.log(`  ✓ Test PASSED`);
} catch (error) {
  console.error(`  ✗ Test FAILED: ${error}`);
  process.exit(1);
}

// Test 20: Parameter Space Sweep
console.log('\n[Test 20] Comprehensive Parameter Space Sweep');
console.log('-'.repeat(70));
try {
  const wellWidths = [0.5, 1.0, 2.0, 3.0];
  const wellDepths = [1.0, 5.0, 10.0, 15.0];
  const barrierWidths = [0.1, 0.5, 1.0, 2.0];

  let sweepTotal = 0;
  let sweepPassed = 0;

  for (const wellWidth of wellWidths) {
    for (const wellDepth of wellDepths) {
      for (const barrierWidth of barrierWidths) {
        sweepTotal++;

        try {
          const result = solveDoubleWell(wellWidth, wellDepth, barrierWidth, 1.0, 10, 800);

          if (result.states.length >= 1) {
            const allFinite = result.states.every(s =>
              isFinite(s.energy) &&
              s.psi.every(p => isFinite(p))
            );

            const allBound = result.states.every(s =>
              s.energy < 0 && s.energy > -wellDepth
            );

            if (allFinite && allBound) {
              sweepPassed++;
            }
          }
        } catch (_error) {
          // Skip failed configurations
        }
      }
    }
  }

  const successRate = sweepPassed / sweepTotal;
  assert(successRate > 0.95, `Parameter sweep success rate should be > 95% (got ${(successRate * 100).toFixed(1)}%)`);

  console.log(`  ✓ Parameter sweep: ${sweepPassed}/${sweepTotal} passed (${(successRate * 100).toFixed(1)}%)`);
  console.log(`  ✓ Test PASSED`);
} catch (error) {
  console.error(`  ✗ Test FAILED: ${error}`);
  process.exit(1);
}

// Final Summary
console.log('\n' + '='.repeat(70));
console.log(' TEST SUMMARY');
console.log('='.repeat(70));
console.log(`Total assertions: ${totalTests}`);
console.log(`Passed: ${passedTests}`);
console.log(`Failed: ${failedTests}`);
console.log('='.repeat(70));

if (failedTests === 0) {
  console.log('\n✓ ALL TESTS PASSED!\n');
  process.exit(0);
} else {
  console.log('\n✗ SOME TESTS FAILED\n');
  process.exit(1);
}
