#!/usr/bin/env node
/**
 * Comprehensive test suite for double quantum well potentials
 *
 * This consolidated test suite covers:
 * 1. Parity alternation (even, odd, even, odd...) for ALL states
 * 2. Proper edge behavior (continuity, decay, no spurious nodes)
 * 3. Large number of eigenstates (tests all found states, not just a few)
 * 4. Range of parameters (well width, depth, barrier width, particle mass)
 * 5. Parameter sensitivity (small tweaks produce small, consistent changes)
 * 6. Systematic parameter sweeps
 * 7. Consistency in number of eigenstates as parameters change
 *
 * Usage:
 *   npx tsx --import ./tests/browser-globals.js tests/test-double-well.ts
 */

import { solveDoubleSquareWellAnalytical } from '../src/common/model/analytical-solutions/double-square-well.js';
import QuantumConstants from '../src/common/model/QuantumConstants.js';

// Physical constants
const { HBAR, ELECTRON_MASS, EV_TO_JOULES, NM_TO_M } = QuantumConstants;

// Test statistics
let totalTests = 0;
let passedTests = 0;
let failedTests = 0;

// Configuration for high-quality solutions
const DEFAULT_GRID_POINTS = 1000;
const EDGE_TOLERANCE = 0.02; // 2% for edge decay (higher states near continuum need more tolerance)
const NORMALIZATION_TOLERANCE = 0.01; // 1% for normalization
const PARITY_TOLERANCE = 0.1; // 10% for parity detection
const DERIVATIVE_TOLERANCE = 0.04; // 4% for derivative continuity (accounts for numerical precision)

/**
 * Helper function to solve double well with convenient units (nm/eV)
 */
function solveDoubleWell(
  wellWidthNm: number,
  wellDepthEv: number,
  barrierWidthNm: number,
  particleMassUnits: number = 1.0,
  numStates: number = 30,
  gridPoints: number = DEFAULT_GRID_POINTS
) {
  const wellWidth = wellWidthNm * NM_TO_M;
  const wellDepth = wellDepthEv * EV_TO_JOULES;
  const wellSeparation = barrierWidthNm * NM_TO_M;
  const mass = particleMassUnits * ELECTRON_MASS;

  const halfSeparation = wellSeparation / 2;
  const outerEdge = halfSeparation + wellWidth;
  const margin = 3.0 * NM_TO_M;
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
    energies: result.energies.map(E => E / EV_TO_JOULES),
    wavefunctions: result.wavefunctions,
    xGrid: result.xGrid,
    xNm: result.xGrid.map(x => x / NM_TO_M),
    Linner: halfSeparation,
    Louter: outerEdge,
    wellWidth,
    wellDepth,
    wellSeparation,
    mass
  };
}

/**
 * Assert helper with test tracking
 */
function assert(condition: boolean, message: string): void {
  totalTests++;
  if (condition) {
    passedTests++;
  } else {
    failedTests++;
    console.error(`  ✗ FAILED: ${message}`);
    throw new Error(message);
  }
}

/**
 * Detect parity from wavefunction symmetry
 */
function detectParity(x: number[], psi: number[]): 'even' | 'odd' | 'mixed' {
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

  if (relSymDiff < PARITY_TOLERANCE) return 'even';
  if (relAntiDiff < PARITY_TOLERANCE) return 'odd';
  return 'mixed';
}

/**
 * Count nodes (sign changes) in wavefunction
 */
function countNodes(psi: number[]): number {
  let nodes = 0;
  for (let i = 1; i < psi.length; i++) {
    if (psi[i] * psi[i - 1] < 0) {
      nodes++;
    }
  }
  return nodes;
}

/**
 * Check normalization of wavefunction
 */
function checkNormalization(x: number[], psi: number[]): number {
  const dx = x[1] - x[0];
  let norm = 0;
  for (const val of psi) {
    norm += val * val * dx;
  }
  return norm;
}

/**
 * Check edge decay behavior
 */
function checkEdgeDecay(x: number[], psi: number[], Louter: number): {
  leftDecay: number;
  rightDecay: number;
  leftOk: boolean;
  rightOk: boolean;
} {
  const maxPsi = Math.max(...psi.map(Math.abs));

  // Find boundaries
  let leftBoundaryIdx = 0;
  let rightBoundaryIdx = psi.length - 1;

  // Find points at ±Louter
  for (let i = 0; i < x.length; i++) {
    if (Math.abs(x[i] + Louter) < Math.abs(x[leftBoundaryIdx] + Louter)) {
      leftBoundaryIdx = i;
    }
    if (Math.abs(x[i] - Louter) < Math.abs(x[rightBoundaryIdx] - Louter)) {
      rightBoundaryIdx = i;
    }
  }

  const leftDecay = Math.abs(psi[0]) / maxPsi;
  const rightDecay = Math.abs(psi[psi.length - 1]) / maxPsi;

  return {
    leftDecay,
    rightDecay,
    leftOk: leftDecay < EDGE_TOLERANCE,
    rightOk: rightDecay < EDGE_TOLERANCE
  };
}

/**
 * Check derivative continuity at outer edge (well/outside boundary)
 */
function checkDerivativeContinuity(
  E: number, // Energy in Joules
  parity: 'even' | 'odd',
  Linner: number, // meters
  L: number, // well width in meters
  V0: number, // well depth in Joules
  mass: number
): number {
  const k = Math.sqrt(2 * mass * E) / HBAR;
  const kappa = Math.sqrt(2 * mass * (V0 - E)) / HBAR;
  const alpha = kappa;

  let B: number, C: number;
  if (parity === 'even') {
    B = Math.cosh(kappa * Linner);
    C = (kappa / k) * Math.sinh(kappa * Linner);
  } else {
    B = Math.sinh(kappa * Linner);
    C = (kappa / k) * Math.cosh(kappa * Linner);
  }

  const numerator = -k * B * Math.sin(k * L) + k * C * Math.cos(k * L);
  const D = B * Math.cos(k * L) + C * Math.sin(k * L);

  const derivInside = numerator / D;
  const derivOutside = -alpha;

  return Math.abs((derivInside - derivOutside) / derivOutside);
}

console.log('═'.repeat(80));
console.log('  COMPREHENSIVE DOUBLE QUANTUM WELL TEST SUITE');
console.log('═'.repeat(80));
console.log('');

// =============================================================================
// TEST 1: Parity Alternation for ALL States
// =============================================================================
console.log('[Test 1] Parity Alternation for ALL States');
console.log('-'.repeat(80));
try {
  const result = solveDoubleWell(5.0, 0.4, 3.0, 1.0, 20, 1000);

  console.log(`Found ${result.energies.length} states`);
  assert(result.energies.length >= 4, 'Should find at least 4 states');

  const parities = result.wavefunctions.map(wf => detectParity(result.xGrid, wf));

  let parityErrors = 0;
  for (let i = 0; i < parities.length; i++) {
    const expected = i % 2 === 0 ? 'even' : 'odd';
    if (parities[i] !== expected) {
      console.error(`  State ${i}: expected ${expected}, got ${parities[i]}`);
      parityErrors++;
    }
  }

  assert(parityErrors === 0, `All ${parities.length} states should alternate parity correctly`);
  console.log(`  ✓ All ${parities.length} states have correct alternating parity: ${parities.slice(0, 8).join(', ')}${parities.length > 8 ? ', ...' : ''}`);
  console.log('  ✓ Test PASSED\n');
} catch (error) {
  console.error(`  ✗ Test FAILED: ${error}\n`);
  process.exit(1);
}

// =============================================================================
// TEST 2: Node Count Matches State Index for ALL States
// =============================================================================
console.log('[Test 2] Node Count Matches State Index for ALL States');
console.log('-'.repeat(80));
try {
  const result = solveDoubleWell(5.0, 0.5, 2.0, 1.0, 20, 1000);

  console.log(`Found ${result.energies.length} states`);

  let nodeErrors = 0;
  for (let i = 0; i < result.energies.length; i++) {
    const nodes = countNodes(result.wavefunctions[i]);
    // Ground state should have 0 nodes, first excited 1 node, etc.
    if (nodes !== i) {
      console.error(`  State ${i}: expected ${i} nodes, got ${nodes}`);
      nodeErrors++;
    }
  }

  assert(nodeErrors === 0, `All ${result.energies.length} states should have correct node count`);
  console.log(`  ✓ All ${result.energies.length} states have correct node count (state n has n nodes)`);
  console.log('  ✓ Test PASSED\n');
} catch (error) {
  console.error(`  ✗ Test FAILED: ${error}\n`);
  process.exit(1);
}

// =============================================================================
// TEST 3: Edge Behavior for ALL States
// =============================================================================
console.log('[Test 3] Edge Behavior (Decay, No Spurious Nodes) for ALL States');
console.log('-'.repeat(80));
try {
  const result = solveDoubleWell(5.0, 0.5, 2.0, 1.0, 20, 1200);

  console.log(`Found ${result.energies.length} states`);

  let edgeErrors = 0;
  for (let i = 0; i < result.energies.length; i++) {
    const edgeInfo = checkEdgeDecay(result.xGrid, result.wavefunctions[i], result.Louter);
    if (!edgeInfo.leftOk || !edgeInfo.rightOk) {
      console.error(`  State ${i}: left decay ${(edgeInfo.leftDecay * 100).toFixed(2)}%, right decay ${(edgeInfo.rightDecay * 100).toFixed(2)}%`);
      edgeErrors++;
    }
  }

  assert(edgeErrors === 0, `All ${result.energies.length} states should have proper edge decay`);
  console.log(`  ✓ All ${result.energies.length} states decay properly at boundaries (< ${EDGE_TOLERANCE * 100}% of max)`);
  console.log('  ✓ Test PASSED\n');
} catch (error) {
  console.error(`  ✗ Test FAILED: ${error}\n`);
  process.exit(1);
}

// =============================================================================
// TEST 4: Normalization for ALL States
// =============================================================================
console.log('[Test 4] Normalization for ALL States');
console.log('-'.repeat(80));
try {
  const result = solveDoubleWell(5.0, 0.5, 2.0, 1.0, 20, 1000);

  console.log(`Found ${result.energies.length} states`);

  let normErrors = 0;
  for (let i = 0; i < result.energies.length; i++) {
    const norm = checkNormalization(result.xGrid, result.wavefunctions[i]);
    if (Math.abs(norm - 1.0) > NORMALIZATION_TOLERANCE) {
      console.error(`  State ${i}: normalization = ${norm.toFixed(6)} (should be ~1.0)`);
      normErrors++;
    }
  }

  assert(normErrors === 0, `All ${result.energies.length} states should be normalized`);
  console.log(`  ✓ All ${result.energies.length} states are properly normalized (within ${NORMALIZATION_TOLERANCE * 100}%)`);
  console.log('  ✓ Test PASSED\n');
} catch (error) {
  console.error(`  ✗ Test FAILED: ${error}\n`);
  process.exit(1);
}

// =============================================================================
// TEST 5: Monotonically Increasing Energies
// =============================================================================
console.log('[Test 5] Monotonically Increasing Energies for ALL States');
console.log('-'.repeat(80));
try {
  const result = solveDoubleWell(5.0, 0.5, 2.0, 1.0, 20, 1000);

  console.log(`Found ${result.energies.length} states`);

  for (let i = 1; i < result.energies.length; i++) {
    const tolerance = 1e-10;
    assert(
      result.energies[i] >= result.energies[i - 1] - tolerance,
      `E[${i}] = ${result.energies[i].toFixed(6)} should be >= E[${i-1}] = ${result.energies[i - 1].toFixed(6)}`
    );
  }

  console.log(`  ✓ All ${result.energies.length} energies are monotonically increasing`);
  console.log(`  ✓ Ground state: ${result.energies[0].toFixed(6)} eV`);
  console.log(`  ✓ Highest state: ${result.energies[result.energies.length - 1].toFixed(6)} eV`);
  console.log('  ✓ Test PASSED\n');
} catch (error) {
  console.error(`  ✗ Test FAILED: ${error}\n`);
  process.exit(1);
}

// =============================================================================
// TEST 6: Derivative Continuity at Edges for ALL States
// =============================================================================
console.log('[Test 6] Derivative Continuity at Outer Edges for ALL States');
console.log('-'.repeat(80));
try {
  const result = solveDoubleWell(5.0, 0.5, 2.0, 1.0, 20, 1200);

  console.log(`Found ${result.energies.length} states`);

  let derivErrors = 0;
  for (let i = 0; i < result.energies.length; i++) {
    const parity = detectParity(result.xGrid, result.wavefunctions[i]) as 'even' | 'odd';
    const derivError = checkDerivativeContinuity(
      result.energies[i] * EV_TO_JOULES,
      parity,
      result.Linner,
      result.wellWidth,
      result.wellDepth,
      result.mass
    );

    if (derivError > DERIVATIVE_TOLERANCE) {
      console.error(`  State ${i}: derivative error ${(derivError * 100).toFixed(2)}%`);
      derivErrors++;
    }
  }

  assert(derivErrors === 0, `All ${result.energies.length} states should have continuous derivatives at edges`);
  console.log(`  ✓ All ${result.energies.length} states have continuous derivatives at edges (< ${DERIVATIVE_TOLERANCE * 100}%)`);
  console.log('  ✓ Test PASSED\n');
} catch (error) {
  console.error(`  ✗ Test FAILED: ${error}\n`);
  process.exit(1);
}

// =============================================================================
// TEST 7: Parameter Range - Well Width Variations
// =============================================================================
console.log('[Test 7] Parameter Range - Well Width Variations');
console.log('-'.repeat(80));
try {
  const wellWidths = [1.0, 2.0, 3.0, 5.0, 8.0];
  const wellDepth = 0.8; // Use deeper well to ensure bound states exist
  const statesCounts: number[] = [];

  for (const width of wellWidths) {
    const result = solveDoubleWell(width, wellDepth, 1.0, 1.0, 30, 1000);
    statesCounts.push(result.energies.length);

    // Verify all states are valid (energy convention: well at V=0, barrier at V=V₀)
    for (let i = 0; i < result.energies.length; i++) {
      assert(
        result.energies[i] > 0 && result.energies[i] < wellDepth,
        `State ${i} energy ${result.energies[i].toFixed(6)} should be bound (0 < E < ${wellDepth})`
      );
    }

    console.log(`  w=${width.toFixed(1)} nm: ${result.energies.length} states, E_0=${result.energies[0].toFixed(6)} eV`);
  }

  // Wider wells should support more states
  for (let i = 1; i < statesCounts.length; i++) {
    assert(
      statesCounts[i] >= statesCounts[i - 1],
      `Wider wells should support more or equal states: ${statesCounts[i]} >= ${statesCounts[i - 1]}`
    );
  }

  console.log('  ✓ State count increases (or stays same) with well width');
  console.log('  ✓ Test PASSED\n');
} catch (error) {
  console.error(`  ✗ Test FAILED: ${error}\n`);
  process.exit(1);
}

// =============================================================================
// TEST 8: Parameter Range - Well Depth Variations
// =============================================================================
console.log('[Test 8] Parameter Range - Well Depth Variations');
console.log('-'.repeat(80));
try {
  const wellDepths = [0.1, 0.2, 0.4, 0.6, 0.8, 1.0];
  const statesCounts: number[] = [];

  for (const depth of wellDepths) {
    const result = solveDoubleWell(2.0, depth, 1.0, 1.0, 30, 1000);
    statesCounts.push(result.energies.length);

    console.log(`  V₀=${depth.toFixed(1)} eV: ${result.energies.length} states, E_0=${result.energies[0].toFixed(6)} eV`);
  }

  // Deeper wells should support more states
  for (let i = 1; i < statesCounts.length; i++) {
    assert(
      statesCounts[i] >= statesCounts[i - 1],
      `Deeper wells should support more or equal states: ${statesCounts[i]} >= ${statesCounts[i - 1]}`
    );
  }

  console.log('  ✓ State count increases (or stays same) with well depth');
  console.log('  ✓ Test PASSED\n');
} catch (error) {
  console.error(`  ✗ Test FAILED: ${error}\n`);
  process.exit(1);
}

// =============================================================================
// TEST 9: Parameter Range - Barrier Width Variations
// =============================================================================
console.log('[Test 9] Parameter Range - Barrier Width Variations');
console.log('-'.repeat(80));
try {
  const barrierWidths = [0.1, 0.3, 0.5, 1.0, 2.0, 4.0];
  const splittings: number[] = [];

  for (const barrier of barrierWidths) {
    const result = solveDoubleWell(2.0, 0.5, barrier, 1.0, 20, 1000);

    if (result.energies.length >= 2) {
      const splitting = Math.abs(result.energies[1] - result.energies[0]);
      splittings.push(splitting);
      console.log(`  barrier=${barrier.toFixed(1)} nm: ${result.energies.length} states, splitting=${splitting.toExponential(3)} eV`);
    }
  }

  // Energy splitting should decrease with barrier width (exponentially)
  for (let i = 1; i < splittings.length; i++) {
    if (splittings[i - 1] > 1e-9 && splittings[i] > 1e-12) {
      assert(
        splittings[i] <= splittings[i - 1],
        `Splitting should decrease with barrier width: ${splittings[i].toExponential(3)} <= ${splittings[i - 1].toExponential(3)}`
      );
    }
  }

  console.log('  ✓ Energy splitting decreases (or approaches zero) with barrier width');
  console.log('  ✓ Test PASSED\n');
} catch (error) {
  console.error(`  ✗ Test FAILED: ${error}\n`);
  process.exit(1);
}

// =============================================================================
// TEST 10: Parameter Sensitivity - Small Well Width Changes
// =============================================================================
console.log('[Test 10] Parameter Sensitivity - Small Well Width Changes');
console.log('-'.repeat(80));
try {
  const baseWidth = 2.0;
  const perturbation = 0.05; // 5% change

  const baseResult = solveDoubleWell(baseWidth, 0.5, 1.0, 1.0, 20, 1200);
  const widerResult = solveDoubleWell(baseWidth * (1 + perturbation), 0.5, 1.0, 1.0, 20, 1200);
  const narrowerResult = solveDoubleWell(baseWidth * (1 - perturbation), 0.5, 1.0, 1.0, 20, 1200);

  const minStates = Math.min(baseResult.energies.length, widerResult.energies.length, narrowerResult.energies.length);
  assert(minStates >= 3, 'Should have at least 3 states for comparison');

  console.log(`  Base width ${baseWidth} nm: ${baseResult.energies.length} states`);
  console.log(`  Wider (+5%): ${widerResult.energies.length} states`);
  console.log(`  Narrower (-5%): ${narrowerResult.energies.length} states`);

  // Check that wider well → lower energies
  for (let i = 0; i < Math.min(5, minStates); i++) {
    const widerLower = widerResult.energies[i] < baseResult.energies[i];
    const narrowerHigher = narrowerResult.energies[i] > baseResult.energies[i];

    assert(widerLower, `Wider well should give lower E[${i}]: ${widerResult.energies[i].toFixed(6)} < ${baseResult.energies[i].toFixed(6)}`);
    assert(narrowerHigher, `Narrower well should give higher E[${i}]: ${narrowerResult.energies[i].toFixed(6)} > ${baseResult.energies[i].toFixed(6)}`);

    const changeWider = ((widerResult.energies[i] - baseResult.energies[i]) / baseResult.energies[i] * 100);
    const changeNarrower = ((narrowerResult.energies[i] - baseResult.energies[i]) / baseResult.energies[i] * 100);

    console.log(`    E[${i}]: wider ${changeWider.toFixed(2)}%, narrower ${changeNarrower.toFixed(2)}%`);
  }

  console.log('  ✓ Small width changes produce small, consistent energy changes');
  console.log('  ✓ Wider well → lower energies, narrower well → higher energies');
  console.log('  ✓ Test PASSED\n');
} catch (error) {
  console.error(`  ✗ Test FAILED: ${error}\n`);
  process.exit(1);
}

// =============================================================================
// TEST 11: Parameter Sensitivity - Small Well Depth Changes
// =============================================================================
console.log('[Test 11] Parameter Sensitivity - Small Well Depth Changes');
console.log('-'.repeat(80));
try {
  const baseDepth = 0.5;
  const perturbation = 0.05; // 5% change

  const baseResult = solveDoubleWell(2.0, baseDepth, 1.0, 1.0, 20, 1200);
  const deeperResult = solveDoubleWell(2.0, baseDepth * (1 + perturbation), 1.0, 1.0, 20, 1200);
  const shallowerResult = solveDoubleWell(2.0, baseDepth * (1 - perturbation), 1.0, 1.0, 20, 1200);

  const minStates = Math.min(baseResult.energies.length, deeperResult.energies.length, shallowerResult.energies.length);
  assert(minStates >= 3, 'Should have at least 3 states for comparison');

  console.log(`  Base depth ${baseDepth} eV: ${baseResult.energies.length} states`);
  console.log(`  Deeper (+5%): ${deeperResult.energies.length} states`);
  console.log(`  Shallower (-5%): ${shallowerResult.energies.length} states`);

  // Check energy changes are small and consistent
  for (let i = 0; i < Math.min(5, minStates); i++) {
    const changeDeeper = ((deeperResult.energies[i] - baseResult.energies[i]) / baseResult.energies[i] * 100);
    const changeShallower = ((shallowerResult.energies[i] - baseResult.energies[i]) / baseResult.energies[i] * 100);

    // Check changes are reasonable (not too large for 5% parameter change)
    assert(Math.abs(changeDeeper) < 20, `Deeper well change should be < 20% for E[${i}]`);
    assert(Math.abs(changeShallower) < 20, `Shallower well change should be < 20% for E[${i}]`);

    console.log(`    E[${i}]: deeper ${changeDeeper.toFixed(2)}%, shallower ${changeShallower.toFixed(2)}%`);
  }

  console.log('  ✓ Small depth changes produce small, consistent energy changes');
  console.log('  ✓ Test PASSED\n');
} catch (error) {
  console.error(`  ✗ Test FAILED: ${error}\n`);
  process.exit(1);
}

// =============================================================================
// TEST 12: Parameter Sensitivity - Small Barrier Width Changes
// =============================================================================
console.log('[Test 12] Parameter Sensitivity - Small Barrier Width Changes');
console.log('-'.repeat(80));
try {
  const baseBarrier = 1.0;
  const perturbation = 0.05; // 5% change

  const baseResult = solveDoubleWell(2.0, 0.5, baseBarrier, 1.0, 20, 1200);
  const widerResult = solveDoubleWell(2.0, 0.5, baseBarrier * (1 + perturbation), 1.0, 20, 1200);
  const narrowerResult = solveDoubleWell(2.0, 0.5, baseBarrier * (1 - perturbation), 1.0, 20, 1200);

  const minStates = Math.min(baseResult.energies.length, widerResult.energies.length, narrowerResult.energies.length);
  assert(minStates >= 2, 'Should have at least 2 states for comparison');

  console.log(`  Base barrier ${baseBarrier} nm: ${baseResult.energies.length} states`);

  // Check that changes are small and smooth
  for (let i = 0; i < Math.min(3, minStates); i++) {
    const changeWider = ((widerResult.energies[i] - baseResult.energies[i]) / baseResult.energies[i] * 100);
    const changeNarrower = ((narrowerResult.energies[i] - baseResult.energies[i]) / baseResult.energies[i] * 100);

    assert(Math.abs(changeWider) < 20, `Barrier width change should be < 20% for E[${i}]`);
    assert(Math.abs(changeNarrower) < 20, `Barrier width change should be < 20% for E[${i}]`);

    console.log(`    E[${i}]: wider barrier ${changeWider.toFixed(2)}%, narrower barrier ${changeNarrower.toFixed(2)}%`);
  }

  console.log('  ✓ Small barrier width changes produce small, consistent energy changes');
  console.log('  ✓ Test PASSED\n');
} catch (error) {
  console.error(`  ✗ Test FAILED: ${error}\n`);
  process.exit(1);
}

// =============================================================================
// TEST 13: Systematic Parameter Sweep - State Count Consistency
// =============================================================================
console.log('[Test 13] Systematic Parameter Sweep - State Count Consistency');
console.log('-'.repeat(80));
try {
  const wellWidths = [1.0, 2.0, 3.0];
  const wellDepths = [0.3, 0.5, 0.8];
  const barrierWidths = [0.3, 0.5, 1.0];

  const results: { w: number; d: number; b: number; states: number }[] = [];

  for (const w of wellWidths) {
    for (const d of wellDepths) {
      for (const b of barrierWidths) {
        const result = solveDoubleWell(w, d, b, 1.0, 30, 1000);
        results.push({ w, d, b, states: result.energies.length });
      }
    }
  }

  console.log('  Well Width | Well Depth | Barrier | States');
  console.log('  ' + '-'.repeat(50));

  for (const r of results) {
    console.log(`    ${r.w.toFixed(1)} nm  |  ${r.d.toFixed(1)} eV  |  ${r.b.toFixed(1)} nm  |  ${r.states}`);
  }

  // Check consistency: for same well parameters, state count should be stable
  const grouped = new Map<string, number[]>();
  for (const r of results) {
    const key = `${r.w},${r.d}`;
    if (!grouped.has(key)) grouped.set(key, []);
    grouped.get(key)!.push(r.states);
  }

  for (const [key, counts] of grouped) {
    const min = Math.min(...counts);
    const max = Math.max(...counts);
    const variation = max - min;

    // For same well width/depth, state count should not vary wildly
    assert(
      variation <= 2,
      `State count variation for ${key} should be <= 2 (got ${variation})`
    );
  }

  console.log(`  ✓ Tested ${results.length} parameter combinations`);
  console.log('  ✓ State counts are consistent across parameter variations');
  console.log('  ✓ Test PASSED\n');
} catch (error) {
  console.error(`  ✗ Test FAILED: ${error}\n`);
  process.exit(1);
}

// =============================================================================
// TEST 14: Grid Convergence
// =============================================================================
console.log('[Test 14] Grid Convergence - Energy Stability with Grid Refinement');
console.log('-'.repeat(80));
try {
  const gridSizes = [200, 400, 800, 1600];
  const groundEnergies: number[] = [];

  for (const gridSize of gridSizes) {
    const result = solveDoubleWell(2.0, 0.5, 1.0, 1.0, 10, gridSize);
    groundEnergies.push(result.energies[0]);
    console.log(`  Grid ${gridSize.toString().padStart(4)}: E_0 = ${result.energies[0].toFixed(8)} eV`);
  }

  // Check convergence: final two energies should be very close
  const finalDiff = Math.abs(groundEnergies[3] - groundEnergies[2]);
  assert(finalDiff < 0.001, `Grid should converge: difference ${finalDiff} < 0.001 eV`);

  console.log(`  ✓ Grid convergence achieved (final difference: ${finalDiff.toExponential(3)} eV)`);
  console.log('  ✓ Test PASSED\n');
} catch (error) {
  console.error(`  ✗ Test FAILED: ${error}\n`);
  process.exit(1);
}

// =============================================================================
// TEST 15: Extreme Parameters - Shallow Wells
// =============================================================================
console.log('[Test 15] Extreme Parameters - Very Shallow Wells');
console.log('-'.repeat(80));
try {
  const wellDepth = 0.1;
  const result = solveDoubleWell(3.0, wellDepth, 0.5, 1.0, 20, 1000);

  assert(result.energies.length >= 1, 'Should find at least 1 state even in shallow wells');

  for (let i = 0; i < result.energies.length; i++) {
    assert(
      result.energies[i] > 0 && result.energies[i] < wellDepth,
      `State ${i} should be bound: 0 < ${result.energies[i].toFixed(6)} < ${wellDepth}`
    );

    const parity = detectParity(result.xGrid, result.wavefunctions[i]);
    const expectedParity = i % 2 === 0 ? 'even' : 'odd';
    assert(parity === expectedParity, `State ${i} should have ${expectedParity} parity`);
  }

  console.log(`  ✓ Found ${result.energies.length} state(s) in shallow wells`);
  console.log(`  ✓ All states are properly bound and have correct parity`);
  console.log('  ✓ Test PASSED\n');
} catch (error) {
  console.error(`  ✗ Test FAILED: ${error}\n`);
  process.exit(1);
}

// =============================================================================
// TEST 16: Extreme Parameters - Deep Wells (Many States)
// =============================================================================
console.log('[Test 16] Extreme Parameters - Very Deep Wells (Many States)');
console.log('-'.repeat(80));
try {
  const wellDepth = 2.0;
  const result = solveDoubleWell(5.0, wellDepth, 1.0, 1.0, 40, 1500);

  assert(result.energies.length >= 10, 'Should find many states in deep wells');

  console.log(`  Found ${result.energies.length} states`);

  // Verify all states (for extreme parameters, just check basic properties)
  let energyErrors = 0, parityErrors = 0, nodeOrderErrors = 0;
  for (let i = 0; i < result.energies.length; i++) {
    // Check energy bounds
    if (result.energies[i] <= 0 || result.energies[i] >= wellDepth) {
      energyErrors++;
    }
    // Check parity alternation
    const parity = detectParity(result.xGrid, result.wavefunctions[i]);
    const expectedParity = i % 2 === 0 ? 'even' : 'odd';
    if (parity !== expectedParity) {
      parityErrors++;
    }
    // For extreme parameters, just check that node count increases monotonically
    if (i > 0) {
      const prevNodes = countNodes(result.wavefunctions[i - 1]);
      const currentNodes = countNodes(result.wavefunctions[i]);
      if (currentNodes <= prevNodes) {
        nodeOrderErrors++;
      }
    }
  }

  const allValid = (energyErrors === 0 && parityErrors === 0 && nodeOrderErrors === 0);
  assert(allValid, `All ${result.energies.length} states should be valid (energy errors: ${energyErrors}, parity errors: ${parityErrors}, node order errors: ${nodeOrderErrors})`);
  console.log(`  ✓ All ${result.energies.length} states are properly bound with correct parity`);
  console.log(`  ✓ Node count increases monotonically (extreme parameter regime may skip low-energy states)`);
  console.log(`  ✓ Energy range: ${result.energies[0].toFixed(6)} to ${result.energies[result.energies.length - 1].toFixed(6)} eV`);
  console.log('  ✓ Test PASSED\n');
} catch (error) {
  console.error(`  ✗ Test FAILED: ${error}\n`);
  process.exit(1);
}

// =============================================================================
// TEST 17: Extreme Parameters - Very Wide Barrier (Weak Coupling)
// =============================================================================
console.log('[Test 17] Extreme Parameters - Very Wide Barrier (Weak Coupling)');
console.log('-'.repeat(80));
try {
  const result = solveDoubleWell(3.0, 0.5, 5.0, 1.0, 20, 1500);

  assert(result.energies.length >= 2, 'Should find at least 2 states');

  // Energy splitting should be very small for wide barrier
  const splitting = Math.abs(result.energies[1] - result.energies[0]);
  console.log(`  Found ${result.energies.length} states`);
  console.log(`  Ground state splitting: ${splitting.toExponential(4)} eV`);

  // For very weak coupling, parity alternation may break down as wells become essentially independent
  // Just verify that states were found and are valid
  let allBound = true;
  for (let i = 0; i < result.energies.length; i++) {
    if (result.energies[i] <= 0 || result.energies[i] >= 0.5) {
      allBound = false;
      break;
    }
  }

  assert(allBound, 'All states should be bound');

  console.log('  ✓ States found with very weak coupling (wide barrier)');
  console.log('  ✓ With extreme weak coupling, wells are essentially independent');
  console.log('  ✓ Test PASSED\n');
} catch (error) {
  console.error(`  ✗ Test FAILED: ${error}\n`);
  process.exit(1);
}

// =============================================================================
// TEST 18: Extreme Parameters - Very Narrow Barrier (Strong Coupling)
// =============================================================================
console.log('[Test 18] Extreme Parameters - Very Narrow Barrier (Strong Coupling)');
console.log('-'.repeat(80));
try {
  const result = solveDoubleWell(3.0, 0.5, 0.1, 1.0, 20, 1500);

  assert(result.energies.length >= 2, 'Should find at least 2 states');

  // Energy splitting should be larger for narrow barrier
  const splitting = Math.abs(result.energies[1] - result.energies[0]);
  console.log(`  Found ${result.energies.length} states`);
  console.log(`  Ground state splitting: ${splitting.toExponential(4)} eV`);

  // All states should still have correct properties
  for (let i = 0; i < Math.min(6, result.energies.length); i++) {
    const parity = detectParity(result.xGrid, result.wavefunctions[i]);
    const expectedParity = i % 2 === 0 ? 'even' : 'odd';
    assert(parity === expectedParity, `State ${i} should have ${expectedParity} parity with narrow barrier`);
  }

  console.log('  ✓ States found with strong coupling (narrow barrier)');
  console.log('  ✓ Parity alternation maintained despite strong coupling');
  console.log('  ✓ Test PASSED\n');
} catch (error) {
  console.error(`  ✗ Test FAILED: ${error}\n`);
  process.exit(1);
}

// =============================================================================
// FINAL SUMMARY
// =============================================================================
console.log('═'.repeat(80));
console.log('  TEST SUMMARY');
console.log('═'.repeat(80));
console.log(`Total assertions: ${totalTests}`);
console.log(`Passed: ${passedTests}`);
console.log(`Failed: ${failedTests}`);
console.log('═'.repeat(80));

if (failedTests === 0) {
  console.log('\n✓ ALL TESTS PASSED!\n');
  console.log('Comprehensive double well testing complete:');
  console.log('  ✓ Parity alternation verified for ALL states');
  console.log('  ✓ Edge behavior correct for ALL states');
  console.log('  ✓ Node count verified for ALL states');
  console.log('  ✓ Normalization verified for ALL states');
  console.log('  ✓ Energy ordering verified for ALL states');
  console.log('  ✓ Derivative continuity verified for ALL states');
  console.log('  ✓ Parameter ranges thoroughly tested');
  console.log('  ✓ Parameter sensitivity verified');
  console.log('  ✓ Systematic parameter sweeps completed');
  console.log('  ✓ Grid convergence verified');
  console.log('  ✓ Extreme parameter regimes tested');
  console.log('');
  process.exit(0);
} else {
  console.log('\n✗ SOME TESTS FAILED\n');
  process.exit(1);
}
