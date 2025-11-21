#!/usr/bin/env node
/**
 * Test suite for QuantumBoundStateSolver improvements:
 * - Logarithmic derivative integration method
 * - Symmetry detection
 * - Parity-based eigenstate finding
 *
 * Usage:
 *   npx tsx --import ./tests/browser-globals.js tests/test-quantum-solver-improvements.ts
 */

import QuantumConstants from '../src/common/model/QuantumConstants.js';
import { QuantumBoundStateSolver } from '../src/common/model/QuantumBoundStateSolver.js';
import { PotentialFunction } from '../src/common/model/PotentialFunction.js';

// Constants
const ELECTRON_MASS = QuantumConstants.ELECTRON_MASS;
const EV_TO_JOULES = QuantumConstants.EV_TO_JOULES;
const NUM_POINTS = 512;

/**
 * Test 1: Symmetry Detection
 */
function testSymmetryDetection() {
  console.log('\n╔════════════════════════════════════════════════════════════════════╗');
  console.log('║   Test 1: Symmetry Detection                                       ║');
  console.log('╚════════════════════════════════════════════════════════════════════╝\n');

  const tests = [
    {
      name: 'Harmonic Oscillator (symmetric)',
      potential: (x: number) => 0.5 * ELECTRON_MASS * (1e15 ** 2) * (x ** 2),
      expectedSymmetric: true
    },
    {
      name: 'Square Well (symmetric)',
      potential: (x: number) => Math.abs(x) < 1e-9 ? 0 : 10 * EV_TO_JOULES,
      expectedSymmetric: true
    },
    {
      name: 'Asymmetric Double Well',
      potential: (x: number) => {
        if (Math.abs(x + 0.8e-9) < 0.5e-9) return 0;
        if (Math.abs(x - 1.2e-9) < 0.5e-9) return 0;
        return 10 * EV_TO_JOULES;
      },
      expectedSymmetric: false
    },
    {
      name: 'Linear Potential (asymmetric)',
      potential: (x: number) => x * 1e10 * EV_TO_JOULES,
      expectedSymmetric: false
    }
  ];

  let passed = 0;
  let failed = 0;

  for (const test of tests) {
    const solver = new QuantumBoundStateSolver(
      ELECTRON_MASS,
      -5e-9,
      5e-9,
      NUM_POINTS,
      test.potential
    );

    const isSymmetric = solver.isPotentialSymmetric();
    const testPassed = isSymmetric === test.expectedSymmetric;

    if (testPassed) {
      console.log(`✓ ${test.name}: ${isSymmetric ? 'symmetric' : 'asymmetric'} (as expected)`);
      passed++;
    } else {
      console.log(`✗ ${test.name}: ${isSymmetric ? 'symmetric' : 'asymmetric'} (expected ${test.expectedSymmetric ? 'symmetric' : 'asymmetric'})`);
      failed++;
    }
  }

  console.log(`\nSymmetry Detection: ${passed}/${passed + failed} tests passed`);
  return { passed, failed };
}

/**
 * Test 2: Inward-Outward Shooting Method
 */
function testInwardOutwardShooting() {
  console.log('\n╔════════════════════════════════════════════════════════════════════╗');
  console.log('║   Test 2: Inward-Outward Shooting Method                           ║');
  console.log('╚════════════════════════════════════════════════════════════════════╝\n');

  // Test with a harmonic oscillator (known analytical solutions)
  const omega = 1e15; // rad/s
  const potential: PotentialFunction = (x: number) =>
    0.5 * ELECTRON_MASS * (omega ** 2) * (x ** 2);

  const xMin = -5e-9;
  const xMax = 5e-9;

  console.log('Testing Harmonic Oscillator (ω = 1e15 rad/s)');
  console.log('Analytical energies: E_n = ℏω(n + 1/2)\n');

  // Create solver
  const solver = new QuantumBoundStateSolver(
    ELECTRON_MASS,
    xMin,
    xMax,
    NUM_POINTS,
    potential
  );

  // Find first 5 states
  const states = solver.findMultipleEigenstates(5);

  // Calculate analytical energies
  const hbar = QuantumConstants.HBAR;
  const analyticalEnergies = [0, 1, 2, 3, 4].map(n =>
    hbar * omega * (n + 0.5) / EV_TO_JOULES
  );

  console.log('State | Numerical (eV) | Analytical (eV) | Error (%) | Match?');
  console.log('─'.repeat(80));

  let passed = 0;
  let failed = 0;

  for (let i = 0; i < Math.min(states.length, analyticalEnergies.length); i++) {
    const eNum = states[i].energy / EV_TO_JOULES;
    const eAna = analyticalEnergies[i];
    const error = Math.abs(eNum - eAna) / eAna * 100;

    // Should match analytical solution to within 1%
    const match = error < 1.0;
    const icon = match ? '✓' : '✗';

    console.log(
      `  ${i}   | ${eNum.toFixed(6).padStart(13)} | ` +
      `${eAna.toFixed(6).padStart(14)} | ` +
      `${error.toFixed(4).padStart(8)} | ${icon}`
    );

    if (match) {
      passed++;
    } else {
      failed++;
    }
  }

  console.log(`\nInward-Outward Shooting: ${passed}/${passed + failed} states match analytical`);
  return { passed, failed };
}

/**
 * Test 3: Matching Point Optimization
 */
function testMatchingPointOptimization() {
  console.log('\n╔════════════════════════════════════════════════════════════════════╗');
  console.log('║   Test 3: Matching Point Optimization                              ║');
  console.log('╚════════════════════════════════════════════════════════════════════╝\n');

  // Double well where matching point matters
  const wellWidth = 1.5e-9;
  const wellSeparation = 1e-9;
  const wellDepth = 10 * EV_TO_JOULES;

  const halfWidth = wellWidth / 2;
  const leftCenter = -wellSeparation / 2;
  const rightCenter = wellSeparation / 2;

  const potential: PotentialFunction = (x: number) => {
    const inLeftWell = Math.abs(x - leftCenter) <= halfWidth;
    const inRightWell = Math.abs(x - rightCenter) <= halfWidth;
    return (inLeftWell || inRightWell) ? 0 : wellDepth;
  };

  const totalExtent = wellSeparation + wellWidth + 4e-9;
  const xMin = -totalExtent / 2;
  const xMax = totalExtent / 2;

  console.log('Testing Double Square Well with Adaptive Matching\n');

  try {
    const solver = new QuantumBoundStateSolver(
      ELECTRON_MASS,
      xMin,
      xMax,
      NUM_POINTS,
      potential,
      {
        adaptiveMatching: true
      }
    );

    const states = solver.findMultipleEigenstates(10);
    console.log(`✓ Found ${states.length} bound states with adaptive matching`);

    // Check that energies are reasonable (below barrier)
    let passed = 0;
    let failed = 0;

    for (const state of states) {
      const energy = state.energy / EV_TO_JOULES;
      const belowBarrier = energy < 10; // Should be below barrier
      if (belowBarrier) {
        passed++;
      } else {
        failed++;
        console.log(`  ✗ State with E=${energy.toFixed(3)} eV is above barrier!`);
      }
    }

    console.log(`\nMatching Optimization: ${passed}/${passed + failed} states below barrier`);
    return { passed, failed };
  } catch (error) {
    console.log(`✗ Error: ${error instanceof Error ? error.message : String(error)}`);
    return { passed: 0, failed: 1 };
  }
}

/**
 * Test 4: Node Counting Accuracy
 */
function testNodeCounting() {
  console.log('\n╔════════════════════════════════════════════════════════════════════╗');
  console.log('║   Test 4: Node Counting Accuracy                                   ║');
  console.log('╚════════════════════════════════════════════════════════════════════╝\n');

  // Create a symmetric square well
  const wellDepth = 10 * EV_TO_JOULES;
  const wellWidth = 2e-9;
  const potential: PotentialFunction = (x: number) =>
    Math.abs(x) < wellWidth / 2 ? 0 : wellDepth;

  const xMin = -5e-9;
  const xMax = 5e-9;

  console.log('Testing Symmetric Square Well (10 eV depth, 2nm width)\n');

  try {
    const solver = new QuantumBoundStateSolver(
      ELECTRON_MASS,
      xMin,
      xMax,
      NUM_POINTS,
      potential
    );

    const states = solver.findMultipleEigenstates(10);
    console.log(`Found ${states.length} bound states\n`);

    console.log('State | Energy (eV) | Nodes | Expected Parity | Correct?');
    console.log('─'.repeat(70));

    let passed = 0;
    let failed = 0;

    for (let i = 0; i < states.length; i++) {
      const state = states[i];
      const energy = state.energy / EV_TO_JOULES;
      const nodes = state.nodes;
      const expectedNodes = i;
      const expectedParity = i % 2 === 0 ? 'Even' : 'Odd';
      const correct = nodes === expectedNodes;
      const icon = correct ? '✓' : '✗';

      console.log(
        `  ${i}   | ${energy.toFixed(6).padStart(10)} | ` +
        `${nodes.toString().padStart(4)} | ${expectedParity.padStart(14)} | ${icon}`
      );

      if (correct) {
        passed++;
      } else {
        failed++;
      }
    }

    console.log(`\nNode Counting: ${passed}/${passed + failed} states have correct node count`);
    return { passed, failed };
  } catch (error) {
    console.log(`✗ Error during node counting test: ${error instanceof Error ? error.message : String(error)}`);
    return { passed: 0, failed: 1 };
  }
}

/**
 * Main test runner
 */
function runTests() {
  console.log('\n╔════════════════════════════════════════════════════════════════════╗');
  console.log('║   QuantumBoundStateSolver - Improvements Test Suite                ║');
  console.log('║   Testing: Inward-Outward, Symmetry Detection, Parity              ║');
  console.log('╚════════════════════════════════════════════════════════════════════╝');

  const results = [
    testSymmetryDetection(),
    testInwardOutwardShooting(),
    testMatchingPointOptimization(),
    testNodeCounting()
  ];

  const totalPassed = results.reduce((sum, r) => sum + r.passed, 0);
  const totalFailed = results.reduce((sum, r) => sum + r.failed, 0);
  const totalTests = totalPassed + totalFailed;

  console.log('\n' + '═'.repeat(80));
  console.log('OVERALL SUMMARY');
  console.log('═'.repeat(80));
  console.log(`Total tests: ${totalTests}`);
  console.log(`Passed: ${totalPassed} (${(totalPassed / totalTests * 100).toFixed(1)}%)`);
  console.log(`Failed: ${totalFailed} (${(totalFailed / totalTests * 100).toFixed(1)}%)`);
  console.log('═'.repeat(80));

  if (totalFailed === 0) {
    console.log('\n🎉 All tests passed! The solver improvements are working correctly.\n');
  } else {
    console.log(`\n⚠ ${totalFailed} test(s) failed. Review the output above for details.\n`);
  }

  // Exit with appropriate code
  process.exit(totalFailed === 0 ? 0 : 1);
}

// Run tests
runTests();
