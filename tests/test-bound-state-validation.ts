/**
 * Tests for QuantumBoundStateSolver bound state validation.
 *
 * Verifies that the solver correctly rejects eigenstates that are not valid
 * bound states (where energy exceeds the maximum potential).
 */

import { QuantumBoundStateSolver, InvalidBoundStateException, solveQuantumBound } from "../src/common/model/QuantumBoundStateSolver.js";
import QuantumConstants from "../src/common/model/QuantumConstants.js";
import { GridConfig, PotentialFunction } from "../src/common/model/PotentialFunction.js";

/**
 * Test result structure
 */
interface TestResult {
  testName: string;
  passed: boolean;
  details: string[];
}

/**
 * Test that the solver correctly identifies and rejects invalid bound states
 * for a finite square well potential (which has a limited number of bound states).
 */
function testFiniteWellBoundStateLimit(): TestResult {
  const details: string[] = [];
  const testName = "Finite Well Bound State Limit";

  const mass = QuantumConstants.ELECTRON_MASS;

  // Create a shallow finite square well that has only a few bound states
  const wellWidth = 1e-9;  // 1 nm
  const wellDepth = 2 * QuantumConstants.EV_TO_JOULES;  // 2 eV - shallow well

  const gridConfig: GridConfig = {
    xMin: -3e-9,
    xMax: 3e-9,
    numPoints: 200
  };

  const potential: PotentialFunction = (x: number) => {
    const halfWidth = wellWidth / 2;
    return (x >= -halfWidth && x <= halfWidth) ? -wellDepth : 0;
  };

  try {
    const solver = new QuantumBoundStateSolver(
      mass,
      gridConfig.xMin,
      gridConfig.xMax,
      gridConfig.numPoints,
      potential
    );

    const maxPotential = solver.getMaxPotential();
    details.push(`Maximum potential: ${maxPotential / QuantumConstants.EV_TO_JOULES} eV`);

    // Request many states - the solver should stop when bound states are exhausted
    const requestedStates = 10;
    const states = solver.findMultipleEigenstates(requestedStates);

    details.push(`Requested ${requestedStates} states, found ${states.length} bound states`);

    // Verify all returned states have energy below max potential
    let allValid = true;
    for (let i = 0; i < states.length; i++) {
      const state = states[i];
      const energyEV = state.energy / QuantumConstants.EV_TO_JOULES;
      const isValid = state.energy < maxPotential;

      if (!isValid) {
        allValid = false;
        details.push(`  State ${i + 1}: E = ${energyEV.toFixed(6)} eV - INVALID (exceeds max potential)`);
      } else {
        details.push(`  State ${i + 1}: E = ${energyEV.toFixed(6)} eV - valid bound state`);
      }
    }

    // A shallow 2 eV well of 1 nm width should have only 2-3 bound states
    const expectedMaxStates = 5;
    const hasReasonableCount = states.length > 0 && states.length <= expectedMaxStates;

    if (!hasReasonableCount) {
      details.push(`Expected 1-${expectedMaxStates} bound states for this shallow well, got ${states.length}`);
    }

    const passed = allValid && hasReasonableCount;

    return {
      testName,
      passed,
      details
    };
  } catch (error) {
    details.push(`ERROR: ${error}`);
    return {
      testName,
      passed: false,
      details
    };
  }
}

/**
 * Test that InvalidBoundStateException is thrown when requesting a specific
 * state that doesn't exist.
 */
function testInvalidStateException(): TestResult {
  const details: string[] = [];
  const testName = "Invalid Bound State Exception";

  const mass = QuantumConstants.ELECTRON_MASS;

  // Create a very shallow well with only 1-2 bound states
  const wellWidth = 0.5e-9;  // 0.5 nm
  const wellDepth = 1 * QuantumConstants.EV_TO_JOULES;  // 1 eV - very shallow

  const gridConfig: GridConfig = {
    xMin: -2e-9,
    xMax: 2e-9,
    numPoints: 200
  };

  const potential: PotentialFunction = (x: number) => {
    const halfWidth = wellWidth / 2;
    return (x >= -halfWidth && x <= halfWidth) ? -wellDepth : 0;
  };

  try {
    const solver = new QuantumBoundStateSolver(
      mass,
      gridConfig.xMin,
      gridConfig.xMax,
      gridConfig.numPoints,
      potential
    );

    const maxPotential = solver.getMaxPotential();
    details.push(`Maximum potential: ${maxPotential / QuantumConstants.EV_TO_JOULES} eV`);

    // First find the ground state (should succeed)
    let _groundStateFound = false;
    try {
      const groundState = solver.findEigenstate(1);
      _groundStateFound = true;
      details.push(`Ground state energy: ${groundState.energy / QuantumConstants.EV_TO_JOULES} eV`);
    } catch (error) {
      if (error instanceof InvalidBoundStateException) {
        details.push(`No ground state exists (very shallow well)`);
      } else {
        throw error;
      }
    }

    // Try to find a high excited state (should fail)
    let exceptionThrown = false;
    try {
      const highState = solver.findEigenstate(10);  // Request 10th state
      details.push(`Unexpectedly found state 10 with energy: ${highState.energy / QuantumConstants.EV_TO_JOULES} eV`);
    } catch (error) {
      if (error instanceof InvalidBoundStateException) {
        exceptionThrown = true;
        details.push(`Correctly threw InvalidBoundStateException for state ${error.stateNumber}`);
        details.push(`  Energy: ${error.energy / QuantumConstants.EV_TO_JOULES} eV`);
        details.push(`  Max potential: ${error.maxPotential / QuantumConstants.EV_TO_JOULES} eV`);
      } else {
        throw error;
      }
    }

    const passed = exceptionThrown;

    return {
      testName,
      passed,
      details
    };
  } catch (error) {
    details.push(`ERROR: ${error}`);
    return {
      testName,
      passed: false,
      details
    };
  }
}

/**
 * Test the functional API (solveQuantumBound) handles bound state limits correctly.
 */
function testFunctionalAPIBoundStateLimit(): TestResult {
  const details: string[] = [];
  const testName = "Functional API Bound State Limit";

  const mass = QuantumConstants.ELECTRON_MASS;

  // Create a double well potential with a barrier
  // This should have limited bound states below the barrier
  const wellWidth = 1e-9;
  const barrierWidth = 0.3e-9;
  const wellDepth = 5 * QuantumConstants.EV_TO_JOULES;
  const barrierHeight = 3 * QuantumConstants.EV_TO_JOULES;  // Barrier is 2 eV below continuum

  const gridConfig: GridConfig = {
    xMin: -3e-9,
    xMax: 3e-9,
    numPoints: 200
  };

  const potential: PotentialFunction = (x: number) => {
    const halfBarrier = barrierWidth / 2;
    const wellEdge = halfBarrier + wellWidth;

    if (Math.abs(x) < halfBarrier) {
      return -wellDepth + barrierHeight;  // Barrier region
    } else if (Math.abs(x) < wellEdge) {
      return -wellDepth;  // Well regions
    } else {
      return 0;  // Outside
    }
  };

  try {
    // Request many states
    const result = solveQuantumBound(potential, mass, 20, gridConfig);

    details.push(`Requested 20 states, found ${result.energies.length} bound states`);

    // Calculate max potential
    const positions = result.xGrid;
    let maxPot = -Infinity;
    for (const x of positions) {
      maxPot = Math.max(maxPot, potential(x));
    }

    details.push(`Maximum potential: ${maxPot / QuantumConstants.EV_TO_JOULES} eV`);

    // Verify all energies are below max potential
    let allValid = true;
    for (let i = 0; i < result.energies.length; i++) {
      const energy = result.energies[i];
      const energyEV = energy / QuantumConstants.EV_TO_JOULES;

      if (energy >= maxPot) {
        allValid = false;
        details.push(`  E_${i}: ${energyEV.toFixed(6)} eV - INVALID`);
      } else {
        details.push(`  E_${i}: ${energyEV.toFixed(6)} eV - valid`);
      }
    }

    // Should find some states but not all 20
    const hasReasonableCount = result.energies.length > 0 && result.energies.length < 20;

    if (result.energies.length === 20) {
      details.push(`Warning: Found all 20 requested states - validation may not be working`);
    }

    const passed = allValid && hasReasonableCount;

    return {
      testName,
      passed,
      details
    };
  } catch (error) {
    details.push(`ERROR: ${error}`);
    return {
      testName,
      passed: false,
      details
    };
  }
}

/**
 * Test with a harmonic oscillator (should have many bound states - tests that
 * validation doesn't reject valid states).
 */
function testHarmonicOscillatorManyStates(): TestResult {
  const details: string[] = [];
  const testName = "Harmonic Oscillator Many States";

  const mass = QuantumConstants.ELECTRON_MASS;
  const omega = 1e15;
  const springConstant = mass * omega * omega;

  const gridConfig: GridConfig = {
    xMin: -5e-9,
    xMax: 5e-9,
    numPoints: 200
  };

  const potential: PotentialFunction = (x: number) => 0.5 * springConstant * x * x;

  try {
    const solver = new QuantumBoundStateSolver(
      mass,
      gridConfig.xMin,
      gridConfig.xMax,
      gridConfig.numPoints,
      potential
    );

    const maxPotential = solver.getMaxPotential();
    details.push(`Maximum potential: ${maxPotential / QuantumConstants.EV_TO_JOULES} eV`);

    // Request 10 states - harmonic oscillator should support many
    const states = solver.findMultipleEigenstates(10);

    details.push(`Requested 10 states, found ${states.length} bound states`);

    // Verify energies follow E_n = (n + 1/2) * hbar * omega pattern
    const hbar = QuantumConstants.HBAR;
    let allValid = true;

    for (let i = 0; i < states.length; i++) {
      const state = states[i];
      const energyEV = state.energy / QuantumConstants.EV_TO_JOULES;
      const expectedEnergy = (i + 0.5) * hbar * omega;
      const expectedEV = expectedEnergy / QuantumConstants.EV_TO_JOULES;
      const error = Math.abs((state.energy - expectedEnergy) / expectedEnergy) * 100;

      if (error > 5) {  // Allow 5% error for shooting method
        allValid = false;
        details.push(`  State ${i + 1}: E = ${energyEV.toFixed(6)} eV (expected ${expectedEV.toFixed(6)} eV) - error ${error.toFixed(2)}%`);
      } else {
        details.push(`  State ${i + 1}: E = ${energyEV.toFixed(6)} eV - OK (error ${error.toFixed(2)}%)`);
      }
    }

    // Should find all 10 states (harmonic oscillator has infinite bound states in principle,
    // but limited by grid extent)
    const foundAllStates = states.length === 10;

    if (!foundAllStates) {
      details.push(`Expected to find 10 states, only found ${states.length}`);
    }

    const passed = allValid && foundAllStates;

    return {
      testName,
      passed,
      details
    };
  } catch (error) {
    details.push(`ERROR: ${error}`);
    return {
      testName,
      passed: false,
      details
    };
  }
}

/**
 * Run all bound state validation tests
 */
export function runBoundStateValidationTests(): void {
  console.log("========================================");
  console.log("Bound State Validation Tests");
  console.log("========================================");
  console.log("Testing that QuantumBoundStateSolver correctly");
  console.log("rejects states with energy above max potential");
  console.log("");

  const results: TestResult[] = [];

  results.push(testFiniteWellBoundStateLimit());
  results.push(testInvalidStateException());
  results.push(testFunctionalAPIBoundStateLimit());
  results.push(testHarmonicOscillatorManyStates());

  // Print results
  for (const result of results) {
    console.log(`\n=== ${result.testName} ===`);
    console.log(`Status: ${result.passed ? "PASSED" : "FAILED"}`);
    for (const detail of result.details) {
      console.log(detail);
    }
  }

  // Summary
  const passed = results.filter(r => r.passed).length;
  const total = results.length;

  console.log("\n========================================");
  console.log("Summary");
  console.log("========================================");
  console.log(`Tests: ${passed}/${total} passed`);

  if (passed === total) {
    console.log("\nAll bound state validation tests passed!");
  } else {
    console.log("\nSome tests failed - review details above.");
  }

  console.log("========================================\n");
}

// Run tests if executed directly
runBoundStateValidationTests();

export default { runBoundStateValidationTests };
