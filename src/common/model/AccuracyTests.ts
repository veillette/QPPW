/**
 * Accuracy validation tests for DVR, Spectral, Matrix Numerov, and FGH methods.
 *
 * This file tests the numerical methods (DVR, Spectral, Matrix Numerov, and FGH) against known
 * analytical solutions to verify they produce results within 1% accuracy.
 *
 * Test cases:
 * 1. Harmonic Oscillator - exact analytical solution: E_n = ℏω(n + 1/2)
 * 2. Infinite Square Well - exact analytical solution: E_n = n²π²ℏ²/(2mL²)
 *
 * To run these tests:
 * import { runAccuracyTests } from './common/model/AccuracyTests.js';
 * runAccuracyTests();
 */

import { solveDVR } from "./DVRSolver.js";
import { solveSpectral } from "./SpectralSolver.js";
import { solveMatrixNumerov } from "./MatrixNumerovSolver.js";
import { solveFGH } from "./FGHSolver.js";
import { solveHarmonicOscillator } from "./analytical-solutions/harmonic-oscillator.js";
import { solveInfiniteWell } from "./analytical-solutions/infinite-square-well.js";
import QuantumConstants from "./QuantumConstants.js";
import { GridConfig, PotentialFunction } from "./PotentialFunction.js";
import Schrodinger1DSolver from "./Schrodinger1DSolver.js";
import qppw from "../../QPPWNamespace.js";

/**
 * Test result structure
 */
interface TestResult {
  testName: string;
  method: string;
  passed: boolean;
  maxError: number;
  details: string[];
}

/**
 * Calculate percentage error between numerical and analytical values
 */
function percentageError(numerical: number, analytical: number): number {
  return Math.abs((numerical - analytical) / analytical) * 100;
}

/**
 * Test DVR method against harmonic oscillator analytical solution
 */
function testDVRHarmonicOscillator(): TestResult {
  const testName = "DVR - Harmonic Oscillator";
  const details: string[] = [];

  // Setup parameters
  const mass = QuantumConstants.ELECTRON_MASS;
  const omega = 1e15; // Angular frequency (rad/s)
  const springConstant = mass * omega * omega; // k = mω²
  const numStates = 5;

  const gridConfig: GridConfig = {
    xMin: -5e-9,
    xMax: 5e-9,
    numPoints: 200,
  };

  // Create harmonic oscillator potential: V(x) = (1/2) * k * x²
  const potential: PotentialFunction = (x: number) => 0.5 * springConstant * x * x;

  // Get numerical solution using DVR
  const numericalResult = solveDVR(potential, mass, numStates, gridConfig);

  // Get analytical solution
  const analyticalResult = solveHarmonicOscillator(
    springConstant,
    mass,
    numStates,
    gridConfig
  );

  // Compare energies
  let maxError = 0;
  let allPassed = true;

  details.push(`Testing ${numStates} energy levels:`);

  for (let n = 0; n < numStates; n++) {
    const E_numerical = numericalResult.energies[n];
    const E_analytical = analyticalResult.energies[n];
    const error = percentageError(E_numerical, E_analytical);

    const E_numerical_eV = Schrodinger1DSolver.joulesToEV(E_numerical);
    const E_analytical_eV = Schrodinger1DSolver.joulesToEV(E_analytical);

    const passed = error < 1.0; // 1% tolerance
    allPassed = allPassed && passed;
    maxError = Math.max(maxError, error);

    const status = passed ? "✓ PASS" : "✗ FAIL";
    details.push(
      `  E_${n}: ${status} - Numerical: ${E_numerical_eV.toFixed(6)} eV, ` +
      `Analytical: ${E_analytical_eV.toFixed(6)} eV, Error: ${error.toFixed(4)}%`
    );
  }

  return {
    testName,
    method: "DVR",
    passed: allPassed,
    maxError,
    details,
  };
}

/**
 * Test Spectral method against harmonic oscillator analytical solution
 */
function testSpectralHarmonicOscillator(): TestResult {
  const testName = "Spectral - Harmonic Oscillator";
  const details: string[] = [];

  // Setup parameters
  const mass = QuantumConstants.ELECTRON_MASS;
  const omega = 1e15; // Angular frequency (rad/s)
  const springConstant = mass * omega * omega; // k = mω²
  const numStates = 5;

  const gridConfig: GridConfig = {
    xMin: -5e-9,
    xMax: 5e-9,
    numPoints: 200,
  };

  // Create harmonic oscillator potential: V(x) = (1/2) * k * x²
  const potential: PotentialFunction = (x: number) => 0.5 * springConstant * x * x;

  // Get numerical solution using Spectral method
  const numericalResult = solveSpectral(potential, mass, numStates, gridConfig);

  // Get analytical solution
  const analyticalResult = solveHarmonicOscillator(
    springConstant,
    mass,
    numStates,
    gridConfig
  );

  // Compare energies
  let maxError = 0;
  let allPassed = true;

  details.push(`Testing ${numStates} energy levels:`);

  for (let n = 0; n < numStates; n++) {
    const E_numerical = numericalResult.energies[n];
    const E_analytical = analyticalResult.energies[n];
    const error = percentageError(E_numerical, E_analytical);

    const E_numerical_eV = Schrodinger1DSolver.joulesToEV(E_numerical);
    const E_analytical_eV = Schrodinger1DSolver.joulesToEV(E_analytical);

    const passed = error < 1.0; // 1% tolerance
    allPassed = allPassed && passed;
    maxError = Math.max(maxError, error);

    const status = passed ? "✓ PASS" : "✗ FAIL";
    details.push(
      `  E_${n}: ${status} - Numerical: ${E_numerical_eV.toFixed(6)} eV, ` +
      `Analytical: ${E_analytical_eV.toFixed(6)} eV, Error: ${error.toFixed(4)}%`
    );
  }

  return {
    testName,
    method: "Spectral",
    passed: allPassed,
    maxError,
    details,
  };
}

/**
 * Test DVR method against infinite square well analytical solution
 */
function testDVRInfiniteWell(): TestResult {
  const testName = "DVR - Infinite Square Well";
  const details: string[] = [];

  // Setup parameters
  const wellWidth = 1e-9; // 1 nm
  const mass = QuantumConstants.ELECTRON_MASS;
  const numStates = 5;

  const gridConfig: GridConfig = {
    xMin: -wellWidth,
    xMax: wellWidth,
    numPoints: 150,
  };

  // Create infinite well potential centered at x=0
  const V_high = 1e10; // Very high potential outside well (effectively infinite)
  const potential: PotentialFunction = (x: number) => {
    const halfWidth = wellWidth / 2;
    return (x >= -halfWidth && x <= halfWidth) ? 0 : V_high;
  };

  // Get numerical solution using DVR
  const numericalResult = solveDVR(potential, mass, numStates, gridConfig);

  // Get analytical solution
  const analyticalResult = solveInfiniteWell(
    wellWidth,
    mass,
    numStates,
    gridConfig
  );

  // Compare energies
  let maxError = 0;
  let allPassed = true;

  details.push(`Testing ${numStates} energy levels:`);

  for (let n = 0; n < Math.min(numStates, numericalResult.energies.length); n++) {
    const E_numerical = numericalResult.energies[n];
    const E_analytical = analyticalResult.energies[n];
    const error = percentageError(E_numerical, E_analytical);

    const E_numerical_eV = Schrodinger1DSolver.joulesToEV(E_numerical);
    const E_analytical_eV = Schrodinger1DSolver.joulesToEV(E_analytical);

    const passed = error < 1.0; // 1% tolerance
    allPassed = allPassed && passed;
    maxError = Math.max(maxError, error);

    const status = passed ? "✓ PASS" : "✗ FAIL";
    details.push(
      `  E_${n + 1}: ${status} - Numerical: ${E_numerical_eV.toFixed(6)} eV, ` +
      `Analytical: ${E_analytical_eV.toFixed(6)} eV, Error: ${error.toFixed(4)}%`
    );
  }

  return {
    testName,
    method: "DVR",
    passed: allPassed,
    maxError,
    details,
  };
}

/**
 * Test FGH method against harmonic oscillator analytical solution
 */
function testFGHHarmonicOscillator(): TestResult {
  const testName = "FGH - Harmonic Oscillator";
  const details: string[] = [];

  // Setup parameters
  const mass = QuantumConstants.ELECTRON_MASS;
  const omega = 1e15; // Angular frequency (rad/s)
  const springConstant = mass * omega * omega; // k = mω²
  const numStates = 5;

  const gridConfig: GridConfig = {
    xMin: -5e-9,
    xMax: 5e-9,
    numPoints: 256, // Must be power of 2 for radix-2 FFT
  };

  // Create harmonic oscillator potential: V(x) = (1/2) * k * x²
  const potential: PotentialFunction = (x: number) => 0.5 * springConstant * x * x;

  // Get numerical solution using FGH method
  const numericalResult = solveFGH(potential, mass, numStates, gridConfig);

  // Get analytical solution
  const analyticalResult = solveHarmonicOscillator(
    springConstant,
    mass,
    numStates,
    gridConfig
  );

  // Compare energies
  let maxError = 0;
  let allPassed = true;

  details.push(`Testing ${numStates} energy levels:`);

  for (let n = 0; n < numStates; n++) {
    const E_numerical = numericalResult.energies[n];
    const E_analytical = analyticalResult.energies[n];
    const error = percentageError(E_numerical, E_analytical);

    const E_numerical_eV = Schrodinger1DSolver.joulesToEV(E_numerical);
    const E_analytical_eV = Schrodinger1DSolver.joulesToEV(E_analytical);

    const passed = error < 1.0; // 1% tolerance
    allPassed = allPassed && passed;
    maxError = Math.max(maxError, error);

    const status = passed ? "✓ PASS" : "✗ FAIL";
    details.push(
      `  E_${n}: ${status} - Numerical: ${E_numerical_eV.toFixed(6)} eV, ` +
      `Analytical: ${E_analytical_eV.toFixed(6)} eV, Error: ${error.toFixed(4)}%`
    );
  }

  return {
    testName,
    method: "FGH",
    passed: allPassed,
    maxError,
    details,
  };
}

/**
 * Test FGH method against infinite square well analytical solution
 */
function testFGHInfiniteWell(): TestResult {
  const testName = "FGH - Infinite Square Well";
  const details: string[] = [];

  // Setup parameters
  const wellWidth = 1e-9; // 1 nm
  const mass = QuantumConstants.ELECTRON_MASS;
  const numStates = 5;

  const gridConfig: GridConfig = {
    xMin: -wellWidth,
    xMax: wellWidth,
    numPoints: 128, // Must be power of 2 for radix-2 FFT
  };

  // Create infinite well potential centered at x=0
  const V_high = 1e10; // Very high potential outside well (effectively infinite)
  const potential: PotentialFunction = (x: number) => {
    const halfWidth = wellWidth / 2;
    return (x >= -halfWidth && x <= halfWidth) ? 0 : V_high;
  };

  // Get numerical solution using FGH method
  const numericalResult = solveFGH(potential, mass, numStates, gridConfig);

  // Get analytical solution
  const analyticalResult = solveInfiniteWell(
    wellWidth,
    mass,
    numStates,
    gridConfig
  );

  // Compare energies
  let maxError = 0;
  let allPassed = true;

  details.push(`Testing ${numStates} energy levels:`);

  for (let n = 0; n < Math.min(numStates, numericalResult.energies.length); n++) {
    const E_numerical = numericalResult.energies[n];
    const E_analytical = analyticalResult.energies[n];
    const error = percentageError(E_numerical, E_analytical);

    const E_numerical_eV = Schrodinger1DSolver.joulesToEV(E_numerical);
    const E_analytical_eV = Schrodinger1DSolver.joulesToEV(E_analytical);

    const passed = error < 1.0; // 1% tolerance
    allPassed = allPassed && passed;
    maxError = Math.max(maxError, error);

    const status = passed ? "✓ PASS" : "✗ FAIL";
    details.push(
      `  E_${n + 1}: ${status} - Numerical: ${E_numerical_eV.toFixed(6)} eV, ` +
      `Analytical: ${E_analytical_eV.toFixed(6)} eV, Error: ${error.toFixed(4)}%`
    );
  }

  return {
    testName,
    method: "FGH",
    passed: allPassed,
    maxError,
    details,
  };
}

/**
 * Test Spectral method against infinite square well analytical solution
 */
function testSpectralInfiniteWell(): TestResult {
  const testName = "Spectral - Infinite Square Well";
  const details: string[] = [];

  // Setup parameters
  const wellWidth = 1e-9; // 1 nm
  const mass = QuantumConstants.ELECTRON_MASS;
  const numStates = 5;

  const gridConfig: GridConfig = {
    xMin: -wellWidth / 2,
    xMax: wellWidth / 2,
    numPoints: 150,
  };

  // Create infinite well potential centered at x=0
  // For spectral method, we use boundary conditions at xMin and xMax
  const potential: PotentialFunction = () => 0;

  // Get numerical solution using Spectral method
  const numericalResult = solveSpectral(potential, mass, numStates, gridConfig);

  // Get analytical solution
  const analyticalResult = solveInfiniteWell(
    wellWidth,
    mass,
    numStates,
    gridConfig
  );

  // Compare energies
  let maxError = 0;
  let allPassed = true;

  details.push(`Testing ${numStates} energy levels:`);

  for (let n = 0; n < Math.min(numStates, numericalResult.energies.length); n++) {
    const E_numerical = numericalResult.energies[n];
    const E_analytical = analyticalResult.energies[n];
    const error = percentageError(E_numerical, E_analytical);

    const E_numerical_eV = Schrodinger1DSolver.joulesToEV(E_numerical);
    const E_analytical_eV = Schrodinger1DSolver.joulesToEV(E_analytical);

    const passed = error < 1.0; // 1% tolerance
    allPassed = allPassed && passed;
    maxError = Math.max(maxError, error);

    const status = passed ? "✓ PASS" : "✗ FAIL";
    details.push(
      `  E_${n + 1}: ${status} - Numerical: ${E_numerical_eV.toFixed(6)} eV, ` +
      `Analytical: ${E_analytical_eV.toFixed(6)} eV, Error: ${error.toFixed(4)}%`
    );
  }

  return {
    testName,
    method: "Spectral",
    passed: allPassed,
    maxError,
    details,
  };
}

/**
 * Test Matrix Numerov method against harmonic oscillator analytical solution
 */
function testMatrixNumerovHarmonicOscillator(): TestResult {
  const testName = "Matrix Numerov - Harmonic Oscillator";
  const details: string[] = [];

  // Setup parameters
  const mass = QuantumConstants.ELECTRON_MASS;
  const omega = 1e15; // Angular frequency (rad/s)
  const springConstant = mass * omega * omega; // k = mω²
  const numStates = 5;

  const gridConfig: GridConfig = {
    xMin: -5e-9,
    xMax: 5e-9,
    numPoints: 200,
  };

  // Create harmonic oscillator potential: V(x) = (1/2) * k * x²
  const potential: PotentialFunction = (x: number) => 0.5 * springConstant * x * x;

  // Get numerical solution using Matrix Numerov
  const numericalResult = solveMatrixNumerov(potential, mass, numStates, gridConfig);

  // Get analytical solution
  const analyticalResult = solveHarmonicOscillator(
    springConstant,
    mass,
    numStates,
    gridConfig
  );

  // Compare energies
  let maxError = 0;
  let allPassed = true;

  details.push(`Testing ${numStates} energy levels:`);

  for (let n = 0; n < numStates; n++) {
    const E_numerical = numericalResult.energies[n];
    const E_analytical = analyticalResult.energies[n];
    const error = percentageError(E_numerical, E_analytical);

    const E_numerical_eV = Schrodinger1DSolver.joulesToEV(E_numerical);
    const E_analytical_eV = Schrodinger1DSolver.joulesToEV(E_analytical);

    const passed = error < 1.0; // 1% tolerance
    allPassed = allPassed && passed;
    maxError = Math.max(maxError, error);

    const status = passed ? "✓ PASS" : "✗ FAIL";
    details.push(
      `  E_${n}: ${status} - Numerical: ${E_numerical_eV.toFixed(6)} eV, ` +
      `Analytical: ${E_analytical_eV.toFixed(6)} eV, Error: ${error.toFixed(4)}%`
    );
  }

  return {
    testName,
    method: "Matrix Numerov",
    passed: allPassed,
    maxError,
    details,
  };
}

/**
 * Test Matrix Numerov method against infinite square well analytical solution
 */
function testMatrixNumerovInfiniteWell(): TestResult {
  const testName = "Matrix Numerov - Infinite Square Well";
  const details: string[] = [];

  // Setup parameters
  const wellWidth = 1e-9; // 1 nm
  const mass = QuantumConstants.ELECTRON_MASS;
  const numStates = 5;

  const gridConfig: GridConfig = {
    xMin: -wellWidth,
    xMax: wellWidth,
    numPoints: 150,
  };

  // Create infinite well potential centered at x=0
  const V_high = 1e10; // Very high potential outside well (effectively infinite)
  const potential: PotentialFunction = (x: number) => {
    const halfWidth = wellWidth / 2;
    return (x >= -halfWidth && x <= halfWidth) ? 0 : V_high;
  };

  // Get numerical solution using Matrix Numerov
  const numericalResult = solveMatrixNumerov(potential, mass, numStates, gridConfig);

  // Get analytical solution
  const analyticalResult = solveInfiniteWell(
    wellWidth,
    mass,
    numStates,
    gridConfig
  );

  // Compare energies
  let maxError = 0;
  let allPassed = true;

  details.push(`Testing ${numStates} energy levels:`);

  for (let n = 0; n < Math.min(numStates, numericalResult.energies.length); n++) {
    const E_numerical = numericalResult.energies[n];
    const E_analytical = analyticalResult.energies[n];
    const error = percentageError(E_numerical, E_analytical);

    const E_numerical_eV = Schrodinger1DSolver.joulesToEV(E_numerical);
    const E_analytical_eV = Schrodinger1DSolver.joulesToEV(E_analytical);

    const passed = error < 1.0; // 1% tolerance
    allPassed = allPassed && passed;
    maxError = Math.max(maxError, error);

    const status = passed ? "✓ PASS" : "✗ FAIL";
    details.push(
      `  E_${n + 1}: ${status} - Numerical: ${E_numerical_eV.toFixed(6)} eV, ` +
      `Analytical: ${E_analytical_eV.toFixed(6)} eV, Error: ${error.toFixed(4)}%`
    );
  }

  return {
    testName,
    method: "Matrix Numerov",
    passed: allPassed,
    maxError,
    details,
  };
}

/**
 * Print test results
 */
function printTestResult(result: TestResult): void {
  const header = `\n=== ${result.testName} ===`;
  const status = result.passed ? "✓ PASSED" : "✗ FAILED";
  const maxErrorStr = `Maximum error: ${result.maxError.toFixed(4)}%`;

  console.log(header);
  console.log(`Status: ${status}`);
  console.log(maxErrorStr);
  console.log("");

  result.details.forEach(detail => console.log(detail));
}

/**
 * Run all accuracy tests
 */
export function runAccuracyTests(): void {
  console.log("========================================");
  console.log("Numerical Method Accuracy Tests");
  console.log("========================================");
  console.log("Tolerance: 1% error from analytical solutions");
  console.log("");

  const results: TestResult[] = [];

  // Run all tests
  results.push(testDVRHarmonicOscillator());
  results.push(testMatrixNumerovHarmonicOscillator());
  results.push(testSpectralHarmonicOscillator());
  results.push(testFGHHarmonicOscillator());
  results.push(testDVRInfiniteWell());
  results.push(testMatrixNumerovInfiniteWell());
  results.push(testSpectralInfiniteWell());
  results.push(testFGHInfiniteWell());

  // Print individual results
  results.forEach(printTestResult);

  // Print summary
  console.log("\n========================================");
  console.log("Test Summary");
  console.log("========================================");

  const totalTests = results.length;
  const passedTests = results.filter(r => r.passed).length;
  const failedTests = totalTests - passedTests;

  console.log(`Total tests: ${totalTests}`);
  console.log(`Passed: ${passedTests}`);
  console.log(`Failed: ${failedTests}`);

  if (failedTests === 0) {
    console.log("\n✓ All tests passed!");
    console.log("All numerical methods (DVR, Matrix Numerov, Spectral, and FGH) produce results within 1% of analytical solutions.");
  } else {
    console.log("\n✗ Some tests failed!");
    console.log("Review the details above for failed tests.");
  }

  console.log("========================================\n");
}

/**
 * Run a quick accuracy check (fewer states for faster testing)
 */
export function runQuickAccuracyCheck(): void {
  console.log("\n=== Quick Accuracy Check ===");
  console.log("Testing DVR, Matrix Numerov, Spectral, and FGH methods with harmonic oscillator...\n");

  const mass = QuantumConstants.ELECTRON_MASS;
  const omega = 1e15;
  const springConstant = mass * omega * omega;
  const numStates = 3;

  const gridConfig: GridConfig = {
    xMin: -5e-9,
    xMax: 5e-9,
    numPoints: 150,
  };

  const gridConfigFGH: GridConfig = {
    xMin: -5e-9,
    xMax: 5e-9,
    numPoints: 128, // Must be power of 2 for radix-2 FFT
  };

  const potential: PotentialFunction = (x: number) => 0.5 * springConstant * x * x;

  // Analytical
  const analytical = solveHarmonicOscillator(springConstant, mass, numStates, gridConfig);

  // DVR
  const dvr = solveDVR(potential, mass, numStates, gridConfig);

  // Matrix Numerov
  const matrixNumerov = solveMatrixNumerov(potential, mass, numStates, gridConfig);

  // Spectral
  const spectral = solveSpectral(potential, mass, numStates, gridConfig);

  // FGH (requires power-of-2 grid)
  const analyticalFGH = solveHarmonicOscillator(springConstant, mass, numStates, gridConfigFGH);
  const fgh = solveFGH(potential, mass, numStates, gridConfigFGH);

  console.log("Ground state energy (E_0):");
  console.log(`  Analytical:     ${Schrodinger1DSolver.joulesToEV(analytical.energies[0]).toFixed(6)} eV`);
  console.log(`  DVR:            ${Schrodinger1DSolver.joulesToEV(dvr.energies[0]).toFixed(6)} eV (${percentageError(dvr.energies[0], analytical.energies[0]).toFixed(4)}% error)`);
  console.log(`  Matrix Numerov: ${Schrodinger1DSolver.joulesToEV(matrixNumerov.energies[0]).toFixed(6)} eV (${percentageError(matrixNumerov.energies[0], analytical.energies[0]).toFixed(4)}% error)`);
  console.log(`  Spectral:       ${Schrodinger1DSolver.joulesToEV(spectral.energies[0]).toFixed(6)} eV (${percentageError(spectral.energies[0], analytical.energies[0]).toFixed(4)}% error)`);
  console.log(`  FGH:            ${Schrodinger1DSolver.joulesToEV(fgh.energies[0]).toFixed(6)} eV (${percentageError(fgh.energies[0], analyticalFGH.energies[0]).toFixed(4)}% error)`);

  const dvrPassed = percentageError(dvr.energies[0], analytical.energies[0]) < 1.0;
  const matrixNumerovPassed = percentageError(matrixNumerov.energies[0], analytical.energies[0]) < 1.0;
  const spectralPassed = percentageError(spectral.energies[0], analytical.energies[0]) < 1.0;
  const fghPassed = percentageError(fgh.energies[0], analyticalFGH.energies[0]) < 1.0;

  if (dvrPassed && matrixNumerovPassed && spectralPassed && fghPassed) {
    console.log("\n✓ Quick check passed! All methods are within 1% accuracy.");
  } else {
    console.log("\n✗ Quick check failed! One or more methods exceed 1% error.");
  }

  console.log("============================\n");
}

qppw.register("AccuracyTests", { runAccuracyTests, runQuickAccuracyCheck });

export default { runAccuracyTests, runQuickAccuracyCheck };
