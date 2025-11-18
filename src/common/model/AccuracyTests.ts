/**
 * Comprehensive accuracy validation tests for DVR, Spectral, Matrix Numerov, and FGH methods.
 *
 * This file tests the numerical methods against known analytical solutions across:
 * 1. Harmonic Oscillator
 * 2. Infinite Square Well
 * 3. Finite Square Wells (multiple heights and widths)
 * 4. 3D Coulomb Potential
 * 5. Double Square Wells (multiple configurations)
 *
 * Tests are performed across different grid sizes to ensure robustness.
 */

import { solveDVR } from "./DVRSolver.js";
import { solveSpectral } from "./SpectralSolver.js";
import { solveMatrixNumerov } from "./MatrixNumerovSolver.js";
import { solveFGH } from "./FGHSolver.js";
import { solveHarmonicOscillator } from "./analytical-solutions/harmonic-oscillator.js";
import { solveFiniteSquareWell } from "./analytical-solutions/finite-square-well.js";
import { solveCoulomb3DPotential } from "./analytical-solutions/coulomb-3d-potential.js";
import { solveMorsePotential } from "./analytical-solutions/morse-potential.js";
import { solvePoschlTellerPotential } from "./analytical-solutions/poschl-teller-potential.js";
import QuantumConstants from "./QuantumConstants.js";
import { BoundStateResult, GridConfig, PotentialFunction } from "./PotentialFunction.js";
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
  executionTime: number; // in milliseconds
  details: string[];
}

/**
 * Test configuration for different grid sizes (all powers of 2)
 * Limited to 128 max to keep test runtime reasonable
 */
const GRID_SIZES = [32, 64, 128];

/**
 * Calculate percentage error between numerical and analytical values
 */
function percentageError(numerical: number, analytical: number): number {
  return Math.abs((numerical - analytical) / analytical) * 100;
}

/**
 * Test a numerical method against analytical solution
 */
function testMethod(
  methodName: string,
  solver: (pot: PotentialFunction, mass: number, numStates: number, grid: GridConfig) => BoundStateResult,
  potential: PotentialFunction,
  analyticalSolution: BoundStateResult,
  mass: number,
  numStates: number,
  gridConfig: GridConfig,
  testName: string,
  tolerance: number = 1.0
): TestResult {
  const details: string[] = [];

  try {
    // Measure execution time
    const startTime = performance.now();
    const numericalResult = solver(potential, mass, numStates, gridConfig);
    const endTime = performance.now();
    const executionTime = endTime - startTime;

    let maxError = 0;
    let allPassed = true;

    details.push(`Testing ${numStates} energy levels with ${gridConfig.numPoints} grid points:`);
    details.push(`Execution time: ${executionTime.toFixed(2)} ms`);

    const numToTest = Math.min(numStates, numericalResult.energies.length, analyticalSolution.energies.length);

    for (let n = 0; n < numToTest; n++) {
      const E_numerical = numericalResult.energies[n];
      const E_analytical = analyticalSolution.energies[n];
      const error = percentageError(E_numerical, E_analytical);

      const E_numerical_eV = Schrodinger1DSolver.joulesToEV(E_numerical);
      const E_analytical_eV = Schrodinger1DSolver.joulesToEV(E_analytical);

      const passed = error < tolerance;
      allPassed = allPassed && passed;
      maxError = Math.max(maxError, error);

      const status = passed ? "✓" : "✗";
      details.push(
        `  E_${n}: ${status} Num: ${E_numerical_eV.toFixed(6)} eV, ` +
        `Ana: ${E_analytical_eV.toFixed(6)} eV, Err: ${error.toFixed(4)}%`
      );
    }

    return {
      testName: `${methodName} - ${testName}`,
      method: methodName,
      passed: allPassed,
      maxError,
      executionTime,
      details,
    };
  } catch (error) {
    details.push(`ERROR: ${error}`);
    return {
      testName: `${methodName} - ${testName}`,
      method: methodName,
      passed: false,
      maxError: Infinity,
      executionTime: 0,
      details,
    };
  }
}

/**
 * Test harmonic oscillator across different grid sizes
 */
function testHarmonicOscillatorComprehensive(): TestResult[] {
  const results: TestResult[] = [];
  const mass = QuantumConstants.ELECTRON_MASS;
  const omega = 1e15;
  const springConstant = mass * omega * omega;
  const numStates = 10;

  for (const numPoints of GRID_SIZES) {
    const gridConfig: GridConfig = {
      xMin: -5e-9,
      xMax: 5e-9,
      numPoints,
    };

    const potential: PotentialFunction = (x: number) => 0.5 * springConstant * x * x;
    const analytical = solveHarmonicOscillator(springConstant, mass, numStates, gridConfig);

    // Test all methods - Harmonic oscillator should be very accurate (0.1% tolerance)
    results.push(testMethod("DVR", solveDVR, potential, analytical, mass, numStates, gridConfig, `Harmonic Oscillator (N=${numPoints})`, 0.1));
    results.push(testMethod("Spectral", solveSpectral, potential, analytical, mass, numStates, gridConfig, `Harmonic Oscillator (N=${numPoints})`, 0.1));
    results.push(testMethod("MatrixNumerov", solveMatrixNumerov, potential, analytical, mass, numStates, gridConfig, `Harmonic Oscillator (N=${numPoints})`, 0.1));

    // FGH requires power of 2
    if (isPowerOfTwo(numPoints)) {
      results.push(testMethod("FGH", solveFGH, potential, analytical, mass, numStates, gridConfig, `Harmonic Oscillator (N=${numPoints})`, 0.1));
    }
  }

  return results;
}

/**
 * Test finite square wells with various heights and widths
 */
function testFiniteSquareWellsComprehensive(): TestResult[] {
  const results: TestResult[] = [];
  const mass = QuantumConstants.ELECTRON_MASS;
  const numStates = 10;

  // Different well configurations: [width in nm, depth in eV]
  // Using values within simulation range (width: 0.1-3 nm, depth: 0.1-15 eV)
  const configurations = [
    { width: 1e-9, depth: 5 },    // Default simulation values
    { width: 1e-9, depth: 10 },   // Medium depth well
    { width: 2e-9, depth: 5 },    // Wide shallow well
    { width: 0.5e-9, depth: 15 }, // Narrow deep well
  ];

  for (const config of configurations) {
    const wellWidth = config.width;
    const wellDepth = config.depth * QuantumConstants.EV_TO_JOULES;

    for (const numPoints of [64, 128]) { // Use moderate grid sizes for finite wells (powers of 2)
      const gridConfig: GridConfig = {
        xMin: -2 * wellWidth,
        xMax: 2 * wellWidth,
        numPoints,
      };

      // Create finite square well potential
      const potential: PotentialFunction = (x: number) => {
        const halfWidth = wellWidth / 2;
        return (x >= -halfWidth && x <= halfWidth) ? -wellDepth : 0;
      };

      const analytical = solveFiniteSquareWell(wellWidth, wellDepth, mass, numStates, gridConfig);

      if (analytical.energies.length > 0) {
        const testName = `Finite Well (W=${(wellWidth * 1e9).toFixed(1)}nm, D=${config.depth}eV, N=${numPoints})`;

        // Finite wells should achieve 0.5% accuracy despite discontinuities
        results.push(testMethod("DVR", solveDVR, potential, analytical, mass, numStates, gridConfig, testName, 0.5));
        results.push(testMethod("MatrixNumerov", solveMatrixNumerov, potential, analytical, mass, numStates, gridConfig, testName, 0.5));

        // FGH for power-of-2 grids
        if (isPowerOfTwo(numPoints)) {
          results.push(testMethod("FGH", solveFGH, potential, analytical, mass, numStates, gridConfig, testName, 0.5));
        }
      }
    }
  }

  return results;
}

/**
 * Test 3D Coulomb potential (hydrogen atom) across different grid sizes
 */
function testCoulomb3DComprehensive(): TestResult[] {
  const results: TestResult[] = [];
  const mass = QuantumConstants.ELECTRON_MASS;
  const numStates = 10;

  // Coulomb strength for hydrogen atom: α = e²/(4πε₀)
  const e = 1.602176634e-19; // Elementary charge (C)
  const epsilon0 = 8.8541878128e-12; // Vacuum permittivity (F/m)
  const coulombStrength = (e * e) / (4 * Math.PI * epsilon0);

  for (const numPoints of [64, 128]) {
    // Use appropriate grid for Coulomb potential (r > 0)
    const gridConfig: GridConfig = {
      xMin: 1e-12, // Small positive value to avoid singularity
      xMax: 10e-9,  // 10 nm
      numPoints,
    };

    const potential: PotentialFunction = (r: number) => {
      const r_abs = Math.abs(r);
      return r_abs > 1e-12 ? -coulombStrength / r_abs : -coulombStrength / 1e-12;
    };

    const analytical = solveCoulomb3DPotential(coulombStrength, mass, numStates, gridConfig);
    const testName = `3D Coulomb (Hydrogen, N=${numPoints})`;

    // Coulomb potential has singularity, but should still achieve 1% accuracy
    results.push(testMethod("DVR", solveDVR, potential, analytical, mass, numStates, gridConfig, testName, 1.0));
    results.push(testMethod("MatrixNumerov", solveMatrixNumerov, potential, analytical, mass, numStates, gridConfig, testName, 1.0));

    if (isPowerOfTwo(numPoints)) {
      results.push(testMethod("FGH", solveFGH, potential, analytical, mass, numStates, gridConfig, testName, 1.0));
    }
  }

  return results;
}

/**
 * Test Morse potential across different configurations
 * Using simulation-realistic values (width: 0.1-6 nm, depth: 0.1-15 eV)
 */
function testMorsePotentialComprehensive(): TestResult[] {
  const results: TestResult[] = [];
  const mass = QuantumConstants.ELECTRON_MASS;
  const numStates = 10;

  // Morse potential configurations: [width parameter a in nm, dissociation energy in eV]
  // V(x) = D_e * (1 - exp(-(x - x_e)/a))^2
  const configurations = [
    { width: 0.5e-9, depth: 5 },   // Default-like: narrow well, moderate depth
    { width: 1.0e-9, depth: 10 },  // Wide well, deep potential
    { width: 0.3e-9, depth: 8 },   // Narrow well, deep potential
    { width: 1.5e-9, depth: 3 },   // Wide well, shallow potential
  ];

  for (const config of configurations) {
    const wellWidth = config.width;
    const dissociationEnergy = config.depth * QuantumConstants.EV_TO_JOULES;
    const equilibriumPosition = 0; // Center the potential

    for (const numPoints of [64, 128]) {
      // Grid needs to extend beyond the equilibrium position
      const gridConfig: GridConfig = {
        xMin: -4e-9,
        xMax: 4e-9,
        numPoints,
      };

      // Create Morse potential function
      const potential: PotentialFunction = (x: number) => {
        const expTerm = 1 - Math.exp(-(x - equilibriumPosition) / wellWidth);
        return dissociationEnergy * expTerm * expTerm - dissociationEnergy;
      };

      try {
        const analytical = solveMorsePotential(
          dissociationEnergy,
          wellWidth,
          equilibriumPosition,
          mass,
          numStates,
          gridConfig
        );

        if (analytical.energies.length > 0) {
          const testName = `Morse (a=${(wellWidth * 1e9).toFixed(1)}nm, D=${config.depth}eV, N=${numPoints})`;

          // Morse potential should achieve 0.5% accuracy
          results.push(testMethod("DVR", solveDVR, potential, analytical, mass, numStates, gridConfig, testName, 0.5));

          // FGH for power-of-2 grids
          if (isPowerOfTwo(numPoints)) {
            results.push(testMethod("FGH", solveFGH, potential, analytical, mass, numStates, gridConfig, testName, 0.5));
          }
        }
      } catch (error) {
        // Skip configurations that don't support bound states
        console.log(`Skipping Morse config (a=${(wellWidth * 1e9).toFixed(1)}nm, D=${config.depth}eV): ${error}`);
      }
    }
  }

  return results;
}

/**
 * Test Pöschl-Teller potential across different configurations
 * Using simulation-realistic values (width: 0.1-6 nm, depth: 0.1-15 eV)
 */
function testPoschlTellerComprehensive(): TestResult[] {
  const results: TestResult[] = [];
  const mass = QuantumConstants.ELECTRON_MASS;
  const numStates = 10;

  // Pöschl-Teller potential configurations: [width parameter a in nm, potential depth in eV]
  // V(x) = -V_0 / cosh²(x/a)
  const configurations = [
    { width: 0.5e-9, depth: 5 },   // Default-like: narrow well, moderate depth
    { width: 1.0e-9, depth: 10 },  // Wide well, deep potential
    { width: 0.3e-9, depth: 12 },  // Narrow well, deep potential
    { width: 1.5e-9, depth: 3 },   // Wide well, shallow potential
  ];

  for (const config of configurations) {
    const wellWidth = config.width;
    const potentialDepth = config.depth * QuantumConstants.EV_TO_JOULES;

    for (const numPoints of [64, 128]) {
      const gridConfig: GridConfig = {
        xMin: -4e-9,
        xMax: 4e-9,
        numPoints,
      };

      // Create Pöschl-Teller potential function
      const potential: PotentialFunction = (x: number) => {
        const sechVal = 1.0 / Math.cosh(x / wellWidth);
        return -potentialDepth * sechVal * sechVal;
      };

      try {
        const analytical = solvePoschlTellerPotential(
          potentialDepth,
          wellWidth,
          mass,
          numStates,
          gridConfig
        );

        if (analytical.energies.length > 0) {
          const testName = `Pöschl-Teller (a=${(wellWidth * 1e9).toFixed(1)}nm, V₀=${config.depth}eV, N=${numPoints})`;

          // Pöschl-Teller potential should achieve 0.5% accuracy
          results.push(testMethod("DVR", solveDVR, potential, analytical, mass, numStates, gridConfig, testName, 0.5));

          // FGH for power-of-2 grids
          if (isPowerOfTwo(numPoints)) {
            results.push(testMethod("FGH", solveFGH, potential, analytical, mass, numStates, gridConfig, testName, 0.5));
          }
        }
      } catch (error) {
        // Skip configurations that don't support bound states
        console.log(`Skipping Pöschl-Teller config (a=${(wellWidth * 1e9).toFixed(1)}nm, V₀=${config.depth}eV): ${error}`);
      }
    }
  }

  return results;
}

/**
 * Test double square wells with various configurations
 */
function testDoubleSquareWellsComprehensive(): TestResult[] {
  const results: TestResult[] = [];
  const mass = QuantumConstants.ELECTRON_MASS;
  const numStates = 10;

  // Double well configurations: [well width, barrier width, well depth, barrier height]
  // Using values within simulation range (width: 0.1-3 nm, depth: 0.1-15 eV, separation: 0.05-0.7 nm)
  const configurations = [
    { wellWidth: 1.0e-9, barrierWidth: 0.2e-9, wellDepth: 5, barrierHeight: 3 },   // Default simulation values
    { wellWidth: 0.5e-9, barrierWidth: 0.3e-9, wellDepth: 10, barrierHeight: 5 },  // Narrow wells, medium barrier
    { wellWidth: 1.5e-9, barrierWidth: 0.5e-9, wellDepth: 8, barrierHeight: 4 },   // Wide wells, medium barrier
  ];

  for (const config of configurations) {
    for (const numPoints of [64, 128]) {
      const totalWidth = 2 * config.wellWidth + config.barrierWidth;
      const gridConfig: GridConfig = {
        xMin: -totalWidth,
        xMax: totalWidth,
        numPoints,
      };

      // Create double square well potential
      const wellDepth = config.wellDepth * QuantumConstants.EV_TO_JOULES;
      const barrierHeight = config.barrierHeight * QuantumConstants.EV_TO_JOULES;

      const potential: PotentialFunction = (x: number) => {
        const halfBarrier = config.barrierWidth / 2;
        const wellEdge = halfBarrier + config.wellWidth;

        if (Math.abs(x) < halfBarrier) {
          return -wellDepth + barrierHeight; // Barrier region
        } else if (Math.abs(x) < wellEdge) {
          return -wellDepth; // Well regions
        } else {
          return 0; // Outside
        }
      };

      const testName = `Double Well (W=${(config.wellWidth * 1e9).toFixed(1)}nm, B=${(config.barrierWidth * 1e9).toFixed(1)}nm, D=${config.wellDepth}eV, H=${config.barrierHeight}eV, N=${numPoints})`;

      // For double wells, we don't have analytical solutions, so we'll compare methods against each other
      // Use DVR as reference (it's generally reliable)
      try {
        const reference = solveDVR(potential, mass, numStates, gridConfig);

        // Test other methods against DVR (with timing) - methods should agree within 1%
        const numerovStart = performance.now();
        const numerovResult = solveMatrixNumerov(potential, mass, numStates, gridConfig);
        const numerovTime = performance.now() - numerovStart;
        results.push(compareResults("MatrixNumerov", numerovResult, reference, testName, 1.0, numerovTime));

        if (isPowerOfTwo(numPoints)) {
          const fghStart = performance.now();
          const fghResult = solveFGH(potential, mass, numStates, gridConfig);
          const fghTime = performance.now() - fghStart;
          results.push(compareResults("FGH", fghResult, reference, testName, 1.0, fghTime));
        }
      } catch (error) {
        results.push({
          testName: `Double Well - ${testName}`,
          method: "Comparison",
          passed: false,
          maxError: Infinity,
          executionTime: 0,
          details: [`ERROR: ${error}`],
        });
      }
    }
  }

  return results;
}

/**
 * Compare two numerical results (when no analytical solution exists)
 */
function compareResults(
  methodName: string,
  result: BoundStateResult,
  reference: BoundStateResult,
  testName: string,
  tolerance: number,
  executionTime: number
): TestResult {
  const details: string[] = [];
  let maxError = 0;
  let allPassed = true;

  const numToTest = Math.min(result.energies.length, reference.energies.length);
  details.push(`Comparing ${numToTest} energy levels vs DVR reference:`);
  details.push(`Execution time: ${executionTime.toFixed(2)} ms`);

  for (let n = 0; n < numToTest; n++) {
    const E_test = result.energies[n];
    const E_ref = reference.energies[n];
    const error = percentageError(E_test, E_ref);

    const E_test_eV = Schrodinger1DSolver.joulesToEV(E_test);
    const E_ref_eV = Schrodinger1DSolver.joulesToEV(E_ref);

    const passed = error < tolerance;
    allPassed = allPassed && passed;
    maxError = Math.max(maxError, error);

    const status = passed ? "✓" : "✗";
    details.push(
      `  E_${n}: ${status} Method: ${E_test_eV.toFixed(6)} eV, ` +
      `DVR: ${E_ref_eV.toFixed(6)} eV, Err: ${error.toFixed(4)}%`
    );
  }

  return {
    testName: `${methodName} - ${testName}`,
    method: methodName,
    passed: allPassed,
    maxError,
    executionTime,
    details,
  };
}

/**
 * Check if a number is a power of 2
 */
function isPowerOfTwo(n: number): boolean {
  return n > 0 && (n & (n - 1)) === 0;
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
 * Run comprehensive accuracy tests
 */
export function runAccuracyTests(): void {
  console.log("========================================");
  console.log("Comprehensive Numerical Method Tests");
  console.log("========================================");
  console.log("Testing: DVR, Spectral, Matrix Numerov, and FGH");
  console.log("Across multiple potentials and grid sizes");
  console.log("");

  const results: TestResult[] = [];

  // Run all test suites
  console.log("\n--- Testing Harmonic Oscillator ---");
  results.push(...testHarmonicOscillatorComprehensive());

  console.log("\n--- Testing Finite Square Wells ---");
  results.push(...testFiniteSquareWellsComprehensive());

  console.log("\n--- Testing 3D Coulomb Potential ---");
  results.push(...testCoulomb3DComprehensive());

  console.log("\n--- Testing Morse Potential ---");
  results.push(...testMorsePotentialComprehensive());

  console.log("\n--- Testing Pöschl-Teller Potential ---");
  results.push(...testPoschlTellerComprehensive());

  console.log("\n--- Testing Double Square Wells ---");
  results.push(...testDoubleSquareWellsComprehensive());

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

  // Timing statistics by method
  console.log("\n--- Performance Summary ---");
  const methodStats = new Map<string, { total: number; count: number; min: number; max: number }>();

  results.forEach(result => {
    if (!methodStats.has(result.method)) {
      methodStats.set(result.method, { total: 0, count: 0, min: Infinity, max: 0 });
    }
    const stats = methodStats.get(result.method)!;
    stats.total += result.executionTime;
    stats.count++;
    stats.min = Math.min(stats.min, result.executionTime);
    stats.max = Math.max(stats.max, result.executionTime);
  });

  methodStats.forEach((stats, method) => {
    const avg = stats.total / stats.count;
    console.log(`${method}:`);
    console.log(`  Average: ${avg.toFixed(2)} ms | Min: ${stats.min.toFixed(2)} ms | Max: ${stats.max.toFixed(2)} ms | Total: ${stats.total.toFixed(2)} ms`);
  });

  if (failedTests === 0) {
    console.log("\n✓ All tests passed!");
    console.log("All numerical methods produce consistent results across different potentials and grid sizes.");
  } else {
    console.log("\n✗ Some tests failed!");
    console.log("Review the details above for failed tests.");
  }

  console.log("========================================\n");
}

/**
 * Run a quick accuracy check (fewer configurations for faster testing)
 */
export function runQuickAccuracyCheck(): void {
  console.log("\n=== Quick Accuracy Check ===");
  console.log("Testing key configurations...\n");

  const results: TestResult[] = [];
  const mass = QuantumConstants.ELECTRON_MASS;

  // 1. Harmonic oscillator
  const omega = 1e15;
  const springConstant = mass * omega * omega;
  const gridConfig: GridConfig = { xMin: -5e-9, xMax: 5e-9, numPoints: 128 };
  const potential: PotentialFunction = (x: number) => 0.5 * springConstant * x * x;
  const analytical = solveHarmonicOscillator(springConstant, mass, 10, gridConfig);

  results.push(testMethod("DVR", solveDVR, potential, analytical, mass, 10, gridConfig, "Harmonic Oscillator", 0.1));
  results.push(testMethod("MatrixNumerov", solveMatrixNumerov, potential, analytical, mass, 10, gridConfig, "Harmonic Oscillator", 0.1));

  // 2. Finite square well (using simulation-realistic values)
  const wellWidth = 1e-9;
  const wellDepth = 5 * QuantumConstants.EV_TO_JOULES;
  const gridConfig2: GridConfig = { xMin: -2e-9, xMax: 2e-9, numPoints: 128 };
  const potential2: PotentialFunction = (x: number) => {
    const halfWidth = wellWidth / 2;
    return (x >= -halfWidth && x <= halfWidth) ? -wellDepth : 0;
  };
  const analytical2 = solveFiniteSquareWell(wellWidth, wellDepth, mass, 10, gridConfig2);

  if (analytical2.energies.length > 0) {
    results.push(testMethod("DVR", solveDVR, potential2, analytical2, mass, 10, gridConfig2, "Finite Square Well", 0.5));
  }

  // Print results
  results.forEach(printTestResult);

  const passed = results.filter(r => r.passed).length;
  const total = results.length;

  if (passed === total) {
    console.log("\n✓ Quick check passed! All methods working correctly.");
  } else {
    console.log(`\n✗ Quick check: ${passed}/${total} tests passed.`);
  }

  console.log("============================\n");
}

qppw.register("AccuracyTests", { runAccuracyTests, runQuickAccuracyCheck });

export default { runAccuracyTests, runQuickAccuracyCheck };
