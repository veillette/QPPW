#!/usr/bin/env node
/**
 * Comprehensive test suite for QuantumBoundStateSolver on double square well potentials
 *
 * This test systematically varies well parameters to diagnose inconsistencies
 * in the number of bound states found.
 *
 * Usage:
 *   npx tsx --import ./tests/browser-globals.js tests/test-double-well-quantum-bound.ts
 */

import QuantumConstants from '../src/common/model/QuantumConstants.js';
import { QuantumBoundStateSolver } from '../src/common/model/QuantumBoundStateSolver.js';
import { PotentialFunction, GridConfig } from '../src/common/model/PotentialFunction.js';

// Constants
const ELECTRON_MASS = QuantumConstants.ELECTRON_MASS;
const EV_TO_JOULES = QuantumConstants.EV_TO_JOULES;
const NUM_POINTS = 512;
const MAX_STATES_TO_FIND = 30;

// Test result interface
interface TestResult {
  wellWidth: number;      // nm
  wellSeparation: number; // nm
  wellDepth: number;      // eV
  statesFound: number;
  energies: number[];     // in eV
  maxPotential: number;   // in eV
  passed: boolean;
  error?: string;
}

/**
 * Creates a double square well potential function
 */
function createDoubleSquareWellPotential(
  wellWidth: number,      // meters
  wellSeparation: number, // meters (center-to-center)
  wellDepth: number       // Joules
): PotentialFunction {
  const halfWidth = wellWidth / 2;
  const leftCenter = -wellSeparation / 2;
  const rightCenter = wellSeparation / 2;

  return (x: number) => {
    const inLeftWell = Math.abs(x - leftCenter) <= halfWidth;
    const inRightWell = Math.abs(x - rightCenter) <= halfWidth;

    if (inLeftWell || inRightWell) {
      return 0; // Bottom of wells at zero energy
    } else {
      return wellDepth; // Barrier at wellDepth energy
    }
  };
}

/**
 * Run a single test configuration
 */
function runSingleTest(
  wellWidthNm: number,
  wellSeparationNm: number,
  wellDepthEv: number
): TestResult {
  // Convert to SI units
  const wellWidth = wellWidthNm * 1e-9;
  const wellSeparation = wellSeparationNm * 1e-9;
  const wellDepth = wellDepthEv * EV_TO_JOULES;

  // Grid configuration - extend well beyond the wells
  const totalExtent = wellSeparation + wellWidth + 4e-9; // 4nm padding on each side
  const xMin = -totalExtent / 2;
  const xMax = totalExtent / 2;

  const gridConfig: GridConfig = {
    xMin,
    xMax,
    numPoints: NUM_POINTS
  };

  const potential = createDoubleSquareWellPotential(wellWidth, wellSeparation, wellDepth);

  try {
    const solver = new QuantumBoundStateSolver(
      ELECTRON_MASS,
      xMin,
      xMax,
      NUM_POINTS,
      potential
    );

    const states = solver.findMultipleEigenstates(MAX_STATES_TO_FIND);
    const maxPotential = solver.getMaxPotential();

    // Convert energies to eV for display
    const energiesEv = states.map(s => s.energy / EV_TO_JOULES);

    return {
      wellWidth: wellWidthNm,
      wellSeparation: wellSeparationNm,
      wellDepth: wellDepthEv,
      statesFound: states.length,
      energies: energiesEv,
      maxPotential: maxPotential / EV_TO_JOULES,
      passed: true
    };
  } catch (error) {
    return {
      wellWidth: wellWidthNm,
      wellSeparation: wellSeparationNm,
      wellDepth: wellDepthEv,
      statesFound: 0,
      energies: [],
      maxPotential: wellDepthEv,
      passed: false,
      error: error instanceof Error ? error.message : String(error)
    };
  }
}

/**
 * Main test runner
 */
function runTests() {
  console.log('╔════════════════════════════════════════════════════════════════════╗');
  console.log('║   QuantumBoundStateSolver - Double Square Well Test Suite         ║');
  console.log('║   Grid: 512 points                                                 ║');
  console.log('╚════════════════════════════════════════════════════════════════════╝\n');

  // Test parameter ranges
  const wellWidths = [0.5, 1.0, 1.5, 2.0];           // nm
  const wellSeparations = [0.2, 0.5, 1.0, 2.0];      // nm (center-to-center)
  const wellDepths = [3, 5, 8, 10];                   // eV

  const results: TestResult[] = [];
  let totalTests = 0;
  let passedTests = 0;

  console.log('Running systematic parameter sweep...\n');
  console.log('Parameters: wellWidth (nm) | wellSeparation (nm) | wellDepth (eV)\n');
  console.log('═'.repeat(80));

  // Run tests for all combinations
  for (const depth of wellDepths) {
    console.log(`\n┌─ Well Depth: ${depth} eV ${'─'.repeat(60)}`);

    for (const width of wellWidths) {
      for (const separation of wellSeparations) {
        totalTests++;
        const result = runSingleTest(width, separation, depth);
        results.push(result);

        if (result.passed) {
          passedTests++;
        }

        // Format output
        const statusIcon = result.passed ? '✓' : '✗';
        const statesStr = result.statesFound.toString().padStart(2);

        // Check for potential issues
        let warning = '';
        if (result.statesFound < 2) {
          warning = ' ⚠ Too few states!';
        } else if (result.statesFound > 20) {
          warning = ' ⚠ Many states';
        }

        console.log(
          `│ ${statusIcon} w=${width.toFixed(1)}nm, sep=${separation.toFixed(1)}nm: ` +
          `${statesStr} states found` + warning
        );

        // Show energies for problematic cases
        if (result.statesFound < 4 || result.error) {
          if (result.energies.length > 0) {
            console.log(`│   Energies (eV): ${result.energies.map(e => e.toFixed(3)).join(', ')}`);
          }
          if (result.error) {
            console.log(`│   Error: ${result.error}`);
          }
        }
      }
    }
    console.log('└' + '─'.repeat(78));
  }

  // Summary statistics
  console.log('\n' + '═'.repeat(80));
  console.log('SUMMARY');
  console.log('═'.repeat(80));
  console.log(`Total tests: ${totalTests}`);
  console.log(`Passed: ${passedTests}`);
  console.log(`Failed: ${totalTests - passedTests}`);

  // Analyze state count distribution
  const stateCounts = new Map<number, number>();
  for (const result of results) {
    const count = stateCounts.get(result.statesFound) || 0;
    stateCounts.set(result.statesFound, count + 1);
  }

  console.log('\nState count distribution:');
  const sortedCounts = [...stateCounts.entries()].sort((a, b) => a[0] - b[0]);
  for (const [states, count] of sortedCounts) {
    const bar = '█'.repeat(Math.min(count, 40));
    console.log(`  ${states.toString().padStart(2)} states: ${bar} (${count} tests)`);
  }

  // Find problematic configurations
  const fewStates = results.filter(r => r.statesFound < 4);
  if (fewStates.length > 0) {
    console.log('\n⚠ Configurations with < 4 states:');
    for (const r of fewStates) {
      console.log(`  - w=${r.wellWidth}nm, sep=${r.wellSeparation}nm, depth=${r.wellDepth}eV: ${r.statesFound} states`);
      if (r.energies.length > 0) {
        console.log(`    Energies: ${r.energies.map(e => e.toFixed(3)).join(', ')} eV`);
        console.log(`    Max potential: ${r.maxPotential.toFixed(3)} eV`);
      }
    }
  }

  // Check for inconsistencies (similar params, very different state counts)
  console.log('\nChecking for inconsistencies...');
  let inconsistencies = 0;
  for (let i = 0; i < results.length; i++) {
    for (let j = i + 1; j < results.length; j++) {
      const r1 = results[i];
      const r2 = results[j];

      // Check if parameters are similar but state counts differ significantly
      const widthSimilar = Math.abs(r1.wellWidth - r2.wellWidth) < 0.3;
      const sepSimilar = Math.abs(r1.wellSeparation - r2.wellSeparation) < 0.3;
      const depthSimilar = Math.abs(r1.wellDepth - r2.wellDepth) < 1;

      if (widthSimilar && sepSimilar && depthSimilar) {
        const stateDiff = Math.abs(r1.statesFound - r2.statesFound);
        if (stateDiff > 5) {
          inconsistencies++;
          console.log(`  ! Large difference (${stateDiff} states):`);
          console.log(`    Config 1: w=${r1.wellWidth}, sep=${r1.wellSeparation}, d=${r1.wellDepth} -> ${r1.statesFound} states`);
          console.log(`    Config 2: w=${r2.wellWidth}, sep=${r2.wellSeparation}, d=${r2.wellDepth} -> ${r2.statesFound} states`);
        }
      }
    }
  }

  if (inconsistencies === 0) {
    console.log('  No major inconsistencies found.');
  }

  // Detailed analysis section
  console.log('\n' + '═'.repeat(80));
  console.log('DETAILED ENERGY ANALYSIS (Selected configurations)');
  console.log('═'.repeat(80));

  // Show detailed results for a few representative configurations
  const selectedConfigs = [
    { w: 1.0, s: 0.5, d: 5 },
    { w: 1.0, s: 1.0, d: 5 },
    { w: 1.5, s: 0.5, d: 8 },
    { w: 2.0, s: 1.0, d: 10 }
  ];

  for (const config of selectedConfigs) {
    const result = results.find(
      r => r.wellWidth === config.w &&
           r.wellSeparation === config.s &&
           r.wellDepth === config.d
    );

    if (result) {
      console.log(`\nConfiguration: w=${config.w}nm, sep=${config.s}nm, depth=${config.d}eV`);
      console.log(`  States found: ${result.statesFound}`);
      console.log(`  Max potential: ${result.maxPotential.toFixed(4)} eV`);
      console.log(`  Energies (eV):`);

      for (let i = 0; i < result.energies.length; i++) {
        const E = result.energies[i];
        const belowBarrier = E < result.maxPotential;
        const marker = belowBarrier ? '  ' : '! ';
        console.log(`    ${marker}E_${i}: ${E.toFixed(6)} eV ${belowBarrier ? '' : '(above barrier!)'}`);
      }
    }
  }

  // Parameter sensitivity tests
  console.log('\n' + '═'.repeat(80));
  console.log('PARAMETER SENSITIVITY TESTS');
  console.log('═'.repeat(80));
  console.log('\nTesting that small parameter changes produce small, consistent energy changes...\n');

  const sensitivityResults: { test: string; passed: boolean; details: string }[] = [];

  // Test configurations for sensitivity
  const baseConfigs = [
    { w: 1.0, s: 0.5, d: 5 },
    { w: 1.5, s: 1.0, d: 8 },
    { w: 2.0, s: 1.0, d: 10 }
  ];

  const perturbation = 0.05; // 5% change

  for (const base of baseConfigs) {
    const baseResult = runSingleTest(base.w, base.s, base.d);
    if (baseResult.statesFound < 3) continue;

    // Test width perturbation
    const widerResult = runSingleTest(base.w * (1 + perturbation), base.s, base.d);
    const narrowerResult = runSingleTest(base.w * (1 - perturbation), base.s, base.d);

    // Test depth perturbation
    const deeperResult = runSingleTest(base.w, base.s, base.d * (1 + perturbation));
    const shallowerResult = runSingleTest(base.w, base.s, base.d * (1 - perturbation));

    // Test separation perturbation
    const widerSepResult = runSingleTest(base.w, base.s * (1 + perturbation), base.d);
    const narrowerSepResult = runSingleTest(base.w, base.s * (1 - perturbation), base.d);

    // Check width sensitivity: wider well → lower energies
    if (widerResult.statesFound >= 3 && narrowerResult.statesFound >= 3) {
      let widthConsistent = true;
      const widthChanges: string[] = [];
      for (let i = 0; i < Math.min(3, baseResult.energies.length, widerResult.energies.length, narrowerResult.energies.length); i++) {
        const widerLower = widerResult.energies[i] < baseResult.energies[i];
        const narrowerHigher = narrowerResult.energies[i] > baseResult.energies[i];
        if (!widerLower || !narrowerHigher) widthConsistent = false;
        const changeWider = ((widerResult.energies[i] - baseResult.energies[i]) / baseResult.energies[i] * 100).toFixed(2);
        const changeNarrower = ((narrowerResult.energies[i] - baseResult.energies[i]) / baseResult.energies[i] * 100).toFixed(2);
        widthChanges.push(`E_${i}: ${changeWider}%/${changeNarrower}%`);
      }
      sensitivityResults.push({
        test: `Width (w=${base.w}nm, s=${base.s}nm, d=${base.d}eV)`,
        passed: widthConsistent,
        details: widthConsistent ? 'Wider→lower, Narrower→higher ✓' : `Changes: ${widthChanges.join(', ')}`
      });
    }

    // Check depth sensitivity: deeper barrier → higher energies (approaching infinite well limit)
    // In finite wells, energies measured from well bottom increase as barrier height increases
    if (deeperResult.statesFound >= 3 && shallowerResult.statesFound >= 3) {
      let depthConsistent = true;
      const depthChanges: string[] = [];
      for (let i = 0; i < Math.min(3, baseResult.energies.length, deeperResult.energies.length, shallowerResult.energies.length); i++) {
        const deeperHigher = deeperResult.energies[i] > baseResult.energies[i];
        const shallowerLower = shallowerResult.energies[i] < baseResult.energies[i];
        if (!deeperHigher || !shallowerLower) depthConsistent = false;
        const changeDeeper = ((deeperResult.energies[i] - baseResult.energies[i]) / baseResult.energies[i] * 100).toFixed(2);
        const changeShallower = ((shallowerResult.energies[i] - baseResult.energies[i]) / baseResult.energies[i] * 100).toFixed(2);
        depthChanges.push(`E_${i}: ${changeDeeper}%/${changeShallower}%`);
      }
      sensitivityResults.push({
        test: `Depth (w=${base.w}nm, s=${base.s}nm, d=${base.d}eV)`,
        passed: depthConsistent,
        details: depthConsistent ? 'Higher barrier→higher E, Lower barrier→lower E ✓' : `Changes: ${depthChanges.join(', ')}`
      });
    }

    // Check separation sensitivity (ground state: larger separation → lower energy due to tunneling)
    if (widerSepResult.statesFound >= 1 && narrowerSepResult.statesFound >= 1) {
      const sepChanges: string[] = [];
      for (let i = 0; i < Math.min(2, baseResult.energies.length, widerSepResult.energies.length, narrowerSepResult.energies.length); i++) {
        const changeWider = ((widerSepResult.energies[i] - baseResult.energies[i]) / baseResult.energies[i] * 100).toFixed(2);
        const changeNarrower = ((narrowerSepResult.energies[i] - baseResult.energies[i]) / baseResult.energies[i] * 100).toFixed(2);
        sepChanges.push(`E_${i}: ${changeWider}%/${changeNarrower}%`);
      }
      // Check that changes are small (< 20% for 5% parameter change)
      const changesSmall = sepChanges.every(c => {
        const match = c.match(/(-?\d+\.?\d*)%\/(-?\d+\.?\d*)%/);
        if (!match) return false;
        return Math.abs(parseFloat(match[1])) < 20 && Math.abs(parseFloat(match[2])) < 20;
      });
      sensitivityResults.push({
        test: `Separation (w=${base.w}nm, s=${base.s}nm, d=${base.d}eV)`,
        passed: changesSmall,
        details: changesSmall ? `Small changes: ${sepChanges.join(', ')} ✓` : `Large changes: ${sepChanges.join(', ')}`
      });
    }
  }

  // Print sensitivity results
  let sensitivityPassed = 0;
  for (const result of sensitivityResults) {
    const icon = result.passed ? '✓' : '✗';
    console.log(`${icon} ${result.test}`);
    console.log(`  ${result.details}`);
    if (result.passed) sensitivityPassed++;
  }

  console.log(`\nSensitivity tests: ${sensitivityPassed}/${sensitivityResults.length} passed`);

  console.log('\n' + '═'.repeat(80));
  console.log('Test suite complete.');
  console.log('═'.repeat(80));

  // Exit with appropriate code
  process.exit(passedTests === totalTests ? 0 : 1);
}

// Run tests
runTests();
