/**
 * Verification script for the improved spectral method implementation
 *
 * This tests the spectral method after removing the empirical correction factor.
 * The method should now work correctly for both harmonic oscillator and infinite square well.
 */

import { solveSpectral } from './src/common/model/SpectralSolver.js';
import { solveHarmonicOscillator } from './src/common/model/analytical-solutions/harmonic-oscillator.js';
import { solveInfiniteWell } from './src/common/model/analytical-solutions/infinite-square-well.js';
import QuantumConstants from './src/common/model/QuantumConstants.js';

function percentageError(numerical: number, analytical: number): number {
  return Math.abs((numerical - analytical) / analytical) * 100;
}

function joulesToEV(joules: number): number {
  return joules / 1.602176634e-19;
}

console.log('========================================');
console.log('Spectral Method Verification');
console.log('After removing empirical correction factor');
console.log('========================================\n');

// Test 1: Harmonic Oscillator
console.log('=== Test 1: Harmonic Oscillator ===');
const mass = QuantumConstants.ELECTRON_MASS;
const omega = 1e15;
const springConstant = mass * omega * omega;
const numStates = 5;

const gridConfig = {
  xMin: -5e-9,
  xMax: 5e-9,
  numPoints: 200,
};

const harmonicPotential = (x: number) => 0.5 * springConstant * x * x;

const spectralHO = solveSpectral(harmonicPotential, mass, numStates, gridConfig);
const analyticalHO = solveHarmonicOscillator(springConstant, mass, numStates, gridConfig);

console.log('Energy levels:');
let hoAllPassed = true;
let hoMaxError = 0;

for (let n = 0; n < numStates; n++) {
  const E_numerical = spectralHO.energies[n];
  const E_analytical = analyticalHO.energies[n];
  const error = percentageError(E_numerical, E_analytical);
  const passed = error < 1.0;
  hoAllPassed = hoAllPassed && passed;
  hoMaxError = Math.max(hoMaxError, error);

  console.log(`  E_${n}: ${passed ? '✓' : '✗'} - Numerical: ${joulesToEV(E_numerical).toFixed(6)} eV, ` +
              `Analytical: ${joulesToEV(E_analytical).toFixed(6)} eV, Error: ${error.toFixed(4)}%`);
}

console.log(`Status: ${hoAllPassed ? '✓ PASSED' : '✗ FAILED'} (Max error: ${hoMaxError.toFixed(4)}%)\n`);

// Test 2: Infinite Square Well
console.log('=== Test 2: Infinite Square Well ===');
const wellWidth = 1e-9;
const wellGridConfig = {
  xMin: -wellWidth / 2,
  xMax: wellWidth / 2,
  numPoints: 150,
};

const infiniteWellPotential = () => 0;

const spectralISW = solveSpectral(infiniteWellPotential, mass, numStates, wellGridConfig);
const analyticalISW = solveInfiniteWell(wellWidth, mass, numStates, wellGridConfig);

console.log('Energy levels:');
let iswAllPassed = true;
let iswMaxError = 0;

for (let n = 0; n < numStates; n++) {
  const E_numerical = spectralISW.energies[n];
  const E_analytical = analyticalISW.energies[n];
  const error = percentageError(E_numerical, E_analytical);
  const passed = error < 1.0;
  iswAllPassed = iswAllPassed && passed;
  iswMaxError = Math.max(iswMaxError, error);

  console.log(`  E_${n + 1}: ${passed ? '✓' : '✗'} - Numerical: ${joulesToEV(E_numerical).toFixed(6)} eV, ` +
              `Analytical: ${joulesToEV(E_analytical).toFixed(6)} eV, Error: ${error.toFixed(4)}%`);
}

console.log(`Status: ${iswAllPassed ? '✓ PASSED' : '✗ FAILED'} (Max error: ${iswMaxError.toFixed(4)}%)\n`);

// Summary
console.log('========================================');
console.log('Summary');
console.log('========================================');
console.log(`Harmonic Oscillator: ${hoAllPassed ? '✓ PASSED' : '✗ FAILED'}`);
console.log(`Infinite Square Well: ${iswAllPassed ? '✓ PASSED' : '✗ FAILED'}`);

if (hoAllPassed && iswAllPassed) {
  console.log('\n✓ All tests passed! The spectral method is working correctly.');
} else {
  console.log('\n✗ Some tests failed. The spectral method may need further investigation.');
}

console.log('========================================\n');
