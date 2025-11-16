/**
 * Simple test runner to check spectral method accuracy
 */

// Import the necessary modules
import { solveSpectral } from './src/common/model/SpectralSolver.js';
import { solveHarmonicOscillator } from './src/common/model/analytical-solutions/harmonic-oscillator.js';
import { solveInfiniteWell } from './src/common/model/analytical-solutions/infinite-square-well.js';
import QuantumConstants from './src/common/model/QuantumConstants.js';

/**
 * Calculate percentage error
 */
function percentageError(numerical, analytical) {
  return Math.abs((numerical - analytical) / analytical) * 100;
}

/**
 * Convert Joules to eV
 */
function joulesToEV(joules) {
  return joules / 1.602176634e-19;
}

console.log('========================================');
console.log('Spectral Method Accuracy Test');
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

const harmonicPotential = (x) => 0.5 * springConstant * x * x;

const spectralHO = solveSpectral(harmonicPotential, mass, numStates, gridConfig);
const analyticalHO = solveHarmonicOscillator(springConstant, mass, numStates, gridConfig);

console.log('Energy levels:');
for (let n = 0; n < numStates; n++) {
  const E_numerical = spectralHO.energies[n];
  const E_analytical = analyticalHO.energies[n];
  const error = percentageError(E_numerical, E_analytical);
  const passed = error < 1.0;

  console.log(`  E_${n}: ${passed ? '✓' : '✗'} - Numerical: ${joulesToEV(E_numerical).toFixed(6)} eV, ` +
              `Analytical: ${joulesToEV(E_analytical).toFixed(6)} eV, Error: ${error.toFixed(4)}%`);
}

// Test 2: Infinite Square Well
console.log('\n=== Test 2: Infinite Square Well ===');
const wellWidth = 1e-9;
const wellGridConfig = {
  xMin: -wellWidth / 2,
  xMax: wellWidth / 2,
  numPoints: 150,
};

const infiniteWellPotential = (_x) => 0;

const spectralISW = solveSpectral(infiniteWellPotential, mass, numStates, wellGridConfig);
const analyticalISW = solveInfiniteWell(wellWidth, mass, numStates, wellGridConfig);

console.log('Energy levels:');
for (let n = 0; n < numStates; n++) {
  const E_numerical = spectralISW.energies[n];
  const E_analytical = analyticalISW.energies[n];
  const error = percentageError(E_numerical, E_analytical);
  const passed = error < 1.0;

  console.log(`  E_${n + 1}: ${passed ? '✓' : '✗'} - Numerical: ${joulesToEV(E_numerical).toFixed(6)} eV, ` +
              `Analytical: ${joulesToEV(E_analytical).toFixed(6)} eV, Error: ${error.toFixed(4)}%`);
}

console.log('\n========================================');
