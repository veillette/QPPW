/**
 * Diagnostic script to understand DVR grid and matrix construction
 */

const HBAR = 1.054571817e-34;
const ELECTRON_MASS = 9.1093837015e-31;
const JOULES_TO_EV = 6.241509074e18;

console.log('\n=== DVR Grid Analysis ===\n');

// Test parameters
const xMin = -5e-9;
const xMax = 5e-9;
const numPoints = 10; // Small number for clarity

// Current implementation (includes boundaries)
const dx_current = (xMax - xMin) / (numPoints - 1);
console.log('Current implementation:');
console.log(`  dx = (xMax - xMin) / (numPoints - 1) = ${dx_current.toExponential(4)}`);
console.log(`  Grid points (first 5):`);
for (let i = 0; i < Math.min(5, numPoints); i++) {
  const x = xMin + i * dx_current;
  console.log(`    x[${i}] = ${x.toExponential(4)}`);
}

// Alternative: Exclude boundaries (typical for finite DVR)
const dx_alternative = (xMax - xMin) / (numPoints + 1);
console.log('\nAlternative (excluding boundaries):');
console.log(`  dx = (xMax - xMin) / (numPoints + 1) = ${dx_alternative.toExponential(4)}`);
console.log(`  Grid points (first 5):`);
for (let i = 1; i <= Math.min(5, numPoints); i++) {
  const x = xMin + i * dx_alternative;
  console.log(`    x[${i-1}] = ${x.toExponential(4)}`);
}

// Check kinetic energy matrix element magnitudes
console.log('\n=== Kinetic Energy Matrix Diagonal Element ===\n');

const mass = ELECTRON_MASS;

const prefactor_current = (HBAR * HBAR) / (2 * mass * dx_current * dx_current);
const T_ii_current = prefactor_current * (Math.PI * Math.PI) / 3;

const prefactor_alternative = (HBAR * HBAR) / (2 * mass * dx_alternative * dx_alternative);
const T_ii_alternative = prefactor_alternative * (Math.PI * Math.PI) / 3;

console.log('Current implementation:');
console.log(`  prefactor = ℏ²/(2m·dx²) = ${prefactor_current.toExponential(4)} J`);
console.log(`  T_ii = prefactor · π²/3 = ${T_ii_current.toExponential(4)} J`);
console.log(`  T_ii = ${(T_ii_current * JOULES_TO_EV).toFixed(6)} eV`);

console.log('\nAlternative (excluding boundaries):');
console.log(`  prefactor = ℏ²/(2m·dx²) = ${prefactor_alternative.toExponential(4)} J`);
console.log(`  T_ii = prefactor · π²/3 = ${T_ii_alternative.toExponential(4)} J`);
console.log(`  T_ii = ${(T_ii_alternative * JOULES_TO_EV).toFixed(6)} eV`);

// Expected ground state energy for harmonic oscillator
const omega = 1e15;
const springConstant = mass * omega * omega;
const E_0_analytical = HBAR * omega * 0.5;

console.log('\n=== Expected Ground State Energy ===\n');
console.log(`Harmonic oscillator (ω = ${omega.toExponential(2)} rad/s):`);
console.log(`  E_0 = ℏω/2 = ${(E_0_analytical * JOULES_TO_EV).toFixed(6)} eV`);

// Ratio analysis
console.log('\n=== Ratio Analysis ===\n');
console.log(`Ratio of T_ii (current) to expected E_0:`);
console.log(`  ${(T_ii_current / E_0_analytical).toFixed(2)}`);
console.log(`Ratio of T_ii (alternative) to expected E_0:`);
console.log(`  ${(T_ii_alternative / E_0_analytical).toFixed(2)}`);

console.log('\n=== Grid Point Comparison ===\n');
const ratio = dx_alternative / dx_current;
console.log(`dx_alternative / dx_current = ${ratio.toFixed(4)}`);
console.log(`(dx_alternative / dx_current)² = ${(ratio * ratio).toFixed(4)}`);
console.log(`This means T matrices differ by factor: ${(1/(ratio*ratio)).toFixed(4)}`);
