/**
 * Check if the grid size is appropriate for the harmonic oscillator
 */

const HBAR = 1.054571817e-34;
const ELECTRON_MASS = 9.1093837015e-31;
const JOULES_TO_EV = 6.241509074e18;

const mass = ELECTRON_MASS;
const omega = 1e15;
const springConstant = mass * omega * omega;

// Characteristic length scale of harmonic oscillator
const x_0 = Math.sqrt(HBAR / (mass * omega));

console.log('\n=== Harmonic Oscillator Length Scale ===\n');
console.log(`Characteristic length x_0 = sqrt(ℏ/mω) = ${(x_0 * 1e9).toFixed(4)} nm`);
console.log(`Wavefunction is significant within ±3·x_0 = ±${(3 * x_0 * 1e9).toFixed(4)} nm`);

// Current grid
const xMin_current = -5e-9;
const xMax_current = 5e-9;
console.log(`\nCurrent grid: ${(xMin_current * 1e9).toFixed(1)} nm to ${(xMax_current * 1e9).toFixed(1)} nm`);
console.log(`Grid extends to ±${(xMax_current / x_0).toFixed(1)} · x_0`);

// Potential at boundaries
const V_boundary = 0.5 * springConstant * xMax_current * xMax_current;
console.log(`Potential at boundaries: ${(V_boundary * JOULES_TO_EV).toFixed(2)} eV`);

// Ground state energy
const E_0 = HBAR * omega * 0.5;
console.log(`Ground state energy: ${(E_0 * JOULES_TO_EV).toFixed(6)} eV`);

console.log(`\nRatio V(boundary)/E_0 = ${(V_boundary / E_0).toFixed(1)}`);

// Suggested grid
const xMin_suggested = -3 * x_0;
const xMax_suggested = 3 * x_0;
console.log(`\nSuggested grid: ${(xMin_suggested * 1e9).toFixed(4)} nm to ${(xMax_suggested * 1e9).toFixed(4)} nm`);

const V_boundary_suggested = 0.5 * springConstant * xMax_suggested * xMax_suggested;
console.log(`Potential at suggested boundaries: ${(V_boundary_suggested * JOULES_TO_EV).toFixed(4)} eV`);
console.log(`Ratio V(suggested)/E_0 = ${(V_boundary_suggested / E_0).toFixed(2)}`);
