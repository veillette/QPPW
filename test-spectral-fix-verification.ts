/**
 * Comprehensive test of the Chebyshev spectral solver fix
 * Verifies:
 * 1. Error is < 1% for various N
 * 2. Error does NOT scale with grid size (should converge)
 * 3. Works for particle in a box (known exact solution)
 */

import { solveSpectral } from "./src/common/model/SpectralSolver.js";
import QuantumConstants from "./src/common/model/QuantumConstants.js";
import type { PotentialFunction, GridConfig } from "./src/common/model/PotentialFunction.js";

const { HBAR, ELECTRON_MASS, EV_TO_JOULES } = QuantumConstants;

console.log("\n=== SPECTRAL SOLVER FIX VERIFICATION ===\n");

// Test 1: Particle in a box
console.log("TEST 1: Particle in a Box");
console.log("=".repeat(70));

const L = 2e-9; // 2 nm box
const mass = ELECTRON_MASS;

// Exact eigenvalues for particle in box of width L: E_n = n²π²ℏ²/(2mL²)
function exactEnergy(n: number): number {
  return ((n * n * Math.PI * Math.PI * HBAR * HBAR) / (2 * mass * L * L));
}

// Infinite square well potential
const potential: PotentialFunction = (x: number) => {
  const halfWidth = L / 2;
  if (x >= -halfWidth && x <= halfWidth) {
    return 0;
  } else {
    return 1e100; // Very large to approximate infinity
  }
};

const gridConfigs: GridConfig[] = [
  { xMin: -L / 2, xMax: L / 2, numPoints: 15 },
  { xMin: -L / 2, xMax: L / 2, numPoints: 21 },
  { xMin: -L / 2, xMax: L / 2, numPoints: 31 },
  { xMin: -L / 2, xMax: L / 2, numPoints: 41 },
  { xMin: -L / 2, xMax: L / 2, numPoints: 51 },
  { xMin: -L / 2, xMax: L / 2, numPoints: 61 },
];

console.log("\nBox width: L = " + (L * 1e9).toFixed(1) + " nm");
console.log("Testing ground state (n=1) energy\n");

console.log("N\tComputed (eV)\tExact (eV)\tError %\t\tStatus");
console.log("-".repeat(70));

let allPass = true;
const errors: number[] = [];

for (const gridConfig of gridConfigs) {
  const result = solveSpectral(potential, mass, 3, gridConfig);

  if (result.energies.length < 1) {
    console.log(`${gridConfig.numPoints}\tFAILED: No bound states found`);
    allPass = false;
    continue;
  }

  const computed = result.energies[0] / EV_TO_JOULES;
  const exact = exactEnergy(1) / EV_TO_JOULES;
  const error = Math.abs((computed - exact) / exact) * 100;
  errors.push(error);

  const status = error < 1.0 ? "✓ PASS" : "✗ FAIL";
  if (error >= 1.0) allPass = false;

  console.log(
    `${gridConfig.numPoints}\t${computed.toFixed(6)}\t${exact.toFixed(6)}\t${error.toFixed(4)}%\t\t${status}`
  );
}

console.log("\n" + "=".repeat(70));

// Check if error decreases with N (convergence)
let converging = true;
for (let i = 1; i < errors.length; i++) {
  if (errors[i] > errors[i - 1] * 1.1) {
    // Allow 10% tolerance for non-monotonic convergence
    converging = false;
    break;
  }
}

console.log("\nRESULTS:");
console.log(`  All errors < 1%: ${allPass ? "✓ YES" : "✗ NO"}`);
console.log(`  Converging with N: ${converging ? "✓ YES" : "✗ NO"}`);
console.log(`  Final error (N=${gridConfigs[gridConfigs.length - 1].numPoints}): ${errors[errors.length - 1].toFixed(4)}%`);

// Test 2: Multiple energy levels
console.log("\n\nTEST 2: Multiple Energy Levels (N=41)");
console.log("=".repeat(70));

const result = solveSpectral(potential, mass, 5, { xMin: -L / 2, xMax: L / 2, numPoints: 41 });

console.log("\nn\tComputed (eV)\tExact (eV)\tError %");
console.log("-".repeat(70));

for (let n = 1; n <= Math.min(5, result.energies.length); n++) {
  const computed = result.energies[n - 1] / EV_TO_JOULES;
  const exact = exactEnergy(n) / EV_TO_JOULES;
  const error = Math.abs((computed - exact) / exact) * 100;

  console.log(`${n}\t${computed.toFixed(6)}\t${exact.toFixed(6)}\t${error.toFixed(4)}%`);
}

// Test 3: Different box size
console.log("\n\nTEST 3: Different Box Size (L = 5 nm, N=41)");
console.log("=".repeat(70));

const L2 = 5e-9;
const potential2: PotentialFunction = (x: number) => {
  const halfWidth = L2 / 2;
  if (x >= -halfWidth && x <= halfWidth) {
    return 0;
  } else {
    return 1e100;
  }
};

const result2 = solveSpectral(potential2, mass, 3, {
  xMin: -L2 / 2,
  xMax: L2 / 2,
  numPoints: 41,
});

function exactEnergy2(n: number): number {
  return ((n * n * Math.PI * Math.PI * HBAR * HBAR) / (2 * mass * L2 * L2));
}

console.log("\nn\tComputed (eV)\tExact (eV)\tError %");
console.log("-".repeat(70));

for (let n = 1; n <= Math.min(3, result2.energies.length); n++) {
  const computed = result2.energies[n - 1] / EV_TO_JOULES;
  const exact = exactEnergy2(n) / EV_TO_JOULES;
  const error = Math.abs((computed - exact) / exact) * 100;

  console.log(`${n}\t${computed.toFixed(6)}\t${exact.toFixed(6)}\t${error.toFixed(4)}%`);
}

console.log("\n" + "=".repeat(70));
console.log("\nOVERALL: " + (allPass && converging ? "✓ ALL TESTS PASSED" : "✗ SOME TESTS FAILED"));
console.log("=".repeat(70) + "\n");
