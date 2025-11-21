/**
 * Diagnostic script to understand the edge discontinuity issue.
 */

import { solveDoubleSquareWellAnalytical } from "../src/common/model/analytical-solutions/double-square-well.js";
import QuantumConstants from "../src/common/model/QuantumConstants.js";

const { HBAR, ELECTRON_MASS, EV_TO_JOULES, NM_TO_M } = QuantumConstants;

// Test parameters
const wellWidth = 5 * NM_TO_M;
const wellDepth = 0.3 * EV_TO_JOULES;
const wellSeparation = 4 * NM_TO_M;
const mass = 0.067 * ELECTRON_MASS;

const Linner = wellSeparation / 2;
const Louter = Linner + wellWidth;
const L = wellWidth;
const V0 = wellDepth;

console.log("Diagnosing Edge Discontinuity");
console.log("=".repeat(70));
console.log(`L_inner = ${Linner / NM_TO_M} nm`);
console.log(`L_outer = ${Louter / NM_TO_M} nm`);
console.log(`L (well width) = ${L / NM_TO_M} nm`);
console.log(`V0 = ${V0 / EV_TO_JOULES} eV\n`);

// Get solution
const gridConfig = {
  xMin: -15 * NM_TO_M,
  xMax: 15 * NM_TO_M,
  numPoints: 2000,
};

const result = solveDoubleSquareWellAnalytical(
  wellWidth,
  wellDepth,
  wellSeparation,
  mass,
  6,
  gridConfig
);

// Check state 2 (the one with the issue)
const stateIdx = 2;
const E = result.energies[stateIdx];
const energyEV = E / EV_TO_JOULES;

console.log(`\nState ${stateIdx}: E = ${energyEV.toFixed(6)} eV`);
console.log("-".repeat(70));

// Calculate wave numbers
const k = Math.sqrt(2 * mass * E) / HBAR;
const kappa = Math.sqrt(2 * mass * (V0 - E)) / HBAR;
const alpha = kappa;

console.log(`k = ${k.toExponential(4)} m^-1`);
console.log(`kappa = ${kappa.toExponential(4)} m^-1`);
console.log(`alpha = ${alpha.toExponential(4)} m^-1`);

// Calculate coefficients (assuming even parity for state 2)
const parity = (stateIdx % 2 === 0) ? "even" : "odd";
console.log(`\nParity: ${parity}`);

const A = 1.0;
let B: number, C: number;

if (parity === "even") {
  B = Math.cosh(kappa * Linner);
  C = (kappa / k) * Math.sinh(kappa * Linner);
} else {
  B = Math.sinh(kappa * Linner);
  C = (kappa / k) * Math.cosh(kappa * Linner);
}

const D = B * Math.cos(k * L) + C * Math.sin(k * L);

console.log(`\nCoefficients:`);
console.log(`A = ${A.toExponential(4)}`);
console.log(`B = ${B.toExponential(4)}`);
console.log(`C = ${C.toExponential(4)}`);
console.log(`D = ${D.toExponential(4)}`);

// Check if D is close to zero (indicating a node at the edge)
if (Math.abs(D) < Math.abs(B) * 0.01) {
  console.log(`\n⚠️  WARNING: D is very small compared to B!`);
  console.log(`   This creates a near-node at the outer edge.`);
}

// Calculate theoretical values at x = L_outer from inside and outside
const psiInside = B * Math.cos(k * L) + C * Math.sin(k * L);
const psiPrimeInside = -B * k * Math.sin(k * L) + C * k * Math.cos(k * L);
const psiOutside = D;
const psiPrimeOutside = -alpha * D;

console.log(`\nAt x = L_outer (theoretical):`);
console.log(`  ψ from inside:  ${psiInside.toExponential(6)}`);
console.log(`  ψ from outside: ${psiOutside.toExponential(6)}`);
console.log(`  ψ' from inside:  ${psiPrimeInside.toExponential(6)}`);
console.log(`  ψ' from outside: ${psiPrimeOutside.toExponential(6)}`);

const psiContinuityError = Math.abs(psiInside - psiOutside);
const psiPrimeContinuityError = Math.abs(psiPrimeInside - psiPrimeOutside) / Math.max(Math.abs(psiPrimeInside), Math.abs(psiPrimeOutside));

console.log(`\n  ψ continuity error: ${psiContinuityError.toExponential(4)}`);
console.log(`  ψ' continuity error: ${(psiPrimeContinuityError * 100).toFixed(2)}%`);

// Now check actual grid values
const xGrid = result.xGrid;
const wf = result.wavefunctions[stateIdx];

// Find points around L_outer
let idxBefore = -1, idxAfter = -1;
for (let i = 0; i < xGrid.length - 1; i++) {
  if (xGrid[i] <= Louter && xGrid[i + 1] > Louter) {
    idxBefore = i;
    idxAfter = i + 1;
    break;
  }
}

if (idxBefore >= 0) {
  console.log(`\nActual grid values around x = L_outer:`);
  console.log(`  x[${idxBefore}] = ${(xGrid[idxBefore] / NM_TO_M).toFixed(6)} nm: ψ = ${wf[idxBefore].toExponential(6)} (inside)`);
  console.log(`  L_outer = ${(Louter / NM_TO_M).toFixed(6)} nm`);
  console.log(`  x[${idxAfter}] = ${(xGrid[idxAfter] / NM_TO_M).toFixed(6)} nm: ψ = ${wf[idxAfter].toExponential(6)} (outside)`);

  const jump = Math.abs(wf[idxAfter] - wf[idxBefore]);
  console.log(`  Jump: ${jump.toExponential(4)}`);

  if (jump > Math.abs(wf[idxBefore]) * 0.1) {
    console.log(`  ⚠️  WARNING: Large discontinuity detected!`);
  }
}

// Check transcendental equation
const numerator = -k * B * Math.sin(k * L) + k * C * Math.cos(k * L);
const denominator = B * Math.cos(k * L) + C * Math.sin(k * L);
const lhs = numerator / denominator;
const rhs = alpha;
const residual = lhs + rhs;

console.log(`\nTranscendental equation check:`);
console.log(`  ψ'/ψ (from matching) = ${lhs.toExponential(6)}`);
console.log(`  -α (expected)        = ${(-rhs).toExponential(6)}`);
console.log(`  Residual (should be ~0): ${residual.toExponential(4)}`);

console.log("\n" + "=".repeat(70));
