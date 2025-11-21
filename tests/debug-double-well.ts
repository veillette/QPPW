/**
 * Debug script for double well analytical solver.
 */

import { solveDoubleSquareWellAnalytical } from "../src/common/model/analytical-solutions/double-square-well.js";
import QuantumConstants from "../src/common/model/QuantumConstants.js";

const { HBAR: _HBAR, ELECTRON_MASS, EV_TO_JOULES, NM_TO_M } = QuantumConstants;

console.log("Debug: Double Well Analytical Solver\n");

const wellWidth = 5 * NM_TO_M;
const wellDepth = 0.3 * EV_TO_JOULES;
const wellSeparation = 4 * NM_TO_M;
const mass = 0.067 * ELECTRON_MASS;

const gridConfig = {
  xMin: -15 * NM_TO_M,
  xMax: 15 * NM_TO_M,
  numPoints: 201,  // Odd number ensures x=0 is a grid point
};

const result = solveDoubleSquareWellAnalytical(
  wellWidth,
  wellDepth,
  wellSeparation,
  mass,
  2,
  gridConfig
);

console.log(`Found ${result.energies.length} states\n`);

// Find x=0 index
let zeroIdx = 0;
for (let i = 0; i < result.xGrid.length; i++) {
  if (Math.abs(result.xGrid[i]) < Math.abs(result.xGrid[zeroIdx])) {
    zeroIdx = i;
  }
}

console.log(`Grid center: x[${zeroIdx}] = ${(result.xGrid[zeroIdx]/NM_TO_M).toFixed(3)} nm\n`);

// Check first state (should be even parity)
console.log("State 0 (ground state):");
console.log(`Energy: ${(result.energies[0]/EV_TO_JOULES).toFixed(6)} eV`);
console.log(`ψ(0) = ${result.wavefunctions[0][zeroIdx].toExponential(3)}`);

// Check symmetry around x=0
console.log("\nSymmetry check (comparing ±x):");
for (let i = 1; i <= 5; i++) {
  const leftIdx = zeroIdx - i;
  const rightIdx = zeroIdx + i;
  const x = result.xGrid[rightIdx] / NM_TO_M;
  const psiLeft = result.wavefunctions[0][leftIdx];
  const psiRight = result.wavefunctions[0][rightIdx];
  const diff = psiLeft - psiRight;
  const sum = psiLeft + psiRight;

  console.log(`x = ±${x.toFixed(2)} nm: ψ(-x) = ${psiLeft.toExponential(3)}, ψ(+x) = ${psiRight.toExponential(3)}`);
  console.log(`  Difference: ${diff.toExponential(3)}, Sum: ${sum.toExponential(3)}`);
}

// Check first excited state (should be odd parity)
if (result.energies.length > 1) {
  console.log("\n\nState 1 (first excited state):");
  console.log(`Energy: ${(result.energies[1]/EV_TO_JOULES).toFixed(6)} eV`);
  console.log(`ψ(0) = ${result.wavefunctions[1][zeroIdx].toExponential(3)}`);

  console.log("\nAntisymmetry check (comparing ±x):");
  for (let i = 1; i <= 5; i++) {
    const leftIdx = zeroIdx - i;
    const rightIdx = zeroIdx + i;
    const x = result.xGrid[rightIdx] / NM_TO_M;
    const psiLeft = result.wavefunctions[1][leftIdx];
    const psiRight = result.wavefunctions[1][rightIdx];
    const diff = psiLeft - psiRight;
    const sum = psiLeft + psiRight;

    console.log(`x = ±${x.toFixed(2)} nm: ψ(-x) = ${psiLeft.toExponential(3)}, ψ(+x) = ${psiRight.toExponential(3)}`);
    console.log(`  Difference: ${diff.toExponential(3)}, Sum: ${sum.toExponential(3)}`);
  }
}
