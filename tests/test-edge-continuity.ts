/**
 * Detailed test of wavefunction continuity and edge behavior.
 * Inspects wavefunctions at boundaries between regions.
 */

import { solveDoubleSquareWellAnalytical } from "../src/common/model/analytical-solutions/double-square-well.js";
import QuantumConstants from "../src/common/model/QuantumConstants.js";

const { ELECTRON_MASS, EV_TO_JOULES, NM_TO_M } = QuantumConstants;

console.log("Testing Wavefunction Continuity at Edges");
console.log("=".repeat(70));

// Test parameters
const wellWidth = 5 * NM_TO_M;
const wellDepth = 0.3 * EV_TO_JOULES;
const wellSeparation = 4 * NM_TO_M;
const mass = 0.067 * ELECTRON_MASS;
const numStates = 6;

const Linner = wellSeparation / 2;
const Louter = Linner + wellWidth;

console.log(`Well width: ${wellWidth / NM_TO_M} nm`);
console.log(`Well depth: ${wellDepth / EV_TO_JOULES} eV`);
console.log(`Well separation: ${wellSeparation / NM_TO_M} nm`);
console.log(`L_inner = ${Linner / NM_TO_M} nm`);
console.log(`L_outer = ${Louter / NM_TO_M} nm\n`);

const gridConfig = {
  xMin: -15 * NM_TO_M,
  xMax: 15 * NM_TO_M,
  numPoints: 2000,  // High resolution for detailed edge inspection
};

const result = solveDoubleSquareWellAnalytical(
  wellWidth,
  wellDepth,
  wellSeparation,
  mass,
  numStates,
  gridConfig
);

console.log(`Found ${result.energies.length} bound states\n`);

// Inspect each state's edge behavior
for (let stateIdx = 0; stateIdx < Math.min(4, result.energies.length); stateIdx++) {
  const energyEV = result.energies[stateIdx] / EV_TO_JOULES;
  const wf = result.wavefunctions[stateIdx];
  const xGrid = result.xGrid;

  console.log(`\nState ${stateIdx}: E = ${energyEV.toFixed(6)} eV`);
  console.log("-".repeat(70));

  // Find indices for key boundaries
  const findClosestIndex = (xTarget: number) => {
    let minDist = Infinity;
    let idx = 0;
    for (let i = 0; i < xGrid.length; i++) {
      const dist = Math.abs(xGrid[i] - xTarget);
      if (dist < minDist) {
        minDist = dist;
        idx = i;
      }
    }
    return idx;
  };

  const idxLouter = findClosestIndex(Louter);
  const idxLinner = findClosestIndex(Linner);
  const idx0 = findClosestIndex(0);
  const _idxMinusLinner = findClosestIndex(-Linner);
  const _idxMinusLouter = findClosestIndex(-Louter);

  // Check continuity at x = L_outer (well/outside boundary on right)
  console.log("\n  Right edge (L_outer = " + (Louter/NM_TO_M).toFixed(3) + " nm):");
  console.log("    x (nm)      | ψ(x)        | Change");
  for (let i = Math.max(0, idxLouter - 5); i <= Math.min(xGrid.length - 1, idxLouter + 10); i++) {
    const x = xGrid[i] / NM_TO_M;
    const psi = wf[i];
    const change = i > 0 ? ((wf[i] - wf[i-1]) / (xGrid[i] - xGrid[i-1]) * NM_TO_M) : 0;
    const marker = i === idxLouter ? " <-- L_outer" : "";
    console.log(
      `    ${x.toFixed(6).padStart(11)} | ${psi.toExponential(6)} | ${change.toExponential(2)}${marker}`
    );
  }

  // Check continuity at x = L_inner (barrier/well boundary on right)
  console.log("\n  Right barrier/well boundary (L_inner = " + (Linner/NM_TO_M).toFixed(3) + " nm):");
  console.log("    x (nm)      | ψ(x)        | Change");
  for (let i = Math.max(0, idxLinner - 5); i <= Math.min(xGrid.length - 1, idxLinner + 5); i++) {
    const x = xGrid[i] / NM_TO_M;
    const psi = wf[i];
    const change = i > 0 ? ((wf[i] - wf[i-1]) / (xGrid[i] - xGrid[i-1]) * NM_TO_M) : 0;
    const marker = i === idxLinner ? " <-- L_inner" : "";
    console.log(
      `    ${x.toFixed(6).padStart(11)} | ${psi.toExponential(6)} | ${change.toExponential(2)}${marker}`
    );
  }

  // Check at x = 0 (center of barrier)
  console.log("\n  Center (x = 0):");
  console.log("    x (nm)      | ψ(x)        | Change");
  for (let i = Math.max(0, idx0 - 5); i <= Math.min(xGrid.length - 1, idx0 + 5); i++) {
    const x = xGrid[i] / NM_TO_M;
    const psi = wf[i];
    const change = i > 0 ? ((wf[i] - wf[i-1]) / (xGrid[i] - xGrid[i-1]) * NM_TO_M) : 0;
    const marker = i === idx0 ? " <-- x=0" : "";
    console.log(
      `    ${x.toFixed(6).padStart(11)} | ${psi.toExponential(6)} | ${change.toExponential(2)}${marker}`
    );
  }

  // Check derivative continuity at L_outer
  const dx = xGrid[1] - xGrid[0];
  const derivInside = (wf[idxLouter] - wf[idxLouter - 1]) / dx;
  const derivOutside = (wf[idxLouter + 1] - wf[idxLouter]) / dx;
  const derivJump = Math.abs(derivOutside - derivInside) / Math.max(Math.abs(derivInside), Math.abs(derivOutside));

  console.log("\n  Derivative continuity at L_outer:");
  console.log(`    Inside:  ${derivInside.toExponential(4)}`);
  console.log(`    Outside: ${derivOutside.toExponential(4)}`);
  console.log(`    Relative jump: ${(derivJump * 100).toFixed(2)}%`);

  if (derivJump > 0.1) {
    console.log(`    WARNING: Large derivative discontinuity!`);
  }
}

console.log("\n" + "=".repeat(70));
console.log("Edge continuity inspection complete");
