/**
 * Test to check if numerical solvers correctly produce odd-parity wavefunctions
 * for the 1D Coulomb potential.
 */

import { solveDVR } from "../src/common/model/DVRSolver.js";
import { solveMatrixNumerov } from "../src/common/model/MatrixNumerovSolver.js";
import { solveFGH } from "../src/common/model/FGHSolver.js";
import { solveCoulomb1DPotential } from "../src/common/model/analytical-solutions/coulomb-1d-potential.js";
import {
  PotentialFunction,
  GridConfig,
} from "../src/common/model/PotentialFunction.js";

// Physical constants (SI units)
const _HBAR = 1.054571817e-34; // Reduced Planck constant (J·s)
const M_E = 9.1093837015e-31; // Electron mass (kg)
const E_CHARGE = 1.602176634e-19; // Elementary charge (C)
const EPSILON_0 = 8.8541878128e-12; // Vacuum permittivity (F/m)
const K_E = 1 / (4 * Math.PI * EPSILON_0); // Coulomb constant

// Hydrogen-like Coulomb strength: α = e²/(4πε₀)
const ALPHA_H = K_E * E_CHARGE * E_CHARGE; // ≈ 2.307e-28 J·m

console.log("=== Testing Numerical Solvers for 1D Coulomb Potential ===\n");

// Create the 1D Coulomb potential function
const coulomb1DPotential: PotentialFunction = (x: number) => {
  const absX = Math.abs(x);
  if (absX < 1e-12) {
    return -ALPHA_H / 1e-12; // Avoid singularity
  }
  return -ALPHA_H / absX;
};

// Grid configuration
const gridConfig: GridConfig = {
  xMin: -30e-10,
  xMax: 30e-10,
  numPoints: 128, // Power of 2 for FGH
};

console.log("Grid configuration:");
console.log(`  xMin = ${(gridConfig.xMin * 1e10).toFixed(1)} Å`);
console.log(`  xMax = ${(gridConfig.xMax * 1e10).toFixed(1)} Å`);
console.log(`  numPoints = ${gridConfig.numPoints}`);
console.log("");

// Get analytical solution
const analytical = solveCoulomb1DPotential(ALPHA_H, M_E, 5, gridConfig);

console.log("Analytical solution (reference):");
console.log("  Energies (eV):");
for (let i = 0; i < analytical.energies.length; i++) {
  const energyEV = analytical.energies[i] / E_CHARGE;
  console.log(`    E_${i} = ${energyEV.toFixed(6)} eV`);
}

// Test numerical solvers
const solvers = [
  { name: "DVR", solver: solveDVR },
  { name: "MatrixNumerov", solver: solveMatrixNumerov },
  { name: "FGH", solver: solveFGH },
];

for (const { name, solver } of solvers) {
  console.log(`\n=== ${name} Solver ===`);

  try {
    const result = solver(coulomb1DPotential, M_E, 5, gridConfig);

    // Check energies
    console.log("Energies (eV):");
    for (let i = 0; i < Math.min(5, result.energies.length); i++) {
      const energyEV = result.energies[i] / E_CHARGE;
      const analyticalEV = analytical.energies[i] / E_CHARGE;
      const error = Math.abs(((energyEV - analyticalEV) / analyticalEV) * 100);
      console.log(
        `  E_${i} = ${energyEV.toFixed(6)} eV (error: ${error.toFixed(2)}%)`,
      );
    }

    // Check parity
    console.log("\nParity check (should be ODD):");
    for (
      let state = 0;
      state < Math.min(3, result.wavefunctions.length);
      state++
    ) {
      const wf = result.wavefunctions[state];
      const midIdx = Math.floor(wf.length / 2);
      const offset = 30;
      const leftVal = wf[midIdx - offset]; // ψ(-x)
      const rightVal = wf[midIdx + offset]; // ψ(x)

      // Check if odd: ψ(-x) = -ψ(x)
      const oddError =
        Math.abs(leftVal + rightVal) /
        Math.max(Math.abs(leftVal), Math.abs(rightVal));
      // Check if even: ψ(-x) = ψ(x)
      const evenError =
        Math.abs(leftVal - rightVal) /
        Math.max(Math.abs(leftVal), Math.abs(rightVal));

      const isOdd = oddError < evenError;
      const parity = isOdd ? "ODD" : "EVEN";
      const correct = isOdd;

      console.log(
        `  n=${state}: ψ(-x) = ${leftVal.toExponential(3)}, ψ(x) = ${rightVal.toExponential(3)}`,
      );
      console.log(
        `        Parity: ${parity} ${correct ? "✓" : "✗ INCORRECT - should be ODD!"}`,
      );
    }
  } catch (error) {
    console.log(`  ERROR: ${error}`);
  }
}

console.log("\n=== End of Test ===");
