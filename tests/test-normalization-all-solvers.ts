/**
 * Comprehensive test for wave function normalization across all analytical solvers.
 * This test verifies that ∫|ψ|² dx = 1 for all wavefunctions from all analytical solvers.
 */

import QuantumConstants from "../src/common/model/QuantumConstants.js";
import {
  solveInfiniteWell,
  solveFiniteSquareWell,
  solveHarmonicOscillator,
  solveMorsePotential,
  solvePoschlTellerPotential,
  solveRosenMorsePotential,
  solveEckartPotential,
  solveAsymmetricTrianglePotential,
  solveCoulomb1DPotential,
  solveCoulomb3DPotential,
  solveTriangularPotential,
} from "../src/common/model/analytical-solutions/index.js";
import { GridConfig } from "../src/common/model/PotentialFunction.js";

const { HBAR, ELECTRON_MASS } = QuantumConstants;

/**
 * Calculate the normalization integral ∫|ψ|² dx using the trapezoidal rule.
 * Should return 1 for properly normalized wavefunctions.
 */
function calculateNormalization(
  wavefunction: number[],
  xGrid: number[],
): number {
  if (wavefunction.length !== xGrid.length) {
    throw new Error("Wavefunction and grid must have same length");
  }

  let integral = 0;
  for (let i = 1; i < xGrid.length; i++) {
    const dx = xGrid[i] - xGrid[i - 1];
    const psi2_left = wavefunction[i - 1] * wavefunction[i - 1];
    const psi2_right = wavefunction[i] * wavefunction[i];
    integral += 0.5 * (psi2_left + psi2_right) * dx;
  }

  return integral;
}

/**
 * Test normalization for a given solver with specified parameters.
 */
function testNormalization(
  solverName: string,
  solver: () => { energies: number[]; wavefunctions: number[][]; xGrid: number[] },
  tolerance: number = 1e-3,
): void {
  console.log(`\nTesting ${solverName}...`);

  try {
    const result = solver();
    const { wavefunctions, xGrid } = result;

    let allPassed = true;
    const errors: string[] = [];

    for (let n = 0; n < wavefunctions.length; n++) {
      const norm = calculateNormalization(wavefunctions[n], xGrid);
      const error = Math.abs(norm - 1.0);

      if (error > tolerance) {
        allPassed = false;
        errors.push(
          `  State ${n}: ∫|ψ|² dx = ${norm.toFixed(6)} (error: ${error.toExponential(3)})`,
        );
      } else {
        console.log(
          `  State ${n}: ∫|ψ|² dx = ${norm.toFixed(6)} ✓`,
        );
      }
    }

    if (!allPassed) {
      console.log(`\n❌ ${solverName} FAILED - Normalization errors:`);
      errors.forEach((err) => console.log(err));
    } else {
      console.log(`✅ ${solverName} PASSED`);
    }
  } catch (error) {
    console.log(`❌ ${solverName} FAILED with exception:`, error);
  }
}

// Run tests for all analytical solvers
console.log("=== TESTING WAVE FUNCTION NORMALIZATION FOR ALL ANALYTICAL SOLVERS ===\n");

// Common grid configuration with high resolution
const highResGrid: GridConfig = {
  xMin: -10e-9,
  xMax: 10e-9,
  numPoints: 2000,
};

// Test 1: Infinite Square Well
testNormalization("Infinite Square Well", () =>
  solveInfiniteWell(
    2e-9, // well width (2 nm)
    ELECTRON_MASS,
    5, // num states
    highResGrid,
  ),
);

// Test 2: Finite Square Well
testNormalization("Finite Square Well", () =>
  solveFiniteSquareWell(
    2e-9, // well width (2 nm)
    1e-19, // well depth (0.625 eV)
    ELECTRON_MASS,
    5, // num states
    highResGrid,
  ),
);

// Test 3: Harmonic Oscillator
// Use a weak spring constant so the characteristic length x_0 = √(ℏ/(mω)) fits in grid
// With k = 1 N/m, ω = √(k/m) ~ 1.05e15 rad/s
// x_0 = √(ℏ/(mω)) ~ 3.3e-10 m = 0.33 nm
// This should fit well within the ±10 nm grid
testNormalization("Harmonic Oscillator", () =>
  solveHarmonicOscillator(
    1.0, // spring constant (weak, suitable for nanoscale)
    ELECTRON_MASS,
    5, // num states
    highResGrid,
  ),
);

// Test 4: Morse Potential
testNormalization("Morse Potential", () =>
  solveMorsePotential(
    1e-19, // dissociation energy
    1e-10, // well width
    0, // equilibrium position
    ELECTRON_MASS,
    5, // num states
    highResGrid,
  ),
);

// Test 5: Pöschl-Teller Potential
testNormalization("Pöschl-Teller Potential", () =>
  solvePoschlTellerPotential(
    1e-19, // potential depth
    1e-9, // well width
    ELECTRON_MASS,
    5, // num states
    highResGrid,
  ),
);

// Test 6: Rosen-Morse Potential
testNormalization("Rosen-Morse Potential", () =>
  solveRosenMorsePotential(
    1e-19, // barrier height
    1e-9, // well width
    ELECTRON_MASS,
    5, // num states
    highResGrid,
  ),
);

// Test 7: Eckart Potential
testNormalization("Eckart Potential", () =>
  solveEckartPotential(
    1e-19, // barrier height
    1e-9, // well width
    ELECTRON_MASS,
    5, // num states
    highResGrid,
  ),
);

// Test 8: Asymmetric Triangle Potential
testNormalization("Asymmetric Triangle Potential", () =>
  solveAsymmetricTrianglePotential(
    1e-19, // potential height
    1e-9, // left width
    1e-9, // right width
    ELECTRON_MASS,
    5, // num states
    highResGrid,
  ),
);

// Test 9: 1D Coulomb Potential
testNormalization("1D Coulomb Potential", () => {
  const coulombStrength = 2.3e-28; // α in J·m
  return solveCoulomb1DPotential(
    coulombStrength,
    ELECTRON_MASS,
    5, // num states
    {
      xMin: -20e-9,
      xMax: 20e-9,
      numPoints: 3000, // Higher resolution for Coulomb
    },
  );
});

// Test 10: 3D Coulomb Potential (Hydrogen atom)
testNormalization("3D Coulomb Potential (Hydrogen)", () => {
  const coulombStrength = 2.3e-28; // α in J·m
  return solveCoulomb3DPotential(
    coulombStrength,
    ELECTRON_MASS,
    5, // num states
    {
      xMin: 0,
      xMax: 20e-9,
      numPoints: 3000, // Higher resolution
    },
  );
});

// Test 11: Triangular Potential
testNormalization("Triangular Potential", () =>
  solveTriangularPotential(
    1e7, // electric field strength
    ELECTRON_MASS,
    5, // num states
    {
      xMin: 0,
      xMax: 10e-9,
      numPoints: 2000,
    },
  ),
);

console.log("\n=== ALL TESTS COMPLETED ===\n");
