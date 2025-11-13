/**
 * Example usage and verification of the 1D Schrödinger equation solver.
 * This file demonstrates how to use the solver for different potential types.
 *
 * To run these examples, import and call the functions in your application:
 * import { runAllExamples } from './common/model/SolverExamples.js';
 * runAllExamples();
 */

import Schrodinger1DSolver, { NumericalMethod } from "./Schrodinger1DSolver.js";
import { PotentialType } from "./PotentialFunction.js";
import QuantumConstants from "./QuantumConstants.js";
import qppw from "../../QPPWNamespace.js";

/**
 * Example 1: Infinite square well using analytical solution
 */
export function exampleInfiniteWellAnalytical(): void {
  console.log("\n=== Example 1: Infinite Square Well (Analytical) ===");

  const solver = new Schrodinger1DSolver();
  const wellWidth = 1e-9; // 1 nm
  const mass = QuantumConstants.ELECTRON_MASS;
  const numStates = 5;

  const gridConfig = {
    xMin: 0,
    xMax: wellWidth,
    numPoints: 100,
  };

  const result = solver.solveAnalyticalIfPossible(
    {
      type: PotentialType.INFINITE_WELL,
      wellWidth: wellWidth,
    },
    mass,
    numStates,
    gridConfig,
  );

  console.log(`Method used: ${result.method}`);
  console.log("\nEnergy eigenvalues (eV):");
  result.energies.forEach((E, n) => {
    const E_eV = Schrodinger1DSolver.joulesToEV(E);
    console.log(`  E_${n + 1} = ${E_eV.toFixed(4)} eV`);
  });

  // Expected values for 1 nm well:
  // E_1 ≈ 0.3760 eV, E_2 ≈ 1.5040 eV, E_3 ≈ 3.3840 eV, etc.
  console.log("\nExpected E_1 ≈ 0.3760 eV (verify above matches)");
}

/**
 * Example 2: Infinite square well using DVR method
 */
export function exampleInfiniteWellDVR(): void {
  console.log("\n=== Example 2: Infinite Square Well (DVR Numerical) ===");

  const solver = new Schrodinger1DSolver(NumericalMethod.DVR);
  const wellWidth = 1e-9; // 1 nm
  const mass = QuantumConstants.ELECTRON_MASS;
  const numStates = 5;

  // Create potential function
  const potential = Schrodinger1DSolver.createInfiniteWellPotential(wellWidth);

  const gridConfig = {
    xMin: 0,
    xMax: wellWidth,
    numPoints: 100,
  };

  const result = solver.solveNumerical(potential, mass, numStates, gridConfig);

  console.log(`Method used: ${result.method}`);
  console.log("\nEnergy eigenvalues (eV):");
  result.energies.forEach((E, n) => {
    const E_eV = Schrodinger1DSolver.joulesToEV(E);
    console.log(`  E_${n + 1} = ${E_eV.toFixed(4)} eV`);
  });

  console.log("\nExpected E_1 ≈ 0.3760 eV (should be close to analytical)");
}

/**
 * Example 3: Harmonic oscillator using analytical solution
 */
export function exampleHarmonicOscillatorAnalytical(): void {
  console.log("\n=== Example 3: Harmonic Oscillator (Analytical) ===");

  const solver = new Schrodinger1DSolver();
  const mass = QuantumConstants.ELECTRON_MASS;
  const omega = 1e15; // Angular frequency (rad/s)
  const springConstant = mass * omega * omega; // k = mω²
  const numStates = 5;

  const gridConfig = {
    xMin: -5e-9,
    xMax: 5e-9,
    numPoints: 200,
  };

  const result = solver.solveAnalyticalIfPossible(
    {
      type: PotentialType.HARMONIC_OSCILLATOR,
      springConstant: springConstant,
    },
    mass,
    numStates,
    gridConfig,
  );

  console.log(`Method used: ${result.method}`);
  console.log("\nEnergy eigenvalues (eV):");
  result.energies.forEach((E, n) => {
    const E_eV = Schrodinger1DSolver.joulesToEV(E);
    console.log(`  E_${n} = ${E_eV.toFixed(6)} eV`);
  });

  // Expected: E_n = ℏω(n + 1/2)
  const E0_expected = (QuantumConstants.HBAR * omega) / 2;
  const E0_eV = Schrodinger1DSolver.joulesToEV(E0_expected);
  console.log(`\nExpected E_0 = ${E0_eV.toFixed(6)} eV (verify above matches)`);
}

/**
 * Example 4: Finite square well using analytical solution
 */
export function exampleFiniteWellAnalytical(): void {
  console.log("\n=== Example 4: Finite Square Well (Analytical) ===");

  const solver = new Schrodinger1DSolver();
  const wellWidth = 1e-9; // 1 nm
  const wellDepth = Schrodinger1DSolver.eVToJoules(10); // 10 eV well depth
  const mass = QuantumConstants.ELECTRON_MASS;
  const numStates = 5;

  const gridConfig = {
    xMin: -2e-9,
    xMax: 2e-9,
    numPoints: 300,
  };

  const result = solver.solveAnalyticalIfPossible(
    {
      type: PotentialType.FINITE_WELL,
      wellWidth: wellWidth,
      wellDepth: wellDepth,
    },
    mass,
    numStates,
    gridConfig,
  );

  console.log(`Method used: ${result.method}`);
  console.log(`Found ${result.energies.length} bound states\n`);
  console.log("Energy eigenvalues (eV) relative to bottom of well:");
  result.energies.forEach((E, n) => {
    const E_eV = Schrodinger1DSolver.joulesToEV(E);
    const E_relative_eV = E_eV + 10; // Shift by well depth to show energy above bottom
    console.log(`  E_${n + 1} = ${E_relative_eV.toFixed(4)} eV (${E_eV.toFixed(4)} eV relative to V=0)`);
  });

  console.log("\nNote: Energies found by solving transcendental equations");
  console.log("States alternate between even and odd parity");
}

/**
 * Example 5: Finite square well using DVR method
 */
export function exampleFiniteWellDVR(): void {
  console.log("\n=== Example 5: Finite Square Well (DVR Numerical) ===");

  const solver = new Schrodinger1DSolver(NumericalMethod.DVR);
  const wellWidth = 1e-9; // 1 nm
  const wellDepth = Schrodinger1DSolver.eVToJoules(5); // 5 eV well depth
  const mass = QuantumConstants.ELECTRON_MASS;
  const numStates = 3;

  // Create finite well potential
  const potential = Schrodinger1DSolver.createFiniteWellPotential(
    wellWidth,
    wellDepth,
    wellWidth / 2, // Center at middle of grid
  );

  const gridConfig = {
    xMin: 0,
    xMax: wellWidth * 2,
    numPoints: 150,
  };

  const result = solver.solveNumerical(potential, mass, numStates, gridConfig);

  console.log(`Method used: ${result.method}`);
  console.log(`Found ${result.energies.length} bound states\n`);
  console.log("Energy eigenvalues (eV):");
  result.energies.forEach((E, n) => {
    const E_eV = Schrodinger1DSolver.joulesToEV(E);
    console.log(`  E_${n + 1} = ${E_eV.toFixed(4)} eV`);
  });
}

/**
 * Example 6: Comparing Numerov and DVR methods
 */
export function exampleCompareMethodsfiniteWell(): void {
  console.log("\n=== Example 6: Compare Numerov vs DVR for Finite Well ===");

  const wellWidth = 1e-9; // 1 nm
  const wellDepth = Schrodinger1DSolver.eVToJoules(5); // 5 eV
  const mass = QuantumConstants.ELECTRON_MASS;
  const numStates = 3;

  const potential = Schrodinger1DSolver.createFiniteWellPotential(
    wellWidth,
    wellDepth,
    wellWidth / 2,
  );

  const gridConfig = {
    xMin: 0,
    xMax: wellWidth * 2,
    numPoints: 150,
  };

  // DVR method
  const solverDVR = new Schrodinger1DSolver(NumericalMethod.DVR);
  const resultDVR = solverDVR.solveNumerical(
    potential,
    mass,
    numStates,
    gridConfig,
  );

  console.log("DVR Method Results:");
  resultDVR.energies.forEach((E, n) => {
    const E_eV = Schrodinger1DSolver.joulesToEV(E);
    console.log(`  E_${n + 1} = ${E_eV.toFixed(6)} eV`);
  });

  // Numerov method
  const solverNumerov = new Schrodinger1DSolver(NumericalMethod.NUMEROV);
  const energyRange: [number, number] = [-wellDepth, 0];
  const resultNumerov = solverNumerov.solveNumerical(
    potential,
    mass,
    numStates,
    gridConfig,
    energyRange,
  );

  console.log("\nNumerov Method Results:");
  resultNumerov.energies.forEach((E, n) => {
    const E_eV = Schrodinger1DSolver.joulesToEV(E);
    console.log(`  E_${n + 1} = ${E_eV.toFixed(6)} eV`);
  });

  console.log("\nBoth methods should give similar results.");
}

/**
 * Example 7: Morse potential using analytical solution
 */
export function exampleMorsePotentialAnalytical(): void {
  console.log("\n=== Example 7: Morse Potential (Analytical) ===");

  const solver = new Schrodinger1DSolver();
  const mass = QuantumConstants.ELECTRON_MASS;

  // Typical values for a diatomic molecule (scaled for electron mass)
  const dissociationEnergy = Schrodinger1DSolver.eVToJoules(5); // 5 eV
  const wellWidth = 1e10; // 1/Å (inverse meters)
  const equilibriumPosition = 0; // Center at origin
  const numStates = 10;

  const gridConfig = {
    xMin: -2e-10,
    xMax: 5e-10,
    numPoints: 200,
  };

  const result = solver.solveAnalyticalIfPossible(
    {
      type: PotentialType.MORSE,
      dissociationEnergy: dissociationEnergy,
      wellWidth: wellWidth,
      equilibriumPosition: equilibriumPosition,
    },
    mass,
    numStates,
    gridConfig,
  );

  console.log(`Method used: ${result.method}`);
  console.log(`Found ${result.energies.length} bound states\n`);
  console.log("Energy eigenvalues (eV) relative to dissociation limit:");
  result.energies.forEach((E, n) => {
    const E_eV = Schrodinger1DSolver.joulesToEV(E);
    console.log(`  E_${n} = ${E_eV.toFixed(6)} eV`);
  });

  console.log("\nNote: Morse potential includes anharmonicity effects");
  console.log("Energy levels get closer together as n increases");
}

/**
 * Example 8: Pöschl-Teller potential using analytical solution
 */
export function examplePoschlTellerPotentialAnalytical(): void {
  console.log("\n=== Example 8: Pöschl-Teller Potential (Analytical) ===");

  const solver = new Schrodinger1DSolver();
  const mass = QuantumConstants.ELECTRON_MASS;

  const potentialDepth = Schrodinger1DSolver.eVToJoules(10); // 10 eV
  const wellWidth = 5e9; // 5/nm (inverse meters)
  const numStates = 5;

  const gridConfig = {
    xMin: -1e-9,
    xMax: 1e-9,
    numPoints: 200,
  };

  const result = solver.solveAnalyticalIfPossible(
    {
      type: PotentialType.POSCHL_TELLER,
      potentialDepth: potentialDepth,
      wellWidth: wellWidth,
    },
    mass,
    numStates,
    gridConfig,
  );

  console.log(`Method used: ${result.method}`);
  console.log(`Found ${result.energies.length} bound states\n`);
  console.log("Energy eigenvalues (eV):");
  result.energies.forEach((E, n) => {
    const E_eV = Schrodinger1DSolver.joulesToEV(E);
    console.log(`  E_${n} = ${E_eV.toFixed(6)} eV`);
  });

  console.log("\nNote: Pöschl-Teller potential is exactly solvable");
  console.log("and useful for modeling quantum wells");
}

/**
 * Run all examples
 */
export function runAllExamples(): void {
  console.log("========================================");
  console.log("1D Schrödinger Equation Solver Examples");
  console.log("========================================");

  exampleInfiniteWellAnalytical();
  exampleInfiniteWellDVR();
  exampleHarmonicOscillatorAnalytical();
  exampleFiniteWellAnalytical();
  exampleFiniteWellDVR();
  exampleCompareMethodsfiniteWell();
  exampleMorsePotentialAnalytical();
  examplePoschlTellerPotentialAnalytical();

  console.log("\n========================================");
  console.log("All examples completed!");
  console.log("========================================\n");
}

qppw.register("SolverExamples", { runAllExamples });

export default { runAllExamples };
