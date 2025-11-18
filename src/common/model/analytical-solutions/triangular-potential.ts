/**
 * Analytical solution for the finite triangular potential well.
 *
 * V(x) = height + offset    for x < 0
 * V(x) = offset             at x = 0
 * V(x) = offset + (height/width) * x    for 0 < x < width
 * V(x) = height + offset    for x > width
 *
 * This creates a triangular well with minimum at x = 0.
 * The solution involves Airy functions Ai(z) and Bi(z).
 *
 * Bound states exist when: offset < E < height + offset
 */

import QuantumConstants from "../QuantumConstants.js";
import { BoundStateResult, GridConfig } from "../PotentialFunction.js";
import { airyAi, airyBi, airyAiPrime, airyBiPrime } from "./math-utilities.js";

/**
 * Analytical solution for the finite triangular potential well.
 *
 * @param height - Height of the potential barrier (Joules)
 * @param width - Width of the triangular region (meters)
 * @param offset - Energy offset/minimum of the well (Joules)
 * @param mass - Particle mass in kg
 * @param numStates - Number of energy levels to calculate
 * @param gridConfig - Grid configuration for wavefunction evaluation
 * @returns Bound state results with energies and wavefunctions
 */
export function solveTriangularPotential(
  height: number,
  width: number,
  offset: number,
  mass: number,
  numStates: number,
  gridConfig: GridConfig,
): BoundStateResult {
  const { HBAR } = QuantumConstants;

  // Slope of the linear potential region: F = height / width
  const F = height / width;

  // Barrier height (maximum potential)
  const V0 = height + offset;

  // Scaling parameter for Airy functions: α = (2mF/ℏ²)^(1/3)
  const alpha = Math.pow((2 * mass * F) / (HBAR * HBAR), 1 / 3);

  // Find bound state energies by solving transcendental equation
  const energies: number[] = [];

  // Energy must be between offset (minimum) and V0 (barrier height)
  // Search for eigenvalues in this range
  const Emin = offset + 1e-6 * Math.abs(height); // Just above the well bottom
  const Emax = V0 - 1e-6 * Math.abs(height); // Just below the barrier

  // Number of search points for initial bracketing
  const numSearchPoints = 1000;
  const dE = (Emax - Emin) / numSearchPoints;

  // Function to evaluate the transcendental equation
  // Returns zero when E is an eigenvalue
  const transcendentalFunction = (E: number): number => {
    // Decay constant in barrier regions
    const kappa = Math.sqrt(2 * mass * (V0 - E)) / HBAR;

    // Classical turning point in linear region
    const x0 = (E - offset) / F;

    // Airy function arguments at boundaries
    const z0 = alpha * (0 - x0); // at x = 0
    const zW = alpha * (width - x0); // at x = width

    // Evaluate Airy functions
    const Ai0 = airyAi(z0);
    const Bi0 = airyBi(z0);
    const AiW = airyAi(zW);
    const BiW = airyBi(zW);

    // Evaluate Airy function derivatives
    const AiPrime0 = airyAiPrime(z0);
    const BiPrime0 = airyBiPrime(z0);
    const AiPrimeW = airyAiPrime(zW);
    const BiPrimeW = airyBiPrime(zW);

    // Left boundary matching coefficients
    const leftAi = kappa * Ai0 - alpha * AiPrime0;
    const leftBi = kappa * Bi0 - alpha * BiPrime0;

    // Right boundary matching coefficients
    const rightAi = kappa * AiW + alpha * AiPrimeW;
    const rightBi = kappa * BiW + alpha * BiPrimeW;

    // Transcendental equation: leftBi * rightAi - leftAi * rightBi = 0
    return leftBi * rightAi - leftAi * rightBi;
  };

  // Find sign changes to bracket roots
  let prevF = transcendentalFunction(Emin);
  const brackets: [number, number][] = [];

  for (let i = 1; i <= numSearchPoints; i++) {
    const E = Emin + i * dE;
    const currF = transcendentalFunction(E);

    if (prevF * currF < 0) {
      // Sign change detected - bracket found
      brackets.push([E - dE, E]);
    }

    prevF = currF;

    // Stop if we have enough brackets
    if (brackets.length >= numStates * 2) {
      break;
    }
  }

  // Refine each bracket using bisection
  for (const [Ea, Eb] of brackets) {
    const energy = refineBisection(transcendentalFunction, Ea, Eb, 1e-12, 100);
    if (energy !== null) {
      energies.push(energy);

      if (energies.length >= numStates) {
        break;
      }
    }
  }

  // Sort energies (should already be sorted, but ensure it)
  energies.sort((a, b) => a - b);

  const actualNumStates = energies.length;

  if (actualNumStates === 0) {
    // No bound states found
    return {
      energies: [],
      wavefunctions: [],
      xGrid: [],
      method: "analytical",
    };
  }

  // Generate grid
  const numPoints = gridConfig.numPoints;
  const xGrid: number[] = [];
  const dx = (gridConfig.xMax - gridConfig.xMin) / (numPoints - 1);
  for (let i = 0; i < numPoints; i++) {
    xGrid.push(gridConfig.xMin + i * dx);
  }

  // Calculate wavefunctions
  const wavefunctions: number[][] = [];

  for (let n = 0; n < actualNumStates; n++) {
    const E = energies[n];
    const wavefunction: number[] = [];

    // Decay constant
    const kappa = Math.sqrt(2 * mass * (V0 - E)) / HBAR;

    // Classical turning point
    const x0 = (E - offset) / F;

    // Airy function arguments at boundaries
    const z0 = alpha * (0 - x0);
    const zW = alpha * (width - x0);

    // Get Airy functions at left boundary to determine A/B ratio
    const Ai0 = airyAi(z0);
    const Bi0 = airyBi(z0);
    const AiPrime0 = airyAiPrime(z0);
    const BiPrime0 = airyBiPrime(z0);

    // From left boundary matching:
    // A/B = -(kappa * Bi0 - alpha * BiPrime0) / (kappa * Ai0 - alpha * AiPrime0)
    const leftAi = kappa * Ai0 - alpha * AiPrime0;
    const leftBi = kappa * Bi0 - alpha * BiPrime0;

    // Choose normalization: set B = 1, A = -leftBi / leftAi
    let A: number, B: number;
    if (Math.abs(leftAi) > 1e-10) {
      B = 1;
      A = -leftBi / leftAi;
    } else {
      // If leftAi is too small, use right boundary instead
      const AiW = airyAi(zW);
      const BiW = airyBi(zW);
      const AiPrimeW = airyAiPrime(zW);
      const BiPrimeW = airyBiPrime(zW);
      const rightAi = kappa * AiW + alpha * AiPrimeW;
      const rightBi = kappa * BiW + alpha * BiPrimeW;

      B = 1;
      A = -rightBi / rightAi;
    }

    // Compute unnormalized wavefunction
    const psiRaw: number[] = [];

    for (const x of xGrid) {
      let psi: number;

      if (x < 0) {
        // Region I: exponential decay to the left
        // ψ = C * exp(κx)
        // C is determined by continuity at x = 0
        const psiAt0 = A * Ai0 + B * Bi0;
        psi = psiAt0 * Math.exp(kappa * x);
      } else if (x <= width) {
        // Region II: linear potential - Airy functions
        // ψ = A * Ai(z) + B * Bi(z)
        const z = alpha * (x - x0);
        psi = A * airyAi(z) + B * airyBi(z);
      } else {
        // Region III: exponential decay to the right
        // ψ = D * exp(-κ(x - width))
        // D is determined by continuity at x = width
        const AiW = airyAi(zW);
        const BiW = airyBi(zW);
        const psiAtW = A * AiW + B * BiW;
        psi = psiAtW * Math.exp(-kappa * (x - width));
      }

      psiRaw.push(psi);
    }

    // Normalize wavefunction
    let normSq = 0;
    for (const psi of psiRaw) {
      normSq += psi * psi * dx;
    }
    const norm = 1 / Math.sqrt(normSq);

    for (const psi of psiRaw) {
      wavefunction.push(norm * psi);
    }

    // Ensure consistent sign convention (positive at first antinode)
    // Find the maximum absolute value and check its sign
    let maxAbsIndex = 0;
    let maxAbsValue = 0;
    for (let i = 0; i < wavefunction.length; i++) {
      if (Math.abs(wavefunction[i]) > maxAbsValue) {
        maxAbsValue = Math.abs(wavefunction[i]);
        maxAbsIndex = i;
      }
    }

    // For ground state (n=0), ensure it's positive
    // For excited states, use alternating convention based on state index
    const shouldBePositive = n % 2 === 0;
    if ((wavefunction[maxAbsIndex] > 0) !== shouldBePositive) {
      for (let i = 0; i < wavefunction.length; i++) {
        wavefunction[i] = -wavefunction[i];
      }
    }

    wavefunctions.push(wavefunction);
  }

  return {
    energies,
    wavefunctions,
    xGrid,
    method: "analytical",
  };
}

/**
 * Refine a root using bisection method.
 */
function refineBisection(
  f: (x: number) => number,
  a: number,
  b: number,
  tolerance: number,
  maxIterations: number
): number | null {
  let fa = f(a);
  let fb = f(b);

  if (fa * fb > 0) {
    return null; // No sign change
  }

  for (let iter = 0; iter < maxIterations; iter++) {
    const c = (a + b) / 2;
    const fc = f(c);

    if (Math.abs(fc) < tolerance || (b - a) / 2 < tolerance) {
      return c;
    }

    if (fa * fc < 0) {
      b = c;
      fb = fc;
    } else {
      a = c;
      fa = fc;
    }
  }

  return (a + b) / 2;
}
