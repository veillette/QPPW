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
import {
  calculateAiryAlpha,
  generateGrid,
  normalizeWavefunction,
  applySignConvention,
  refineBisection,
} from "./airy-utilities.js";

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

  // Scaling parameter for Airy functions
  const alpha = calculateAiryAlpha(mass, F);

  // Find bound state energies by solving transcendental equation
  const energies: number[] = [];

  // Energy must be between offset (minimum) and V0 (barrier height)
  const Emin = offset + 1e-6 * Math.abs(height);
  const Emax = V0 - 1e-6 * Math.abs(height);

  // Number of search points for initial bracketing
  const numSearchPoints = 1000;
  const dE = (Emax - Emin) / numSearchPoints;

  // Function to evaluate the transcendental equation
  const transcendentalFunction = (E: number): number => {
    // Decay constant in barrier regions
    const kappa = Math.sqrt(2 * mass * (V0 - E)) / HBAR;

    // Classical turning point in linear region
    const x0 = (E - offset) / F;

    // Airy function arguments at boundaries
    const z0 = alpha * (0 - x0);
    const zW = alpha * (width - x0);

    // Evaluate Airy functions and derivatives
    const Ai0 = airyAi(z0);
    const Bi0 = airyBi(z0);
    const AiW = airyAi(zW);
    const BiW = airyBi(zW);
    const AiPrime0 = airyAiPrime(z0);
    const BiPrime0 = airyBiPrime(z0);
    const AiPrimeW = airyAiPrime(zW);
    const BiPrimeW = airyBiPrime(zW);

    // Boundary matching coefficients
    const leftAi = kappa * Ai0 - alpha * AiPrime0;
    const leftBi = kappa * Bi0 - alpha * BiPrime0;
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
      brackets.push([E - dE, E]);
    }

    prevF = currF;

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

  // Sort energies
  energies.sort((a, b) => a - b);

  const actualNumStates = energies.length;

  if (actualNumStates === 0) {
    return {
      energies: [],
      wavefunctions: [],
      xGrid: [],
      method: "analytical",
    };
  }

  // Generate grid
  const { xGrid, dx } = generateGrid(gridConfig);

  // Calculate wavefunctions
  const wavefunctions: number[][] = [];

  for (let n = 0; n < actualNumStates; n++) {
    const E = energies[n];

    // Decay constant
    const kappa = Math.sqrt(2 * mass * (V0 - E)) / HBAR;

    // Classical turning point
    const x0 = (E - offset) / F;

    // Airy function arguments at boundaries
    const z0 = alpha * (0 - x0);
    const zW = alpha * (width - x0);

    // Get Airy functions at left boundary
    const Ai0 = airyAi(z0);
    const Bi0 = airyBi(z0);
    const AiPrime0 = airyAiPrime(z0);
    const BiPrime0 = airyBiPrime(z0);

    // Boundary matching coefficients
    const leftAi = kappa * Ai0 - alpha * AiPrime0;
    const leftBi = kappa * Bi0 - alpha * BiPrime0;

    // Get Airy functions at right boundary
    const AiW = airyAi(zW);
    const BiW = airyBi(zW);
    const AiPrimeW = airyAiPrime(zW);
    const BiPrimeW = airyBiPrime(zW);
    const rightAi = kappa * AiW + alpha * AiPrimeW;
    const rightBi = kappa * BiW + alpha * BiPrimeW;

    // Choose normalization based on numerical stability
    let A: number, B: number;
    if (Math.abs(leftAi) > Math.abs(rightAi) && Math.abs(leftAi) > 1e-10) {
      B = 1;
      A = -leftBi / leftAi;
    } else if (Math.abs(rightAi) > 1e-10) {
      B = 1;
      A = -rightBi / rightAi;
    } else {
      B = 1;
      A = Math.abs(leftAi) > 1e-15 ? -leftBi / leftAi : 0;
    }

    // Compute unnormalized wavefunction
    const psiRaw: number[] = [];

    // Compute wavefunction at x=0 for left region
    const psiAt0 = A * Ai0 + B * Bi0;

    // For the finite triangular well, Bi(z) grows exponentially for z > 0.
    // The classical turning point is at x0 where z = 0.
    // For x > x0 (classically forbidden in linear region), we're in the
    // exponentially decaying part of the Airy solution.
    //
    // Key insight: Only use Airy functions where they're well-behaved (z < some limit).
    // For large positive z, switch to WKB/exponential approximation.

    // Maximum z value where Airy functions are numerically stable
    // Bi(z) ~ exp(2/3 * z^(3/2)) / sqrt(pi * z^(1/2)) for large z
    // Keep z below ~5 to avoid overflow
    const zMax = 5.0;

    for (const x of xGrid) {
      let psi: number;

      if (x < 0) {
        // Region I: exponential decay to the left
        psi = psiAt0 * Math.exp(kappa * x);
      } else if (x >= width) {
        // Region III: exponential decay to the right (barrier region)
        // Use value at classical turning point and decay from there
        const psiAtX0 = A * airyAi(0) + B * airyBi(0);
        psi = psiAtX0 * Math.exp(-kappa * (x - x0));
      } else {
        // Region II: linear potential
        const z = alpha * (x - x0);

        if (z < zMax) {
          // Safe to use Airy functions
          psi = A * airyAi(z) + B * airyBi(z);
        } else {
          // z is too large, use exponential decay from the point where z = zMax
          const xTransition = x0 + zMax / alpha;
          const psiTransition = A * airyAi(zMax) + B * airyBi(zMax);
          psi = psiTransition * Math.exp(-kappa * (x - xTransition));
        }
      }

      psiRaw.push(psi);
    }

    // Normalize and apply sign convention
    const normalizedPsi = normalizeWavefunction(psiRaw, dx);
    const wavefunction = applySignConvention(normalizedPsi, n);
    wavefunctions.push(wavefunction);
  }

  return {
    energies,
    wavefunctions,
    xGrid,
    method: "analytical",
  };
}
