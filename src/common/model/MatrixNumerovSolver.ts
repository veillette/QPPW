/**
 * Matrix Numerov method for solving the 1D time-independent Schrödinger equation.
 *
 * This method reformulates the Numerov finite difference scheme as a matrix
 * eigenvalue problem, combining the accuracy of the Numerov formula with the
 * robustness of matrix diagonalization.
 *
 * The TISE is: -ℏ²/(2m) d²ψ/dx² + V(x)ψ = Eψ
 *
 * The Numerov formula gives:
 * ψ_{i-1} - 2ψ_i + ψ_{i+1} = (h²/12)[f_{i-1}ψ_{i-1} + 10f_iψ_i + f_{i+1}ψ_{i+1}]
 * where f_i = h²·2m(E-V_i)/ℏ²
 *
 * Rearranging yields a generalized eigenvalue problem: H·ψ = E·S·ψ
 * where H includes kinetic and potential terms, and S is the Numerov overlap matrix.
 *
 * References:
 * - T. G. Walker, "Matrix Numerov method for solving Schrödinger's equation"
 *   https://pages.physics.wisc.edu/~tgwalker/106.Numerov.pdf
 */

import QuantumConstants from "./QuantumConstants.js";
import { BoundStateResult, GridConfig, PotentialFunction } from "./PotentialFunction.js";
import { Matrix, diagonalize, normalizeWavefunction, matrixToArray } from "./LinearAlgebraUtils.js";
import qppw from "../../QPPWNamespace.js";

/**
 * Solve the 1D Schrödinger equation using the matrix Numerov method.
 *
 * This method combines the O(h⁴) accuracy of the Numerov formula with matrix
 * diagonalization, making it more robust than the shooting method for finding
 * all eigenvalues simultaneously.
 *
 * @param potential - Function V(x) that returns potential energy in Joules
 * @param mass - Particle mass in kg
 * @param numStates - Number of bound states to find
 * @param gridConfig - Grid configuration
 * @returns Bound state results
 */
export function solveMatrixNumerov(
  potential: PotentialFunction,
  mass: number,
  numStates: number,
  gridConfig: GridConfig,
): BoundStateResult {
  const { xMin, xMax, numPoints } = gridConfig;
  const dx = (xMax - xMin) / (numPoints - 1);

  // Generate grid
  const xGrid: number[] = [];
  for (let i = 0; i < numPoints; i++) {
    xGrid.push(xMin + i * dx);
  }

  // Evaluate potential on grid
  const V = xGrid.map(potential);

  const N = numPoints;
  const { HBAR } = QuantumConstants;

  // Construct the generalized eigenvalue problem: H·ψ = E·S·ψ
  //
  // From the Numerov formula, we derive:
  // H is the Hamiltonian matrix (kinetic + potential terms)
  // S is the overlap matrix (accounts for Numerov corrections)

  // Physical constants
  const coeff = (2 * mass) / (HBAR * HBAR);
  const h2 = dx * dx;

  // Initialize matrices using dot's Matrix class
  const H = new Matrix(N, N);
  const S = new Matrix(N, N);

  // Construct matrices using Numerov formulation
  // The equation (after rearranging) becomes:
  //
  // [1 + h²/12·2m/ℏ²·V_{i-1}]·ψ_{i-1} - [2 - 10h²/12·2m/ℏ²·V_i]·ψ_i + [1 + h²/12·2m/ℏ²·V_{i+1}]·ψ_{i+1}
  //   = E·h²·2m/ℏ²·[ψ_{i-1}/12 + 10ψ_i/12 + ψ_{i+1}/12]
  //
  // This gives us H·ψ = E·S·ψ where:

  const factor = (h2 / 12) * coeff;

  for (let i = 0; i < N; i++) {
    // Hamiltonian matrix H (kinetic + potential energy)
    // Diagonal term
    H.set(i, i, (2 / h2) - (10 * factor * V[i]));

    // Off-diagonal terms (kinetic energy coupling)
    if (i > 0) {
      const offDiagValue = -(1 / h2) - (factor * V[i - 1]);
      H.set(i, i - 1, offDiagValue);
      H.set(i - 1, i, offDiagValue); // Symmetric
    }

    // Overlap matrix S (Numerov correction terms)
    S.set(i, i, 10 * factor);

    if (i > 0) {
      S.set(i, i - 1, factor);
      S.set(i - 1, i, factor); // Symmetric
    }
  }

  // Solve generalized eigenvalue problem: H·v = λ·S·v
  // Transform to standard form: (S^(-1)·H)·v = λ·v
  // Use dot's Matrix.inverse() and Matrix.times()
  const Sinv = S.inverse();
  const SinvH = Sinv.times(H) as Matrix;

  // Convert to array and diagonalize to get eigenvalues and eigenvectors
  const SinvHArray = matrixToArray(SinvH);
  const eigen = diagonalize(SinvHArray);

  // Sort eigenvalues by energy
  const sortedIndices = eigen.eigenvalues
    .map((e, i) => ({ energy: e, index: i }))
    .sort((a, b) => a.energy - b.energy)
    .map((item) => item.index);

  // Extract bound states
  const energies: number[] = [];
  const wavefunctions: number[][] = [];

  // Estimate boundary potential
  const V_boundary = Math.max(V[0], V[N - 1]);

  for (let i = 0; i < Math.min(numStates, N); i++) {
    const idx = sortedIndices[i];
    const energy = eigen.eigenvalues[idx];

    // Only include bound states (E < V at boundaries)
    // Also filter out negative energies that are too large (numerical artifacts)
    if (energy < V_boundary && isFinite(energy)) {
      energies.push(energy);

      // Extract and normalize wavefunction
      const wavefunction = [...eigen.eigenvectors[idx]];

      // Apply boundary conditions (force ψ=0 at boundaries)
      wavefunction[0] = 0;
      wavefunction[N - 1] = 0;

      const normalizedPsi = normalizeWavefunction(wavefunction, dx);
      wavefunctions.push(normalizedPsi);
    }
  }

  return {
    energies,
    wavefunctions,
    xGrid,
    method: "numerov", // Use same method tag as traditional Numerov
  };
}

qppw.register("MatrixNumerovSolver", { solveMatrixNumerov });
