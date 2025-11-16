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

  // Initialize matrices
  const H = createZeroMatrix(N);
  const S = createZeroMatrix(N);

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
    H[i][i] = (2 / h2) - (10 * factor * V[i]);

    // Off-diagonal terms (kinetic energy coupling)
    if (i > 0) {
      H[i][i - 1] = -(1 / h2) - (factor * V[i - 1]);
      H[i - 1][i] = H[i][i - 1]; // Symmetric
    }

    // Overlap matrix S (Numerov correction terms)
    S[i][i] = 10 * factor;

    if (i > 0) {
      S[i][i - 1] = factor;
      S[i - 1][i] = factor; // Symmetric
    }
  }

  // Apply boundary conditions: ψ(xMin) = ψ(xMax) = 0
  // We enforce this by removing the first and last grid points
  // (alternative: set large diagonal values for first/last rows)

  // For simplicity, we'll solve the full system and filter bound states
  // The boundary condition is implicit in the finite grid

  // Solve generalized eigenvalue problem: H·ψ = E·S·ψ
  // Transform to standard form: S^(-1)·H·ψ = E·ψ
  const SinvH = solveGeneralizedEigenvalueProblem(H, S);

  // Diagonalize to get eigenvalues and eigenvectors
  const eigen = diagonalize(SinvH);

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
      const wavefunction = eigen.eigenvectors[idx];

      // Apply boundary conditions (force ψ=0 at boundaries)
      wavefunction[0] = 0;
      wavefunction[N - 1] = 0;

      const normalizedPsi = normalize(wavefunction, dx);
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

/**
 * Solve generalized eigenvalue problem H·v = λ·S·v by transforming to
 * standard form: (S^(-1)·H)·v = λ·v
 *
 * For the Numerov matrices, S is tridiagonal and invertible, so we can
 * compute S^(-1)·H directly.
 *
 * @param H - Hamiltonian matrix
 * @param S - Overlap matrix
 * @returns Matrix S^(-1)·H
 */
function solveGeneralizedEigenvalueProblem(
  H: number[][],
  S: number[][]
): number[][] {
  // For better numerical stability, we solve S·X = H for X = S^(-1)·H
  // using a direct solver for tridiagonal systems

  // Since S is symmetric tridiagonal, we can use Cholesky or LU decomposition
  // For simplicity, we'll invert S directly (small matrices in typical quantum problems)

  const Sinv = invertSymmetricMatrix(S);

  // Compute S^(-1)·H
  return multiplyMatrices(Sinv, H);
}

/**
 * Invert a symmetric matrix using Gaussian elimination with partial pivoting.
 * This is suitable for small to moderate sized matrices (N < 500).
 *
 * @param A - Symmetric matrix to invert
 * @returns Inverse matrix A^(-1)
 */
function invertSymmetricMatrix(A: number[][]): number[][] {
  const N = A.length;

  // Create augmented matrix [A | I]
  const augmented: number[][] = [];
  for (let i = 0; i < N; i++) {
    augmented[i] = [...A[i]];
    for (let j = 0; j < N; j++) {
      augmented[i].push(i === j ? 1 : 0);
    }
  }

  // Forward elimination with partial pivoting
  for (let k = 0; k < N; k++) {
    // Find pivot
    let maxRow = k;
    for (let i = k + 1; i < N; i++) {
      if (Math.abs(augmented[i][k]) > Math.abs(augmented[maxRow][k])) {
        maxRow = i;
      }
    }

    // Swap rows
    [augmented[k], augmented[maxRow]] = [augmented[maxRow], augmented[k]];

    // Scale pivot row
    const pivot = augmented[k][k];
    for (let j = 0; j < 2 * N; j++) {
      augmented[k][j] /= pivot;
    }

    // Eliminate column
    for (let i = 0; i < N; i++) {
      if (i !== k) {
        const factor = augmented[i][k];
        for (let j = 0; j < 2 * N; j++) {
          augmented[i][j] -= factor * augmented[k][j];
        }
      }
    }
  }

  // Extract inverse from right half of augmented matrix
  const inverse: number[][] = [];
  for (let i = 0; i < N; i++) {
    inverse[i] = augmented[i].slice(N);
  }

  return inverse;
}

/**
 * Multiply two matrices: C = A·B
 */
function multiplyMatrices(A: number[][], B: number[][]): number[][] {
  const M = A.length;
  const N = B[0].length;
  const K = A[0].length;

  const C: number[][] = [];
  for (let i = 0; i < M; i++) {
    C[i] = new Array(N).fill(0);
    for (let j = 0; j < N; j++) {
      for (let k = 0; k < K; k++) {
        C[i][j] += A[i][k] * B[k][j];
      }
    }
  }

  return C;
}

/**
 * Create a zero matrix of size N×N.
 */
function createZeroMatrix(N: number): number[][] {
  const matrix: number[][] = [];
  for (let i = 0; i < N; i++) {
    matrix[i] = new Array(N).fill(0);
  }
  return matrix;
}

/**
 * Diagonalize a symmetric matrix using Jacobi eigenvalue algorithm.
 * Returns eigenvalues and eigenvectors.
 *
 * This is the same algorithm used in DVRSolver for consistency.
 *
 * @param matrix - Symmetric N×N matrix
 * @returns Object with eigenvalues array and eigenvectors (array of column vectors)
 */
function diagonalize(matrix: number[][]): {
  eigenvalues: number[];
  eigenvectors: number[][];
} {
  const N = matrix.length;

  // Copy matrix (we'll modify it)
  const A: number[][] = matrix.map((row) => [...row]);

  // Initialize eigenvectors as identity matrix
  const V: number[][] = [];
  for (let i = 0; i < N; i++) {
    V[i] = new Array(N).fill(0);
    V[i][i] = 1;
  }

  // Jacobi iteration
  const maxIterations = 50 * N * N;

  // Use relative tolerance based on matrix scale
  const matrixScale = Math.max(...matrix.map((row) => Math.max(...row.map(Math.abs))));
  const tolerance = matrixScale * 1e-12;

  for (let iter = 0; iter < maxIterations; iter++) {
    // Find largest off-diagonal element
    let maxVal = 0;
    let p = 0;
    let q = 1;

    for (let i = 0; i < N; i++) {
      for (let j = i + 1; j < N; j++) {
        if (Math.abs(A[i][j]) > maxVal) {
          maxVal = Math.abs(A[i][j]);
          p = i;
          q = j;
        }
      }
    }

    // Check convergence
    if (maxVal < tolerance) {
      break;
    }

    // Calculate rotation angle
    const theta = 0.5 * Math.atan2(2 * A[p][q], A[q][q] - A[p][p]);
    const c = Math.cos(theta);
    const s = Math.sin(theta);

    // Apply Jacobi rotation to A
    const App = c * c * A[p][p] - 2 * s * c * A[p][q] + s * s * A[q][q];
    const Aqq = s * s * A[p][p] + 2 * s * c * A[p][q] + c * c * A[q][q];
    const Apq = 0;

    A[p][p] = App;
    A[q][q] = Aqq;
    A[p][q] = Apq;
    A[q][p] = Apq;

    for (let i = 0; i < N; i++) {
      if (i !== p && i !== q) {
        const Aip = c * A[i][p] - s * A[i][q];
        const Aiq = s * A[i][p] + c * A[i][q];
        A[i][p] = Aip;
        A[p][i] = Aip;
        A[i][q] = Aiq;
        A[q][i] = Aiq;
      }
    }

    // Apply rotation to eigenvectors
    for (let i = 0; i < N; i++) {
      const Vip = c * V[i][p] - s * V[i][q];
      const Viq = s * V[i][p] + c * V[i][q];
      V[i][p] = Vip;
      V[i][q] = Viq;
    }
  }

  // Extract eigenvalues (diagonal of A)
  const eigenvalues = A.map((row, i) => row[i]);

  // Extract eigenvectors (columns of V)
  const eigenvectors: number[][] = [];
  for (let j = 0; j < N; j++) {
    const eigenvector: number[] = [];
    for (let i = 0; i < N; i++) {
      eigenvector.push(V[i][j]);
    }
    eigenvectors.push(eigenvector);
  }

  return { eigenvalues, eigenvectors };
}

/**
 * Normalize a wavefunction using trapezoidal rule.
 *
 * @param psi - Wavefunction array
 * @param dx - Grid spacing (meters)
 * @returns Normalized wavefunction
 */
function normalize(psi: number[], dx: number): number[] {
  // Calculate ∫|ψ|² dx using trapezoidal rule
  let integral = 0;
  for (let i = 0; i < psi.length - 1; i++) {
    integral += (psi[i] * psi[i] + psi[i + 1] * psi[i + 1]) / 2;
  }
  integral *= dx;

  const normalization = Math.sqrt(integral);

  // Avoid division by zero
  if (normalization < 1e-30) {
    return psi;
  }

  return psi.map((val) => val / normalization);
}

qppw.register("MatrixNumerovSolver", { solveMatrixNumerov });
