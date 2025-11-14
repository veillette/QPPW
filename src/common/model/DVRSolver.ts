/**
 * Discrete Variable Representation (DVR) method for solving the 1D
 * time-independent Schrödinger equation.
 *
 * This method constructs the Hamiltonian matrix H = T + V and diagonalizes it:
 * - Potential energy V is diagonal: V_ii = V(x_i)
 * - Kinetic energy T uses Colbert-Miller formula for uniform grids
 *
 * Reference: Colbert & Miller, J. Chem. Phys. 96, 1982 (1992)
 */

import QuantumConstants from "./QuantumConstants.js";
import { BoundStateResult, GridConfig, PotentialFunction } from "./PotentialFunction.js";
import qppw from "../../QPPWNamespace.js";

/**
 * Solve the 1D Schrödinger equation using DVR method with Colbert-Miller kinetic energy.
 *
 * @param potential - Function V(x) that returns potential energy in Joules
 * @param mass - Particle mass in kg
 * @param numStates - Number of lowest bound states to return
 * @param gridConfig - Grid configuration
 * @returns Bound state results
 */
export function solveDVR(
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

  const N = numPoints;

  // Construct potential energy matrix V (diagonal)
  const potentialEnergyMatrix = createDiagonalMatrix(N);
  for (let i = 0; i < N; i++) {
    potentialEnergyMatrix[i][i] = potential(xGrid[i]);
  }

  // Construct kinetic energy matrix T using Colbert-Miller formula
  const kineticEnergyMatrix = createKineticEnergyMatrix(N, dx, mass);

  // Construct Hamiltonian H = T + V
  const hamiltonianMatrix = addMatrices(kineticEnergyMatrix, potentialEnergyMatrix);

  // Diagonalize Hamiltonian to get eigenvalues (energies) and eigenvectors (wavefunctions)
  const eigen = diagonalize(hamiltonianMatrix);

  // Sort eigenvalues and eigenvectors by energy
  const sortedIndices = eigen.eigenvalues
    .map((e, i) => ({ energy: e, index: i }))
    .sort((a, b) => a.energy - b.energy)
    .map((item) => item.index);

  // Extract the lowest numStates bound states
  const energies: number[] = [];
  const wavefunctions: number[][] = [];

  for (let i = 0; i < Math.min(numStates, N); i++) {
    const idx = sortedIndices[i];
    const energy = eigen.eigenvalues[idx];

    // Only include bound states (energy < V at boundaries)
    const V_boundary = Math.max(potential(xMin), potential(xMax));
    if (energy < V_boundary) {
      energies.push(energy);

      // Normalize wavefunction
      const wavefunction = eigen.eigenvectors[idx];
      const normalizedPsi = normalize(wavefunction, dx);
      wavefunctions.push(normalizedPsi);
    }
  }

  return {
    energies,
    wavefunctions,
    xGrid,
    method: "dvr",
  };
}

/**
 * Create kinetic energy matrix using Colbert-Miller formula:
 * T_ij = (ℏ²/2m) * (1/Δx²) * { π²/3 for i=j, 2*(-1)^(i-j)/(i-j)² for i≠j }
 *
 * @param N - Matrix dimension (number of grid points)
 * @param dx - Grid spacing (meters)
 * @param mass - Particle mass (kg)
 * @returns N×N kinetic energy matrix
 */
function createKineticEnergyMatrix(N: number, dx: number, mass: number): number[][] {
  const { HBAR } = QuantumConstants;
  const prefactor = (HBAR * HBAR) / (2 * mass * dx * dx);

  const T: number[][] = [];
  for (let i = 0; i < N; i++) {
    T[i] = [];
    for (let j = 0; j < N; j++) {
      if (i === j) {
        // Diagonal elements
        T[i][j] = prefactor * (Math.PI * Math.PI) / 3;
      } else {
        // Off-diagonal elements
        const diff = i - j;
        const sign = Math.pow(-1, diff);
        T[i][j] = prefactor * (2 * sign) / (diff * diff);
      }
    }
  }

  return T;
}

/**
 * Create a diagonal matrix of size N×N initialized to zero.
 */
function createDiagonalMatrix(N: number): number[][] {
  const matrix: number[][] = [];
  for (let i = 0; i < N; i++) {
    matrix[i] = new Array(N).fill(0);
  }
  return matrix;
}

/**
 * Add two matrices element-wise.
 */
function addMatrices(A: number[][], B: number[][]): number[][] {
  const N = A.length;
  const C: number[][] = [];
  for (let i = 0; i < N; i++) {
    C[i] = [];
    for (let j = 0; j < N; j++) {
      C[i][j] = A[i][j] + B[i][j];
    }
  }
  return C;
}

/**
 * Diagonalize a symmetric matrix using Jacobi eigenvalue algorithm.
 * Returns eigenvalues and eigenvectors.
 *
 * Note: For production code, consider using a more robust library like numeric.js or ml-matrix.
 * This implementation is sufficient for moderate matrix sizes.
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
  const tolerance = 1e-12;

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
    const theta =
      0.5 * Math.atan2(2 * A[p][q], A[q][q] - A[p][p]);
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
  return psi.map((val) => val / normalization);
}

qppw.register("DVRSolver", { solveDVR });
