/**
 * Fourier Grid Hamiltonian (FGH) method for solving the 1D
 * time-independent Schrödinger equation.
 *
 * This method uses plane wave basis (Fourier basis):
 * - Kinetic energy is diagonal in momentum space: T_k = ℏ²k²/(2m)
 * - Potential energy is diagonal in position space: V(x)
 * - Hamiltonian is constructed by transforming between spaces using FFT
 *
 * Natural for periodic systems with spectral accuracy.
 *
 * Reference: Marston & Balint-Kurti, J. Chem. Phys. 91, 3571 (1989)
 */

import QuantumConstants from "./QuantumConstants.js";
import { BoundStateResult, GridConfig, PotentialFunction } from "./PotentialFunction.js";
import qppw from "../../QPPWNamespace.js";

/**
 * Complex number representation
 */
interface Complex {
  real: number;
  imag: number;
}

/**
 * Solve the 1D Schrödinger equation using FGH method.
 *
 * @param potential - Function V(x) that returns potential energy in Joules
 * @param mass - Particle mass in kg
 * @param numStates - Number of lowest bound states to return
 * @param gridConfig - Grid configuration
 * @returns Bound state results
 */
export function solveFGH(
  potential: PotentialFunction,
  mass: number,
  numStates: number,
  gridConfig: GridConfig,
): BoundStateResult {
  const { xMin, xMax, numPoints } = gridConfig;
  // L = (xMax - xMin) / 2 for periodic domain [-L, L]
  const dx = (xMax - xMin) / numPoints; // Note: periodic grid, no endpoint
  const N = numPoints;

  // Generate grid (periodic, exclude endpoint)
  const xGrid: number[] = [];
  for (let i = 0; i < N; i++) {
    xGrid.push(xMin + i * dx);
  }

  // Generate momentum space grid
  const waveNumberGrid = fftFreq(N, dx);

  // Kinetic energy in momentum space: T_k = ℏ²k²/(2m)
  const { HBAR } = QuantumConstants;
  const T_k = waveNumberGrid.map((waveNumber) => (HBAR * HBAR * waveNumber * waveNumber) / (2 * mass));

  // Potential energy in position space
  const V_x = xGrid.map(potential);

  // Build Hamiltonian matrix using FGH method
  // H_ij = ∫ φ_i*(x) H φ_j(x) dx
  // where φ_j(x) are plane waves
  const H = buildFGHHamiltonian(N, T_k, V_x);

  // Diagonalize Hamiltonian
  const eigen = diagonalize(H);

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
    method: "dvr", // Using "dvr" for compatibility with existing code
  };
}

/**
 * Build the FGH Hamiltonian matrix.
 * H = T + V, where T is kinetic energy (non-diagonal due to FFT transform)
 * and V is potential energy (diagonal).
 *
 * @param N - Matrix dimension
 * @param T_k - Kinetic energy in momentum space
 * @param V_x - Potential energy in position space
 * @returns N×N Hamiltonian matrix
 */
function buildFGHHamiltonian(N: number, T_k: number[], V_x: number[]): number[][] {
  const H: number[][] = [];

  // Initialize matrix
  for (let i = 0; i < N; i++) {
    H[i] = new Array(N).fill(0);
  }

  // Build kinetic energy matrix by applying T_k in momentum space
  // For each basis function (column j), apply kinetic energy operator
  for (let j = 0; j < N; j++) {
    // Create basis vector in position space (delta function at grid point j)
    const psi_x: Complex[] = [];
    for (let i = 0; i < N; i++) {
      psi_x.push({ real: i === j ? 1.0 : 0.0, imag: 0.0 });
    }

    // Transform to momentum space
    const psi_k = fft(psi_x);

    // Apply kinetic energy in momentum space
    for (let i = 0; i < N; i++) {
      psi_k[i].real *= T_k[i];
      psi_k[i].imag *= T_k[i];
    }

    // Transform back to position space
    const T_psi_x = ifft(psi_k);

    // Store in kinetic energy contribution to H
    for (let i = 0; i < N; i++) {
      H[i][j] = T_psi_x[i].real;
    }
  }

  // Add potential energy (diagonal)
  for (let i = 0; i < N; i++) {
    H[i][i] += V_x[i];
  }

  return H;
}

/**
 * Generate FFT frequency array similar to numpy.fft.fftfreq
 * Returns k = 2π * freq / dx
 *
 * @param N - Number of points
 * @param dx - Grid spacing
 * @returns Array of wave numbers k
 */
function fftFreq(N: number, dx: number): number[] {
  const waveNumberArray: number[] = [];
  const dk = (2 * Math.PI) / (N * dx);

  for (let i = 0; i < N; i++) {
    if (i < N / 2) {
      waveNumberArray.push(i * dk);
    } else {
      waveNumberArray.push((i - N) * dk);
    }
  }

  return waveNumberArray;
}

/**
 * Fast Fourier Transform (FFT)
 * Cooley-Tukey radix-2 decimation-in-time algorithm
 *
 * @param x - Input array of complex numbers
 * @returns FFT of input
 */
function fft(x: Complex[]): Complex[] {
  const N = x.length;

  // Base case
  if (N <= 1) {
    return x;
  }

  // Divide
  const even: Complex[] = [];
  const odd: Complex[] = [];
  for (let i = 0; i < N; i++) {
    if (i % 2 === 0) {
      even.push(x[i]);
    } else {
      odd.push(x[i]);
    }
  }

  // Conquer
  const fftEven = fft(even);
  const fftOdd = fft(odd);

  // Combine
  const result: Complex[] = new Array(N);
  for (let frequencyIndex = 0; frequencyIndex < N / 2; frequencyIndex++) {
    const angle = (-2 * Math.PI * frequencyIndex) / N;
    const twiddle: Complex = {
      real: Math.cos(angle),
      imag: Math.sin(angle),
    };

    const temp = complexMultiply(twiddle, fftOdd[frequencyIndex]);
    result[frequencyIndex] = complexAdd(fftEven[frequencyIndex], temp);
    result[frequencyIndex + N / 2] = complexSubtract(fftEven[frequencyIndex], temp);
  }

  return result;
}

/**
 * Inverse Fast Fourier Transform (IFFT)
 *
 * @param X - Input array of complex numbers
 * @returns IFFT of input
 */
function ifft(X: Complex[]): Complex[] {
  const N = X.length;

  // Conjugate input
  const XConj = X.map((x) => ({ real: x.real, imag: -x.imag }));

  // Apply FFT
  const result = fft(XConj);

  // Conjugate output and scale
  return result.map((x) => ({
    real: x.real / N,
    imag: -x.imag / N,
  }));
}

/**
 * Complex number addition
 */
function complexAdd(a: Complex, b: Complex): Complex {
  return {
    real: a.real + b.real,
    imag: a.imag + b.imag,
  };
}

/**
 * Complex number subtraction
 */
function complexSubtract(a: Complex, b: Complex): Complex {
  return {
    real: a.real - b.real,
    imag: a.imag - b.imag,
  };
}

/**
 * Complex number multiplication
 */
function complexMultiply(a: Complex, b: Complex): Complex {
  return {
    real: a.real * b.real - a.imag * b.imag,
    imag: a.real * b.imag + a.imag * b.real,
  };
}

/**
 * Diagonalize a symmetric matrix using Jacobi eigenvalue algorithm.
 * (Reused from DVRSolver - could be extracted to shared utility)
 *
 * @param matrix - Symmetric N×N matrix
 * @returns Object with eigenvalues array and eigenvectors
 */
function diagonalize(matrix: number[][]): {
  eigenvalues: number[];
  eigenvectors: number[][];
} {
  const N = matrix.length;

  // Copy matrix
  const A: number[][] = matrix.map((row) => [...row]);

  // Initialize eigenvectors as identity matrix
  const V: number[][] = [];
  for (let i = 0; i < N; i++) {
    V[i] = new Array(N).fill(0);
    V[i][i] = 1;
  }

  // Jacobi iteration
  const maxIterations = 50 * N * N;

  // CRITICAL FIX: Use relative tolerance based on matrix scale
  // For quantum mechanical energies (~1e-20 J), absolute tolerance of 1e-12
  // causes premature convergence. Use relative tolerance instead.
  const matrixScale = Math.max(...matrix.map((row) => Math.max(...row.map(Math.abs))));
  const tolerance = matrixScale * 1e-12;

  for (let iter = 0; iter < maxIterations; iter++) {
    // Find largest off-diagonal element
    let maxVal = 0;
    let rowIndex = 0;
    let colIndex = 1;

    for (let i = 0; i < N; i++) {
      for (let j = i + 1; j < N; j++) {
        if (Math.abs(A[i][j]) > maxVal) {
          maxVal = Math.abs(A[i][j]);
          rowIndex = i;
          colIndex = j;
        }
      }
    }

    // Check convergence
    if (maxVal < tolerance) {
      break;
    }

    // Calculate rotation angle
    const theta = 0.5 * Math.atan2(2 * A[rowIndex][colIndex], A[colIndex][colIndex] - A[rowIndex][rowIndex]);
    const cosineTheta = Math.cos(theta);
    const sineTheta = Math.sin(theta);

    // Apply Jacobi rotation
    const App = cosineTheta * cosineTheta * A[rowIndex][rowIndex] - 2 * sineTheta * cosineTheta * A[rowIndex][colIndex] + sineTheta * sineTheta * A[colIndex][colIndex];
    const Aqq = sineTheta * sineTheta * A[rowIndex][rowIndex] + 2 * sineTheta * cosineTheta * A[rowIndex][colIndex] + cosineTheta * cosineTheta * A[colIndex][colIndex];
    const Apq = 0;

    A[rowIndex][rowIndex] = App;
    A[colIndex][colIndex] = Aqq;
    A[rowIndex][colIndex] = Apq;
    A[colIndex][rowIndex] = Apq;

    for (let i = 0; i < N; i++) {
      if (i !== rowIndex && i !== colIndex) {
        const Aip = cosineTheta * A[i][rowIndex] - sineTheta * A[i][colIndex];
        const Aiq = sineTheta * A[i][rowIndex] + cosineTheta * A[i][colIndex];
        A[i][rowIndex] = Aip;
        A[rowIndex][i] = Aip;
        A[i][colIndex] = Aiq;
        A[colIndex][i] = Aiq;
      }
    }

    // Apply rotation to eigenvectors
    for (let i = 0; i < N; i++) {
      const Vip = cosineTheta * V[i][rowIndex] - sineTheta * V[i][colIndex];
      const Viq = sineTheta * V[i][rowIndex] + cosineTheta * V[i][colIndex];
      V[i][rowIndex] = Vip;
      V[i][colIndex] = Viq;
    }
  }

  // Extract eigenvalues
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
 * @param dx - Grid spacing
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

qppw.register("FGHSolver", { solveFGH });
