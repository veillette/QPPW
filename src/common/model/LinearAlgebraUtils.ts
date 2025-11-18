/**
 * Linear algebra utilities compatible with scenerystack/dot API.
 *
 * This module provides shared matrix operations and eigenvalue decomposition
 * for the quantum mechanics solvers. The Matrix class API mirrors scenerystack/dot's
 * Matrix class for compatibility.
 *
 * In browser: Uses scenerystack/dot Matrix class directly
 * In Node: Uses compatible implementation for testing
 */

import qppw from "../../QPPWNamespace.js";

/**
 * Matrix class compatible with scenerystack/dot Matrix API.
 * Provides the same methods used by the quantum mechanics solvers.
 */
export class DotMatrix {
  private m: number; // rows
  private n: number; // columns
  private entries: Float64Array;

  /**
   * Construct a matrix.
   *
   * @param m - Number of rows
   * @param n - Number of columns
   * @param entries - Optional initial entries (row-major order)
   */
  constructor(m: number, n: number, entries?: number[] | Float64Array) {
    this.m = m;
    this.n = n;
    if (entries) {
      this.entries = new Float64Array(entries);
    } else {
      this.entries = new Float64Array(m * n);
    }
  }

  getRowDimension(): number {
    return this.m;
  }

  getColumnDimension(): number {
    return this.n;
  }

  /**
   * Get element at (i, j).
   */
  get(i: number, j: number): number {
    return this.entries[i * this.n + j];
  }

  /**
   * Set element at (i, j).
   */
  set(i: number, j: number, value: number): void {
    this.entries[i * this.n + j] = value;
  }

  /**
   * Get a submatrix.
   *
   * @param i0 - Initial row index
   * @param i1 - Final row index
   * @param j0 - Initial column index
   * @param j1 - Final column index
   * @returns Submatrix
   */
  getMatrix(i0: number, i1: number, j0: number, j1: number): DotMatrix {
    const rows = i1 - i0 + 1;
    const cols = j1 - j0 + 1;
    const result = new DotMatrix(rows, cols);
    for (let i = 0; i < rows; i++) {
      for (let j = 0; j < cols; j++) {
        result.set(i, j, this.get(i0 + i, j0 + j));
      }
    }
    return result;
  }

  /**
   * Create a copy of this matrix.
   */
  copy(): DotMatrix {
    return new DotMatrix(this.m, this.n, this.entries);
  }

  /**
   * Matrix transpose.
   */
  transpose(): DotMatrix {
    const result = new DotMatrix(this.n, this.m);
    for (let i = 0; i < this.m; i++) {
      for (let j = 0; j < this.n; j++) {
        result.set(j, i, this.get(i, j));
      }
    }
    return result;
  }

  /**
   * Matrix addition.
   */
  plus(B: DotMatrix): DotMatrix {
    const result = new DotMatrix(this.m, this.n);
    for (let i = 0; i < this.m; i++) {
      for (let j = 0; j < this.n; j++) {
        result.set(i, j, this.get(i, j) + B.get(i, j));
      }
    }
    return result;
  }

  /**
   * Multiply matrix by scalar in place.
   */
  timesEquals(scalar: number): DotMatrix {
    for (let i = 0; i < this.entries.length; i++) {
      this.entries[i] *= scalar;
    }
    return this;
  }

  /**
   * Matrix multiplication.
   */
  times(B: DotMatrix): DotMatrix {
    const result = new DotMatrix(this.m, B.n);
    for (let i = 0; i < this.m; i++) {
      for (let j = 0; j < B.n; j++) {
        let sum = 0;
        for (let k = 0; k < this.n; k++) {
          sum += this.get(i, k) * B.get(k, j);
        }
        result.set(i, j, sum);
      }
    }
    return result;
  }

  /**
   * Matrix inverse using Gaussian elimination with partial pivoting.
   */
  inverse(): DotMatrix {
    if (this.m !== this.n) {
      throw new Error("Matrix must be square for inverse");
    }
    const N = this.m;

    // Create augmented matrix [A | I]
    const aug = new DotMatrix(N, 2 * N);
    for (let i = 0; i < N; i++) {
      for (let j = 0; j < N; j++) {
        aug.set(i, j, this.get(i, j));
        aug.set(i, j + N, i === j ? 1 : 0);
      }
    }

    // Gaussian elimination with partial pivoting
    for (let col = 0; col < N; col++) {
      // Find pivot
      let maxVal = Math.abs(aug.get(col, col));
      let maxRow = col;
      for (let row = col + 1; row < N; row++) {
        const absVal = Math.abs(aug.get(row, col));
        if (absVal > maxVal) {
          maxVal = absVal;
          maxRow = row;
        }
      }

      // Swap rows
      if (maxRow !== col) {
        for (let j = 0; j < 2 * N; j++) {
          const temp = aug.get(col, j);
          aug.set(col, j, aug.get(maxRow, j));
          aug.set(maxRow, j, temp);
        }
      }

      // Scale pivot row
      const pivot = aug.get(col, col);
      if (Math.abs(pivot) < 1e-14) {
        throw new Error("Matrix is singular");
      }
      for (let j = 0; j < 2 * N; j++) {
        aug.set(col, j, aug.get(col, j) / pivot);
      }

      // Eliminate column
      for (let row = 0; row < N; row++) {
        if (row !== col) {
          const factor = aug.get(row, col);
          for (let j = 0; j < 2 * N; j++) {
            aug.set(row, j, aug.get(row, j) - factor * aug.get(col, j));
          }
        }
      }
    }

    // Extract inverse from augmented matrix
    const inv = new DotMatrix(N, N);
    for (let i = 0; i < N; i++) {
      for (let j = 0; j < N; j++) {
        inv.set(i, j, aug.get(i, j + N));
      }
    }
    return inv;
  }

  /**
   * Create identity matrix.
   */
  static identity(m: number, n: number): DotMatrix {
    const result = new DotMatrix(m, n);
    const min = Math.min(m, n);
    for (let i = 0; i < min; i++) {
      result.set(i, i, 1);
    }
    return result;
  }

  /**
   * Create diagonal matrix from array.
   */
  static diagonalMatrix(values: number[]): DotMatrix {
    const N = values.length;
    const result = new DotMatrix(N, N);
    for (let i = 0; i < N; i++) {
      result.set(i, i, values[i]);
    }
    return result;
  }
}

/**
 * Complex number interface for FFT operations
 */
export interface ComplexNumber {
  real: number;
  imaginary: number;
}

/**
 * Result of eigenvalue decomposition
 */
export interface EigenDecompositionResult {
  eigenvalues: number[];
  eigenvectors: number[][];
}

/**
 * Convert a 2D number array to a dot Matrix.
 *
 * @param array - 2D array in row-major order
 * @returns Matrix object
 */
export function arrayToMatrix(array: number[][]): DotMatrix {
  const m = array.length;
  const n = array[0]?.length ?? 0;
  const flat: number[] = [];

  for (let i = 0; i < m; i++) {
    for (let j = 0; j < n; j++) {
      flat.push(array[i][j]);
    }
  }

  return new DotMatrix(m, n, flat);
}

/**
 * Convert a dot Matrix to a 2D number array.
 *
 * @param matrix - Matrix object
 * @returns 2D array in row-major order
 */
export function matrixToArray(matrix: DotMatrix): number[][] {
  const m = matrix.getRowDimension();
  const n = matrix.getColumnDimension();
  const array: number[][] = [];

  for (let i = 0; i < m; i++) {
    array[i] = [];
    for (let j = 0; j < n; j++) {
      array[i][j] = matrix.get(i, j);
    }
  }

  return array;
}

/**
 * Create a zero matrix of size N×N.
 *
 * @param N - Matrix dimension
 * @returns N×N zero matrix
 */
export function createZeroMatrix(N: number): DotMatrix {
  return new DotMatrix(N, N);
}

/**
 * Create an identity matrix of size N×N.
 *
 * @param N - Matrix dimension
 * @returns N×N identity matrix
 */
export function createIdentityMatrix(N: number): DotMatrix {
  return DotMatrix.identity(N, N);
}

/**
 * Create a diagonal matrix from an array of values.
 *
 * @param values - Array of diagonal values
 * @returns Diagonal matrix
 */
export function createDiagonalMatrix(values: number[]): DotMatrix {
  return DotMatrix.diagonalMatrix(values);
}

/**
 * Add two matrices using dot's Matrix.plus().
 *
 * @param A - First matrix
 * @param B - Second matrix
 * @returns Sum A + B
 */
export function addMatrices(A: DotMatrix, B: DotMatrix): DotMatrix {
  return A.plus(B);
}

/**
 * Multiply two matrices using dot's Matrix.times().
 *
 * @param A - First matrix
 * @param B - Second matrix
 * @returns Product A × B
 */
export function multiplyMatrices(A: DotMatrix, B: DotMatrix): DotMatrix {
  return A.times(B);
}

/**
 * Compute the inverse of a matrix using dot's Matrix.inverse().
 *
 * @param A - Matrix to invert
 * @returns Inverse matrix A^(-1)
 */
export function invertMatrix(A: DotMatrix): DotMatrix {
  return A.inverse();
}

/**
 * Symmetrize a matrix by computing (M + M^T)/2.
 * This ensures the matrix is exactly symmetric, which is required for
 * eigenvalue problems involving self-adjoint operators.
 *
 * @param M - Matrix to symmetrize
 * @returns Symmetric matrix
 */
export function symmetrizeMatrix(M: DotMatrix): DotMatrix {
  const transposed = M.transpose();
  const sum = M.plus(transposed);
  return sum.timesEquals(0.5);
}

/**
 * Extract interior matrix (remove first and last rows/columns)
 * for implementing Dirichlet boundary conditions.
 *
 * @param H - Full matrix
 * @returns Interior matrix with boundary rows/columns removed
 */
export function extractInteriorMatrix(H: DotMatrix): DotMatrix {
  const N = H.getRowDimension();
  return H.getMatrix(1, N - 2, 1, N - 2);
}

/**
 * Diagonalize a symmetric matrix using Jacobi eigenvalue algorithm.
 * Returns eigenvalues and eigenvectors.
 *
 * This implementation uses cyclic Jacobi rotations to find all eigenvalues
 * and eigenvectors of a real symmetric matrix.
 *
 * @param matrix - Symmetric N×N matrix (can be Matrix or number[][])
 * @returns Object with eigenvalues array and eigenvectors (array of column vectors)
 */
export function diagonalize(
  matrix: DotMatrix | number[][]
): EigenDecompositionResult {
  // Convert to 2D array for Jacobi algorithm
  const A = Array.isArray(matrix) ? matrix.map(row => [...row]) : matrixToArray(matrix);
  const N = A.length;

  // Initialize eigenvector matrix as identity
  const V: number[][] = [];
  for (let i = 0; i < N; i++) {
    V[i] = [];
    for (let j = 0; j < N; j++) {
      V[i][j] = i === j ? 1 : 0;
    }
  }

  // Jacobi rotation parameters
  const maxIterations = 50;
  const tolerance = 1e-12;

  for (let iter = 0; iter < maxIterations; iter++) {
    // Find largest off-diagonal element
    let maxOffDiag = 0;
    let p = 0;
    let q = 1;

    for (let i = 0; i < N; i++) {
      for (let j = i + 1; j < N; j++) {
        const absVal = Math.abs(A[i][j]);
        if (absVal > maxOffDiag) {
          maxOffDiag = absVal;
          p = i;
          q = j;
        }
      }
    }

    // Check convergence
    if (maxOffDiag < tolerance) {
      break;
    }

    // Compute rotation angle
    const App = A[p][p];
    const Aqq = A[q][q];
    const Apq = A[p][q];

    let theta: number;
    if (Math.abs(App - Aqq) < 1e-30) {
      theta = Math.PI / 4;
    } else {
      theta = 0.5 * Math.atan2(2 * Apq, Aqq - App);
    }

    const c = Math.cos(theta);
    const s = Math.sin(theta);

    // Apply rotation to A
    const newApp = c * c * App - 2 * s * c * Apq + s * s * Aqq;
    const newAqq = s * s * App + 2 * s * c * Apq + c * c * Aqq;

    A[p][p] = newApp;
    A[q][q] = newAqq;
    A[p][q] = 0;
    A[q][p] = 0;

    for (let i = 0; i < N; i++) {
      if (i !== p && i !== q) {
        const Aip = A[i][p];
        const Aiq = A[i][q];
        A[i][p] = c * Aip - s * Aiq;
        A[p][i] = A[i][p];
        A[i][q] = s * Aip + c * Aiq;
        A[q][i] = A[i][q];
      }
    }

    // Update eigenvector matrix
    for (let i = 0; i < N; i++) {
      const Vip = V[i][p];
      const Viq = V[i][q];
      V[i][p] = c * Vip - s * Viq;
      V[i][q] = s * Vip + c * Viq;
    }
  }

  // Extract eigenvalues and eigenvectors
  const eigenvalues: number[] = [];
  const eigenvectors: number[][] = [];

  for (let i = 0; i < N; i++) {
    eigenvalues.push(A[i][i]);

    // Extract column i of V as eigenvector i
    const eigenvector: number[] = [];
    for (let j = 0; j < N; j++) {
      eigenvector.push(V[j][i]);
    }
    eigenvectors.push(eigenvector);
  }

  return { eigenvalues, eigenvectors };
}

/**
 * Normalize a wavefunction using trapezoidal rule.
 * Ensures ∫|ψ|² dx = 1.
 *
 * @param psi - Wavefunction array
 * @param dx - Grid spacing (meters)
 * @returns Normalized wavefunction
 */
export function normalizeWavefunction(psi: number[], dx: number): number[] {
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

/**
 * Normalize a wavefunction using Clenshaw-Curtis quadrature.
 * Used for Chebyshev spectral methods where the grid is non-uniform.
 *
 * @param psi - Wavefunction at Chebyshev points
 * @param xMin - Minimum value of physical domain
 * @param xMax - Maximum value of physical domain
 * @returns Normalized wavefunction
 */
export function normalizeWavefunctionChebyshev(
  psi: number[],
  xMin: number,
  xMax: number
): number[] {
  const N = psi.length;

  // Clenshaw-Curtis quadrature weights for ξ ∈ [-1,1]
  const weights: number[] = [];
  for (let j = 0; j < N; j++) {
    if (j === 0 || j === N - 1) {
      weights.push(Math.PI / (2 * (N - 1)));
    } else {
      weights.push(Math.PI / (N - 1));
    }
  }

  // Jacobian for coordinate transformation: dx = (xMax - xMin)/2 * dξ
  const jacobian = (xMax - xMin) / 2;

  // Compute ∫|ψ|² dx
  let integral = 0;
  for (let i = 0; i < N; i++) {
    integral += weights[i] * psi[i] * psi[i] * jacobian;
  }

  const normalization = Math.sqrt(integral);
  return psi.map((val) => val / normalization);
}

/**
 * Complex number addition.
 *
 * @param a - First complex number
 * @param b - Second complex number
 * @returns Sum a + b
 */
export function complexAdd(a: ComplexNumber, b: ComplexNumber): ComplexNumber {
  return {
    real: a.real + b.real,
    imaginary: a.imaginary + b.imaginary,
  };
}

/**
 * Complex number subtraction.
 *
 * @param a - First complex number
 * @param b - Second complex number
 * @returns Difference a - b
 */
export function complexSubtract(a: ComplexNumber, b: ComplexNumber): ComplexNumber {
  return {
    real: a.real - b.real,
    imaginary: a.imaginary - b.imaginary,
  };
}

/**
 * Complex number multiplication.
 *
 * @param a - First complex number
 * @param b - Second complex number
 * @returns Product a × b
 */
export function complexMultiply(a: ComplexNumber, b: ComplexNumber): ComplexNumber {
  return {
    real: a.real * b.real - a.imaginary * b.imaginary,
    imaginary: a.real * b.imaginary + a.imaginary * b.real,
  };
}

/**
 * Fast Fourier Transform (FFT)
 * Cooley-Tukey radix-2 decimation-in-time algorithm
 *
 * @param x - Input array of complex numbers
 * @returns FFT of input
 */
export function fft(x: ComplexNumber[]): ComplexNumber[] {
  const N = x.length;

  // Base case
  if (N <= 1) {
    return x;
  }

  // Divide
  const even: ComplexNumber[] = [];
  const odd: ComplexNumber[] = [];
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
  const result: ComplexNumber[] = new Array(N);
  for (let frequencyIndex = 0; frequencyIndex < N / 2; frequencyIndex++) {
    const angle = (-2 * Math.PI * frequencyIndex) / N;
    const twiddle: ComplexNumber = {
      real: Math.cos(angle),
      imaginary: Math.sin(angle),
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
export function ifft(X: ComplexNumber[]): ComplexNumber[] {
  const N = X.length;

  // Conjugate input
  const XConj = X.map((x) => ({ real: x.real, imaginary: -x.imaginary }));

  // Apply FFT
  const result = fft(XConj);

  // Conjugate output and scale
  return result.map((x) => ({
    real: x.real / N,
    imaginary: -x.imaginary / N,
  }));
}

/**
 * Generate FFT frequency array similar to numpy.fft.fftfreq
 * Returns k = 2π * freq / dx
 *
 * @param N - Number of points
 * @param dx - Grid spacing
 * @returns Array of wave numbers k
 */
export function fftFreq(N: number, dx: number): number[] {
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

qppw.register("LinearAlgebraUtils", {
  arrayToMatrix,
  matrixToArray,
  createZeroMatrix,
  createIdentityMatrix,
  createDiagonalMatrix,
  addMatrices,
  multiplyMatrices,
  invertMatrix,
  symmetrizeMatrix,
  extractInteriorMatrix,
  diagonalize,
  normalizeWavefunction,
  normalizeWavefunctionChebyshev,
  complexAdd,
  complexSubtract,
  complexMultiply,
  fft,
  ifft,
  fftFreq,
});
