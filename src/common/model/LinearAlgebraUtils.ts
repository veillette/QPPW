/**
 * Linear algebra utilities using the dot library from scenerystack.
 *
 * This module provides shared matrix operations and eigenvalue decomposition
 * for the quantum mechanics solvers, leveraging the dot library's Matrix class
 * for optimized operations.
 *
 * The eigenvalue decomposition uses the Jacobi eigenvalue algorithm, which is
 * well-suited for the symmetric matrices that arise in quantum mechanics.
 */

import qppw from "../../QPPWNamespace.js";

/**
 * Complex number interface for FFT operations
 */
export interface ComplexNumber {
  real: number;
  imaginary: number;
}

/**
 * Matrix class for linear algebra operations.
 * This implementation mirrors the API of scenerystack/dot Matrix class
 * but works in Node.js environments for testing.
 */
export class Matrix {
  private entries: number[];
  private m: number;
  private n: number;

  /**
   * Creates an m×n matrix.
   * @param m - Number of rows
   * @param n - Number of columns
   */
  constructor(m: number, n: number) {
    this.m = m;
    this.n = n;
    this.entries = new Array(m * n).fill(0);
  }

  /**
   * Get number of rows.
   */
  getRowDimension(): number {
    return this.m;
  }

  /**
   * Get number of columns.
   */
  getColumnDimension(): number {
    return this.n;
  }

  /**
   * Get element at row i, column j.
   */
  get(i: number, j: number): number {
    return this.entries[i * this.n + j];
  }

  /**
   * Set element at row i, column j.
   */
  set(i: number, j: number, value: number): void {
    this.entries[i * this.n + j] = value;
  }

  /**
   * Matrix addition.
   */
  plus(matrix: Matrix): Matrix {
    const result = new Matrix(this.m, this.n);
    for (let i = 0; i < this.m; i++) {
      for (let j = 0; j < this.n; j++) {
        result.set(i, j, this.get(i, j) + matrix.get(i, j));
      }
    }
    return result;
  }

  /**
   * Matrix or scalar multiplication.
   */
  times(matrixOrScalar: Matrix | number): Matrix {
    if (typeof matrixOrScalar === 'number') {
      const result = new Matrix(this.m, this.n);
      for (let i = 0; i < this.m; i++) {
        for (let j = 0; j < this.n; j++) {
          result.set(i, j, this.get(i, j) * matrixOrScalar);
        }
      }
      return result;
    } else {
      const B = matrixOrScalar;
      const result = new Matrix(this.m, B.getColumnDimension());
      for (let i = 0; i < this.m; i++) {
        for (let j = 0; j < B.getColumnDimension(); j++) {
          let sum = 0;
          for (let k = 0; k < this.n; k++) {
            sum += this.get(i, k) * B.get(k, j);
          }
          result.set(i, j, sum);
        }
      }
      return result;
    }
  }

  /**
   * In-place scalar multiplication.
   */
  timesEquals(scalar: number): Matrix {
    for (let i = 0; i < this.entries.length; i++) {
      this.entries[i] *= scalar;
    }
    return this;
  }

  /**
   * Matrix transpose.
   */
  transpose(): Matrix {
    const result = new Matrix(this.n, this.m);
    for (let i = 0; i < this.m; i++) {
      for (let j = 0; j < this.n; j++) {
        result.set(j, i, this.get(i, j));
      }
    }
    return result;
  }

  /**
   * Matrix inversion using Gaussian elimination with partial pivoting.
   */
  inverse(): Matrix {
    const N = this.m;

    // Create augmented matrix [A | I]
    const augmented: number[][] = [];
    for (let i = 0; i < N; i++) {
      augmented[i] = [];
      for (let j = 0; j < N; j++) {
        augmented[i][j] = this.get(i, j);
      }
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
    const result = new Matrix(N, N);
    for (let i = 0; i < N; i++) {
      for (let j = 0; j < N; j++) {
        result.set(i, j, augmented[i][N + j]);
      }
    }

    return result;
  }

  /**
   * Creates a diagonal matrix from an array of values.
   */
  static diagonalMatrix(values: number[]): Matrix {
    const N = values.length;
    const result = new Matrix(N, N);
    for (let i = 0; i < N; i++) {
      result.set(i, i, values[i]);
    }
    return result;
  }

  /**
   * Creates an identity matrix.
   */
  static identity(m: number, n: number): Matrix {
    const result = new Matrix(m, n);
    const min = Math.min(m, n);
    for (let i = 0; i < min; i++) {
      result.set(i, i, 1);
    }
    return result;
  }
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
export function arrayToMatrix(array: number[][]): Matrix {
  const m = array.length;
  const n = array[0]?.length ?? 0;
  const matrix = new Matrix(m, n);

  for (let i = 0; i < m; i++) {
    for (let j = 0; j < n; j++) {
      matrix.set(i, j, array[i][j]);
    }
  }

  return matrix;
}

/**
 * Convert a dot Matrix to a 2D number array.
 *
 * @param matrix - Matrix object
 * @returns 2D array in row-major order
 */
export function matrixToArray(matrix: Matrix): number[][] {
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
export function createZeroMatrix(N: number): Matrix {
  return new Matrix(N, N);
}

/**
 * Create an identity matrix of size N×N.
 *
 * @param N - Matrix dimension
 * @returns N×N identity matrix
 */
export function createIdentityMatrix(N: number): Matrix {
  return Matrix.identity(N, N);
}

/**
 * Create a diagonal matrix from an array of values.
 *
 * @param values - Array of diagonal values
 * @returns Diagonal matrix
 */
export function createDiagonalMatrix(values: number[]): Matrix {
  return Matrix.diagonalMatrix(values);
}

/**
 * Add two matrices using dot's Matrix.plus().
 *
 * @param A - First matrix
 * @param B - Second matrix
 * @returns Sum A + B
 */
export function addMatrices(A: Matrix, B: Matrix): Matrix {
  return A.plus(B);
}

/**
 * Multiply two matrices using dot's Matrix.times().
 *
 * @param A - First matrix
 * @param B - Second matrix
 * @returns Product A × B
 */
export function multiplyMatrices(A: Matrix, B: Matrix): Matrix {
  return A.times(B) as Matrix;
}

/**
 * Compute the inverse of a matrix using dot's Matrix.inverse().
 *
 * @param A - Matrix to invert
 * @returns Inverse matrix A^(-1)
 */
export function invertMatrix(A: Matrix): Matrix {
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
export function symmetrizeMatrix(M: Matrix): Matrix {
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
export function extractInteriorMatrix(H: Matrix): Matrix {
  const N = H.getRowDimension();
  const interior = new Matrix(N - 2, N - 2);

  for (let i = 1; i < N - 1; i++) {
    for (let j = 1; j < N - 1; j++) {
      interior.set(i - 1, j - 1, H.get(i, j));
    }
  }

  return interior;
}

/**
 * Diagonalize a symmetric matrix using Jacobi eigenvalue algorithm.
 * Returns eigenvalues and eigenvectors.
 *
 * The Jacobi method is particularly well-suited for:
 * - Symmetric matrices (guaranteed by Hermitian quantum operators)
 * - Dense matrices where all eigenvalues are needed
 * - Situations requiring high accuracy
 *
 * @param matrix - Symmetric N×N matrix (can be Matrix or number[][])
 * @returns Object with eigenvalues array and eigenvectors (array of column vectors)
 */
export function diagonalize(
  matrix: Matrix | number[][]
): EigenDecompositionResult {
  // Convert to array if needed for internal computation
  const inputArray = matrix instanceof Matrix ? matrixToArray(matrix) : matrix;
  const N = inputArray.length;

  // Copy matrix (we'll modify it)
  const A: number[][] = inputArray.map((row: number[]) => [...row]);

  // Initialize eigenvectors as identity matrix
  const V: number[][] = [];
  for (let i = 0; i < N; i++) {
    V[i] = new Array(N).fill(0);
    V[i][i] = 1;
  }

  // Jacobi iteration
  const maxIterations = 50 * N * N;

  // Use relative tolerance based on matrix scale
  // For quantum mechanical energies (~1e-20 J), absolute tolerance
  // causes premature convergence. Use relative tolerance instead.
  const matrixScale = Math.max(
    ...inputArray.map((row: number[]) => Math.max(...row.map(Math.abs)))
  );
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
    const App =
      c * c * A[p][p] - 2 * s * c * A[p][q] + s * s * A[q][q];
    const Aqq =
      s * s * A[p][p] + 2 * s * c * A[p][q] + c * c * A[q][q];
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
