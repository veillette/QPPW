/**
 * Standalone accuracy test runner for spectral method
 * This file can run directly without building the full project
 */

// Physical constants
const HBAR = 1.054571817e-34; // Reduced Planck constant (J·s)
const ELECTRON_MASS = 9.1093837015e-31; // Electron mass (kg)

/**
 * Calculate percentage error
 */
function percentageError(numerical, analytical) {
  return Math.abs((numerical - analytical) / analytical) * 100;
}

/**
 * Convert Joules to eV
 */
function joulesToEV(joules) {
  return joules / 1.602176634e-19;
}

/**
 * Chebyshev-Gauss-Lobatto points
 */
function chebyshevPoints(N) {
  const points = [];
  for (let j = 0; j < N; j++) {
    points.push(-Math.cos((Math.PI * j) / (N - 1)));
  }
  return points;
}

/**
 * Chebyshev differentiation matrix
 */
function chebyshevDifferentiationMatrix(N) {
  const x = chebyshevPoints(N);
  const c = new Array(N).fill(1);
  c[0] = 2;
  c[N - 1] = 2;

  const D = [];
  for (let i = 0; i < N; i++) {
    D[i] = [];
    for (let j = 0; j < N; j++) {
      if (i === j) {
        if (i === 0) {
          D[i][j] = (2 * (N - 1) * (N - 1) + 1) / 6;
        } else if (i === N - 1) {
          D[i][j] = -(2 * (N - 1) * (N - 1) + 1) / 6;
        } else {
          D[i][j] = -x[j] / (2 * (1 - x[j] * x[j]));
        }
      } else {
        const sign = Math.pow(-1, i + j);
        D[i][j] = (c[i] / c[j]) * sign / (x[i] - x[j]);
      }
    }
  }
  return D;
}

/**
 * Matrix multiplication
 */
function matrixMultiply(A, B) {
  const N = A.length;
  const C = [];
  for (let i = 0; i < N; i++) {
    C[i] = [];
    for (let j = 0; j < N; j++) {
      let sum = 0;
      for (let k = 0; k < N; k++) {
        sum += A[i][k] * B[k][j];
      }
      C[i][j] = sum;
    }
  }
  return C;
}

/**
 * Extract interior matrix (remove boundaries)
 */
function extractInterior(H) {
  const N = H.length;
  const interior = [];
  for (let i = 1; i < N - 1; i++) {
    interior[i - 1] = [];
    for (let j = 1; j < N - 1; j++) {
      interior[i - 1].push(H[i][j]);
    }
  }
  return interior;
}

/**
 * Symmetrize matrix: (M + M^T)/2
 */
function symmetrize(M) {
  const N = M.length;
  const S = [];
  for (let i = 0; i < N; i++) {
    S[i] = [];
    for (let j = 0; j < N; j++) {
      S[i][j] = (M[i][j] + M[j][i]) / 2;
    }
  }
  return S;
}

/**
 * Jacobi diagonalization
 */
function diagonalize(matrix) {
  const N = matrix.length;
  const A = matrix.map(row => [...row]);
  const V = [];
  for (let i = 0; i < N; i++) {
    V[i] = new Array(N).fill(0);
    V[i][i] = 1;
  }

  const matrixScale = Math.max(...matrix.map(row => Math.max(...row.map(Math.abs))));
  const tolerance = matrixScale * 1e-12;
  const maxIterations = 50 * N * N;

  for (let iter = 0; iter < maxIterations; iter++) {
    let maxVal = 0;
    let p = 0, q = 1;

    for (let i = 0; i < N; i++) {
      for (let j = i + 1; j < N; j++) {
        if (Math.abs(A[i][j]) > maxVal) {
          maxVal = Math.abs(A[i][j]);
          p = i;
          q = j;
        }
      }
    }

    if (maxVal < tolerance) break;

    const theta = 0.5 * Math.atan2(2 * A[p][q], A[q][q] - A[p][p]);
    const c = Math.cos(theta);
    const s = Math.sin(theta);

    const App = c * c * A[p][p] - 2 * s * c * A[p][q] + s * s * A[q][q];
    const Aqq = s * s * A[p][p] + 2 * s * c * A[p][q] + c * c * A[q][q];

    A[p][p] = App;
    A[q][q] = Aqq;
    A[p][q] = 0;
    A[q][p] = 0;

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

    for (let i = 0; i < N; i++) {
      const Vip = c * V[i][p] - s * V[i][q];
      const Viq = s * V[i][p] + c * V[i][q];
      V[i][p] = Vip;
      V[i][q] = Viq;
    }
  }

  const eigenvalues = A.map((row, i) => row[i]);
  const eigenvectors = [];
  for (let j = 0; j < N; j++) {
    const eigenvector = [];
    for (let i = 0; i < N; i++) {
      eigenvector.push(V[i][j]);
    }
    eigenvectors.push(eigenvector);
  }

  return { eigenvalues, eigenvectors };
}

/**
 * Solve using spectral method
 */
function solveSpectral(potential, mass, numStates, xMin, xMax, N) {
  // Generate Chebyshev points
  const xiGrid = chebyshevPoints(N);
  const xGrid = xiGrid.map(xi => ((xMax - xMin) * xi + (xMax + xMin)) / 2);

  // Get differentiation matrix and compute D²
  const D = chebyshevDifferentiationMatrix(N);
  const D2 = matrixMultiply(D, D);

  // Apply domain scaling (this is the fix - no empirical correction!)
  const scale = 2 / (xMax - xMin);
  for (let i = 0; i < N; i++) {
    for (let j = 0; j < N; j++) {
      D2[i][j] *= scale * scale;
    }
  }

  // Build Hamiltonian
  const kineticPrefactor = -(HBAR * HBAR) / (2 * mass);
  const H = [];
  for (let i = 0; i < N; i++) {
    H[i] = [];
    for (let j = 0; j < N; j++) {
      H[i][j] = kineticPrefactor * D2[i][j];
    }
    H[i][i] += potential(xGrid[i]);
  }

  // Extract interior, symmetrize, and diagonalize
  const H_interior = extractInterior(H);
  const H_symmetric = symmetrize(H_interior);
  const { eigenvalues } = diagonalize(H_symmetric);

  // Sort by energy
  const sorted = eigenvalues
    .map((e, i) => ({ energy: e, index: i }))
    .sort((a, b) => a.energy - b.energy);

  const energies = sorted.slice(0, numStates).map(s => s.energy);
  return { energies, xGrid };
}

/**
 * Analytical solution: Harmonic oscillator
 */
function harmonicOscillatorEnergies(springConstant, mass, numStates) {
  const omega = Math.sqrt(springConstant / mass);
  const energies = [];
  for (let n = 0; n < numStates; n++) {
    energies.push(HBAR * omega * (n + 0.5));
  }
  return energies;
}

/**
 * Analytical solution: Infinite square well
 */
function infiniteWellEnergies(wellWidth, mass, numStates) {
  const energies = [];
  for (let n = 1; n <= numStates; n++) {
    energies.push((n * n * Math.PI * Math.PI * HBAR * HBAR) / (2 * mass * wellWidth * wellWidth));
  }
  return energies;
}

// ============ RUN TESTS ============

console.log('========================================');
console.log('Spectral Method Accuracy Tests');
console.log('Testing the FIX (no empirical correction)');
console.log('========================================\n');

// Test 1: Harmonic Oscillator
console.log('=== Test 1: Harmonic Oscillator ===');
const omega = 1e15;
const springConstant = ELECTRON_MASS * omega * omega;
const numStates = 5;

const hoNumerical = solveSpectral(
  x => 0.5 * springConstant * x * x,
  ELECTRON_MASS,
  numStates,
  -5e-9,
  5e-9,
  200
);

const hoAnalytical = harmonicOscillatorEnergies(springConstant, ELECTRON_MASS, numStates);

let hoAllPassed = true;
let hoMaxError = 0;

console.log('Energy levels:');
for (let n = 0; n < numStates; n++) {
  const error = percentageError(hoNumerical.energies[n], hoAnalytical[n]);
  const passed = error < 1.0;
  hoAllPassed = hoAllPassed && passed;
  hoMaxError = Math.max(hoMaxError, error);

  const status = passed ? '✓ PASS' : '✗ FAIL';
  console.log(`  E_${n}: ${status} - Numerical: ${joulesToEV(hoNumerical.energies[n]).toFixed(6)} eV, ` +
              `Analytical: ${joulesToEV(hoAnalytical[n]).toFixed(6)} eV, Error: ${error.toFixed(4)}%`);
}

console.log(`Status: ${hoAllPassed ? '✓ PASSED' : '✗ FAILED'} (Max error: ${hoMaxError.toFixed(4)}%)\n`);

// Test 2: Infinite Square Well
console.log('=== Test 2: Infinite Square Well ===');
const wellWidth = 1e-9;

const iswNumerical = solveSpectral(
  () => 0,
  ELECTRON_MASS,
  numStates,
  -wellWidth / 2,
  wellWidth / 2,
  200  // Increased from 150 to 200 for better accuracy
);

const iswAnalytical = infiniteWellEnergies(wellWidth, ELECTRON_MASS, numStates);

let iswAllPassed = true;
let iswMaxError = 0;

console.log('Energy levels:');
for (let n = 0; n < numStates; n++) {
  const error = percentageError(iswNumerical.energies[n], iswAnalytical[n]);
  const passed = error < 1.0;
  iswAllPassed = iswAllPassed && passed;
  iswMaxError = Math.max(iswMaxError, error);

  const status = passed ? '✓ PASS' : '✗ FAIL';
  console.log(`  E_${n + 1}: ${status} - Numerical: ${joulesToEV(iswNumerical.energies[n]).toFixed(6)} eV, ` +
              `Analytical: ${joulesToEV(iswAnalytical[n]).toFixed(6)} eV, Error: ${error.toFixed(4)}%`);
}

console.log(`Status: ${iswAllPassed ? '✓ PASSED' : '✗ FAILED'} (Max error: ${iswMaxError.toFixed(4)}%)\n`);

// Summary
console.log('========================================');
console.log('Summary');
console.log('========================================');
console.log(`Harmonic Oscillator: ${hoAllPassed ? '✓ PASSED' : '✗ FAILED'}`);
console.log(`Infinite Square Well: ${iswAllPassed ? '✓ PASSED' : '✗ FAILED'}`);

if (hoAllPassed && iswAllPassed) {
  console.log('\n✓✓✓ ALL TESTS PASSED! ✓✓✓');
  console.log('The spectral method fix is WORKING CORRECTLY!');
  console.log('Removing the empirical correction factor solved the problem.');
} else {
  console.log('\n✗ Some tests failed.');
  console.log('The spectral method may need further investigation.');
}

console.log('========================================\n');
