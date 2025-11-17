/**
 * Standalone test for FGH accuracy
 * Run with: node test-fgh-accuracy.js
 */

// Quantum constants
const HBAR = 1.054571817e-34; // Planck's constant / 2π (J·s)
const ELECTRON_MASS = 9.1093837015e-31; // Electron mass (kg)

/**
 * Complex number operations
 */
function complexAdd(a, b) {
  return { real: a.real + b.real, imag: a.imag + b.imag };
}

function complexSubtract(a, b) {
  return { real: a.real - b.real, imag: a.imag - b.imag };
}

function complexMultiply(a, b) {
  return {
    real: a.real * b.real - a.imag * b.imag,
    imag: a.real * b.imag + a.imag * b.real,
  };
}

/**
 * FFT implementation
 */
function fft(x) {
  const N = x.length;
  if (N <= 1) return x;

  const even = [];
  const odd = [];
  for (let i = 0; i < N; i++) {
    if (i % 2 === 0) even.push(x[i]);
    else odd.push(x[i]);
  }

  const fftEven = fft(even);
  const fftOdd = fft(odd);

  const result = new Array(N);
  for (let k = 0; k < N / 2; k++) {
    const angle = (-2 * Math.PI * k) / N;
    const twiddle = { real: Math.cos(angle), imag: Math.sin(angle) };
    const temp = complexMultiply(twiddle, fftOdd[k]);
    result[k] = complexAdd(fftEven[k], temp);
    result[k + N / 2] = complexSubtract(fftEven[k], temp);
  }
  return result;
}

/**
 * IFFT implementation
 */
function ifft(X) {
  const N = X.length;
  const XConj = X.map((x) => ({ real: x.real, imag: -x.imag }));
  const result = fft(XConj);
  return result.map((x) => ({ real: x.real / N, imag: -x.imag / N }));
}

/**
 * FFT frequency array
 */
function fftFreq(N, dx) {
  const waveNumberArray = [];
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
 * Build FGH Hamiltonian matrix
 */
function buildFGHHamiltonian(N, T_k, V_x) {
  const H = [];
  for (let i = 0; i < N; i++) {
    H[i] = new Array(N).fill(0);
  }

  for (let j = 0; j < N; j++) {
    const psi_x = [];
    for (let i = 0; i < N; i++) {
      psi_x.push({ real: i === j ? 1.0 : 0.0, imag: 0.0 });
    }

    const psi_k = fft(psi_x);
    for (let i = 0; i < N; i++) {
      psi_k[i].real *= T_k[i];
      psi_k[i].imag *= T_k[i];
    }

    const T_psi_x = ifft(psi_k);
    for (let i = 0; i < N; i++) {
      H[i][j] = T_psi_x[i].real;
    }
  }

  for (let i = 0; i < N; i++) {
    H[i][i] += V_x[i];
  }

  return H;
}

/**
 * Jacobi diagonalization
 */
function diagonalize(matrix) {
  const N = matrix.length;
  const A = matrix.map((row) => [...row]);
  const V = [];
  for (let i = 0; i < N; i++) {
    V[i] = new Array(N).fill(0);
    V[i][i] = 1;
  }

  const maxIterations = 50 * N * N;
  const matrixScale = Math.max(...matrix.map((row) => Math.max(...row.map(Math.abs))));
  const tolerance = matrixScale * 1e-12;

  for (let iter = 0; iter < maxIterations; iter++) {
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

    if (maxVal < tolerance) break;

    const theta = 0.5 * Math.atan2(2 * A[rowIndex][colIndex], A[colIndex][colIndex] - A[rowIndex][rowIndex]);
    const cosTheta = Math.cos(theta);
    const sinTheta = Math.sin(theta);

    const App = cosTheta * cosTheta * A[rowIndex][rowIndex] - 2 * sinTheta * cosTheta * A[rowIndex][colIndex] + sinTheta * sinTheta * A[colIndex][colIndex];
    const Aqq = sinTheta * sinTheta * A[rowIndex][rowIndex] + 2 * sinTheta * cosTheta * A[rowIndex][colIndex] + cosTheta * cosTheta * A[colIndex][colIndex];

    A[rowIndex][rowIndex] = App;
    A[colIndex][colIndex] = Aqq;
    A[rowIndex][colIndex] = 0;
    A[colIndex][rowIndex] = 0;

    for (let i = 0; i < N; i++) {
      if (i !== rowIndex && i !== colIndex) {
        const Aip = cosTheta * A[i][rowIndex] - sinTheta * A[i][colIndex];
        const Aiq = sinTheta * A[i][rowIndex] + cosTheta * A[i][colIndex];
        A[i][rowIndex] = Aip;
        A[rowIndex][i] = Aip;
        A[i][colIndex] = Aiq;
        A[colIndex][i] = Aiq;
      }
    }

    for (let i = 0; i < N; i++) {
      const Vip = cosTheta * V[i][rowIndex] - sinTheta * V[i][colIndex];
      const Viq = sinTheta * V[i][rowIndex] + cosTheta * V[i][colIndex];
      V[i][rowIndex] = Vip;
      V[i][colIndex] = Viq;
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
 * Normalize wavefunction
 */
function normalize(psi, dx) {
  let integral = 0;
  for (let i = 0; i < psi.length - 1; i++) {
    integral += (psi[i] * psi[i] + psi[i + 1] * psi[i + 1]) / 2;
  }
  integral *= dx;
  const normalization = Math.sqrt(integral);
  return psi.map((val) => val / normalization);
}

/**
 * FGH solver
 */
function solveFGH(potential, mass, numStates, gridConfig) {
  const { xMin, xMax, numPoints } = gridConfig;
  const dx = (xMax - xMin) / numPoints;
  const N = numPoints;

  const xGrid = [];
  for (let i = 0; i < N; i++) {
    xGrid.push(xMin + i * dx);
  }

  const waveNumberGrid = fftFreq(N, dx);
  const T_k = waveNumberGrid.map((k) => (HBAR * HBAR * k * k) / (2 * mass));
  const V_x = xGrid.map(potential);

  const H = buildFGHHamiltonian(N, T_k, V_x);
  const eigen = diagonalize(H);

  const sortedIndices = eigen.eigenvalues
    .map((e, i) => ({ energy: e, index: i }))
    .sort((a, b) => a.energy - b.energy)
    .map((item) => item.index);

  const energies = [];
  const wavefunctions = [];

  for (let i = 0; i < Math.min(numStates, N); i++) {
    const idx = sortedIndices[i];
    const energy = eigen.eigenvalues[idx];

    const V_boundary = Math.max(potential(xMin), potential(xMax));
    if (energy < V_boundary) {
      energies.push(energy);
      const wavefunction = eigen.eigenvectors[idx];
      const normalizedPsi = normalize(wavefunction, dx);
      wavefunctions.push(normalizedPsi);
    }
  }

  return { energies, wavefunctions, xGrid };
}

/**
 * Test FGH with harmonic oscillator
 */
function testFGHHarmonicOscillator() {
  console.log('\n=== FGH - Harmonic Oscillator Test ===\n');

  const mass = ELECTRON_MASS;
  const omega = 1e15; // rad/s
  const springConstant = mass * omega * omega;
  const numStates = 5;

  const gridConfig = {
    xMin: -5e-9,
    xMax: 5e-9,
    numPoints: 256, // Must be power of 2 for radix-2 FFT
  };

  const potential = (x) => 0.5 * springConstant * x * x;

  // Get FGH solution
  const result = solveFGH(potential, mass, numStates, gridConfig);

  // Analytical energies: E_n = ℏω(n + 1/2)
  console.log('Testing energy levels:');
  let maxError = 0;
  let allPassed = true;

  for (let n = 0; n < numStates; n++) {
    const E_numerical = result.energies[n];
    const E_analytical = HBAR * omega * (n + 0.5);
    const error = Math.abs((E_numerical - E_analytical) / E_analytical) * 100;

    const E_numerical_eV = E_numerical / 1.602176634e-19;
    const E_analytical_eV = E_analytical / 1.602176634e-19;

    const passed = error < 1.0;
    allPassed = allPassed && passed;
    maxError = Math.max(maxError, error);

    const status = passed ? '✓ PASS' : '✗ FAIL';
    console.log(`  E_${n}: ${status} - Numerical: ${E_numerical_eV.toFixed(6)} eV, Analytical: ${E_analytical_eV.toFixed(6)} eV, Error: ${error.toFixed(4)}%`);
  }

  console.log(`\nMaximum error: ${maxError.toFixed(4)}%`);
  if (allPassed) {
    console.log('✓ Test PASSED! FGH method is within 1% accuracy.\n');
  } else {
    console.log('✗ Test FAILED! FGH method exceeds 1% error.\n');
  }

  return { passed: allPassed, maxError };
}

// Run the test
console.log('========================================');
console.log('FGH Solver Accuracy Test');
console.log('========================================');

const result = testFGHHarmonicOscillator();

console.log('========================================');
if (result.passed) {
  console.log('✓ FGH solver passes accuracy test');
} else {
  console.log('✗ FGH solver fails accuracy test');
  console.log(`Maximum error: ${result.maxError.toFixed(4)}%`);
}
console.log('========================================\n');
