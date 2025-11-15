/**
 * Test various DVR formula variations to find the correct one
 */

const HBAR = 1.054571817e-34;
const ELECTRON_MASS = 9.1093837015e-31;
const JOULES_TO_EV = 6.241509074e18;

function createKineticEnergyMatrix_Variant(N, dx, mass, variant) {
  const T = [];
  let prefactor;
  let diag_factor;
  let offdiag_numerator;

  switch(variant) {
    case 1: // Original
      prefactor = (HBAR * HBAR) / (2 * mass * dx * dx);
      diag_factor = (Math.PI * Math.PI) / 3;
      offdiag_numerator = 2;
      break;
    case 2: // Try 2π instead of π
      prefactor = (HBAR * HBAR) / (2 * mass * dx * dx);
      diag_factor = (4 * Math.PI * Math.PI) / 3;
      offdiag_numerator = 2;
      break;
    case 3: // Try without the dx² (maybe it's already in the formula)
      prefactor = (HBAR * HBAR) / (2 * mass);
      diag_factor = (Math.PI * Math.PI) / (3 * dx * dx);
      offdiag_numerator = 2 / (dx * dx);
      break;
    case 4: // Try dx instead of dx²
      prefactor = (HBAR * HBAR) / (2 * mass * dx);
      diag_factor = (Math.PI * Math.PI) / 3;
      offdiag_numerator = 2;
      break;
    case 5: // Try multiplying by N or N²
      prefactor = (HBAR * HBAR) / (2 * mass * dx * dx) * N * N;
      diag_factor = (Math.PI * Math.PI) / 3;
      offdiag_numerator = 2;
      break;
    case 6: // Different normalization
      prefactor = (HBAR * HBAR) / (2 * mass * dx * dx) / (N * N);
      diag_factor = (Math.PI * Math.PI) / 3;
      offdiag_numerator = 2;
      break;
    default:
      prefactor = (HBAR * HBAR) / (2 * mass * dx * dx);
      diag_factor = (Math.PI * Math.PI) / 3;
      offdiag_numerator = 2;
  }

  for (let i = 0; i < N; i++) {
    T[i] = [];
    for (let j = 0; j < N; j++) {
      if (i === j) {
        T[i][j] = prefactor * diag_factor;
      } else {
        const diff = i - j;
        const sign = Math.pow(-1, diff);
        T[i][j] = prefactor * (offdiag_numerator * sign) / (diff * diff);
      }
    }
  }
  return T;
}

function diagonalizeSimple(matrix) {
  const N = matrix.length;
  const A = matrix.map(row => [...row]);
  const V = [];
  for (let i = 0; i < N; i++) {
    V[i] = new Array(N).fill(0);
    V[i][i] = 1;
  }

  const maxIterations = 50 * N * N;
  const tolerance = 1e-12;

  for (let iter = 0; iter < maxIterations; iter++) {
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
  return { eigenvalues };
}

console.log('\n=== Testing Formula Variations ===\n');

const mass = ELECTRON_MASS;
const omega = 1e15;
const springConstant = mass * omega * omega;

const xMin = -5e-9;
const xMax = 5e-9;
const numPoints = 200;
const dx = (xMax - xMin) / (numPoints - 1);

const xGrid = [];
for (let i = 0; i < numPoints; i++) {
  xGrid.push(xMin + i * dx);
}

// Analytical ground state
const E_0_analytical = (HBAR * omega * 0.5) * JOULES_TO_EV;

console.log(`Analytical ground state: ${E_0_analytical.toFixed(6)} eV\n`);

for (let variant = 1; variant <= 6; variant++) {
  const V = [];
  for (let i = 0; i < numPoints; i++) {
    V[i] = new Array(numPoints).fill(0);
    V[i][i] = 0.5 * springConstant * xGrid[i] * xGrid[i];
  }

  const T = createKineticEnergyMatrix_Variant(numPoints, dx, mass, variant);

  const H = [];
  for (let i = 0; i < numPoints; i++) {
    H[i] = [];
    for (let j = 0; j < numPoints; j++) {
      H[i][j] = T[i][j] + V[i][j];
    }
  }

  const eigen = diagonalizeSimple(H);
  const sorted = eigen.eigenvalues.sort((a, b) => a - b);
  const E_0_numerical = sorted[0] * JOULES_TO_EV;
  const error = Math.abs((E_0_numerical - E_0_analytical) / E_0_analytical) * 100;

  console.log(`Variant ${variant}: E_0 = ${E_0_numerical.toFixed(6)} eV, Error = ${error.toFixed(2)}%`);
}
