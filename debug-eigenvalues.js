/**
 * Debug: Print ALL eigenvalues to see what's happening
 */

const HBAR = 1.054571817e-34;
const ELECTRON_MASS = 9.1093837015e-31;
const JOULES_TO_EV = 6.241509074e18;

function createKineticEnergyMatrix(N, dx, mass) {
  const prefactor = (HBAR * HBAR) / (2 * mass * dx * dx);
  const T = [];
  for (let i = 0; i < N; i++) {
    T[i] = [];
    for (let j = 0; j < N; j++) {
      if (i === j) {
        T[i][j] = prefactor * (Math.PI * Math.PI) / 3;
      } else {
        const diff = i - j;
        const sign = Math.pow(-1, diff);
        T[i][j] = prefactor * (2 * sign) / (diff * diff);
      }
    }
  }
  return T;
}

function diagonalize(matrix) {
  const N = matrix.length;
  const A = matrix.map(row => [...row]);
  const V = [];
  for (let i = 0; i < N; i++) {
    V[i] = new Array(N).fill(0);
    V[i][i] = 1;
  }

  const maxIterations = 50 * N * N;
  const tolerance = 1e-12;
  let actualIterations = 0;

  for (let iter = 0; iter < maxIterations; iter++) {
    actualIterations = iter + 1;
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
  return { eigenvalues, actualIterations };
}

console.log('\n=== Eigenvalue Debug ===\n');

const mass = ELECTRON_MASS;
const omega = 1e15;
const springConstant = mass * omega * omega;

const xMin = -5e-9;
const xMax = 5e-9;
const numPoints = 50;  // Use fewer points for debugging

const dx = (xMax - xMin) / (numPoints - 1);

const xGrid = [];
for (let i = 0; i < numPoints; i++) {
  xGrid.push(xMin + i * dx);
}

console.log(`N = ${numPoints} points\n`);

// Test 1: Just kinetic energy (V=0)
console.log('TEST 1: Pure kinetic energy (V=0)');
const T = createKineticEnergyMatrix(numPoints, dx, mass);
const eigen_T = diagonalize(T);
const sorted_T = eigen_T.eigenvalues.sort((a, b) => a - b);

console.log(`Diagonalization took ${eigen_T.actualIterations} iterations`);
console.log(`Lowest 5 eigenvalues (eV):`);
for (let i = 0; i < 5; i++) {
  console.log(`  λ_${i}: ${(sorted_T[i] * JOULES_TO_EV).toFixed(6)}`);
}
console.log(`Highest 5 eigenvalues (eV):`);
for (let i = numPoints - 5; i < numPoints; i++) {
  console.log(`  λ_${i}: ${(sorted_T[i] * JOULES_TO_EV).toFixed(6)}`);
}

// Test 2: Kinetic + Potential
console.log('\nTEST 2: Kinetic + Harmonic Potential');
const V = [];
for (let i = 0; i < numPoints; i++) {
  V[i] = new Array(numPoints).fill(0);
  V[i][i] = 0.5 * springConstant * xGrid[i] * xGrid[i];
}

const H = [];
for (let i = 0; i < numPoints; i++) {
  H[i] = [];
  for (let j = 0; j < numPoints; j++) {
    H[i][j] = T[i][j] + V[i][j];
  }
}

const eigen_H = diagonalize(H);
const sorted_H = eigen_H.eigenvalues.sort((a, b) => a - b);

console.log(`Diagonalization took ${eigen_H.actualIterations} iterations`);
console.log(`Lowest 5 eigenvalues (eV):`);
for (let i = 0; i < 5; i++) {
  console.log(`  E_${i}: ${(sorted_H[i] * JOULES_TO_EV).toFixed(6)}`);
}

console.log('\n Expected E_0 = 0.329106 eV');

// Check matrix diagonal
console.log('\nSample potential values at grid points (eV):');
for (let i = 0; i < Math.min(5, numPoints); i++) {
  console.log(`  V(x[${i}] = ${(xGrid[i] * 1e9).toFixed(3)} nm) = ${(V[i][i] * JOULES_TO_EV).toFixed(3)}`);
}
