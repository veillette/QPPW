/**
 * Test DVR with FIXED TOLERANCE (relative, not absolute)
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

  // FIX: Use RELATIVE tolerance based on matrix scale
  const matrixScale = Math.max(...matrix.map(row => Math.max(...row.map(Math.abs))));
  const tolerance = matrixScale * 1e-12;

  console.log(`Matrix scale: ${matrixScale.toExponential(4)}, Tolerance: ${tolerance.toExponential(4)}`);

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

    if (maxVal < tolerance) {
      console.log(`Converged in ${actualIterations} iterations`);
      break;
    }

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

console.log('\n=== DVR with FIXED TOLERANCE ===\n');

const mass = ELECTRON_MASS;
const omega = 1e15;
const springConstant = mass * omega * omega;
const numStates = 5;

const xMin = -5e-9;
const xMax = 5e-9;
const numPoints = 200;
const dx = (xMax - xMin) / (numPoints - 1);

const xGrid = [];
for (let i = 0; i < numPoints; i++) {
  xGrid.push(xMin + i * dx);
}

const V = [];
for (let i = 0; i < numPoints; i++) {
  V[i] = new Array(numPoints).fill(0);
  V[i][i] = 0.5 * springConstant * xGrid[i] * xGrid[i];
}

const T = createKineticEnergyMatrix(numPoints, dx, mass);

const H = [];
for (let i = 0; i < numPoints; i++) {
  H[i] = [];
  for (let j = 0; j < numPoints; j++) {
    H[i][j] = T[i][j] + V[i][j];
  }
}

const eigen = diagonalize(H);
const sorted = eigen.eigenvalues.sort((a, b) => a - b);

// Analytical
const E_analytical = [];
for (let n = 0; n < numStates; n++) {
  E_analytical.push(HBAR * omega * (n + 0.5));
}

console.log('\nEnergy Levels:');
console.log('State | Numerical (eV) | Analytical (eV) | Error (%)');
console.log('------|----------------|-----------------|----------');

let allPassed = true;
for (let n = 0; n < numStates; n++) {
  const E_num = sorted[n] * JOULES_TO_EV;
  const E_ana = E_analytical[n] * JOULES_TO_EV;
  const error = Math.abs((sorted[n] - E_analytical[n]) / E_analytical[n]) * 100;

  const passed = error < 1.0;
  allPassed = allPassed && passed;

  const status = passed ? '✓' : '✗';
  console.log(`E_${n}   ${status} | ${E_num.toFixed(6).padStart(13)} | ${E_ana.toFixed(6).padStart(15)} | ${error.toFixed(4).padStart(8)}`);
}

console.log('\n' + (allPassed ? '🎉🎉🎉 SUCCESS!!! DVR WORKS!!! 🎉🎉🎉' : '✗ Still failing...'));
process.exit(allPassed ? 0 : 1);
