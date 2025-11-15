/**
 * Test DVR on infinite square well - simpler than harmonic oscillator
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

console.log('\n=== DVR: Infinite Square Well Test ===\n');

const mass = ELECTRON_MASS;
const wellWidth = 1e-9; // 1 nm well
const numStates = 5;

// Grid: from -wellWidth to +wellWidth (well is centered, width is wellWidth)
const xMin = -wellWidth;
const xMax = wellWidth;
const numPoints = 150;

const dx = (xMax - xMin) / (numPoints - 1);

const xGrid = [];
for (let i = 0; i < numPoints; i++) {
  xGrid.push(xMin + i * dx);
}

console.log(`Well width: ${(wellWidth * 1e9).toFixed(1)} nm`);
console.log(`Grid: ${(xMin * 1e9).toFixed(2)} to ${(xMax * 1e9).toFixed(2)} nm`);
console.log(`Grid spacing: ${(dx * 1e12).toFixed(2)} pm`);
console.log(`Number of points: ${numPoints}\n`);

// Potential: 0 inside well, very high outside
const V_high = 1e10; // Very high but not infinite
const halfWidth = wellWidth / 2;
const V = [];
for (let i = 0; i < numPoints; i++) {
  V[i] = new Array(numPoints).fill(0);
  const x = xGrid[i];
  V[i][i] = (x >= -halfWidth && x <= halfWidth) ? 0 : V_high;
}

// Kinetic energy
const T = createKineticEnergyMatrix(numPoints, dx, mass);

// Hamiltonian
const H = [];
for (let i = 0; i < numPoints; i++) {
  H[i] = [];
  for (let j = 0; j < numPoints; j++) {
    H[i][j] = T[i][j] + V[i][j];
  }
}

const eigen = diagonalize(H);
const sorted = eigen.eigenvalues.sort((a, b) => a - b);

// Analytical solution for infinite square well: E_n = n²π²ℏ²/(2m·L²)
const L = wellWidth; // well width
const E_analytical = [];
for (let n = 1; n <= numStates; n++) {
  E_analytical.push((n * n * Math.PI * Math.PI * HBAR * HBAR) / (2 * mass * L * L));
}

console.log('Energy Levels:');
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
  console.log(`E_${n + 1}   ${status} | ${E_num.toFixed(6).padStart(13)} | ${E_ana.toFixed(6).padStart(15)} | ${error.toFixed(4).padStart(8)}`);
}

console.log('\n' + (allPassed ? '✓✓✓ WORKS FOR SQUARE WELL! ✓✓✓' : '✗ Also fails for square well'));
process.exit(allPassed ? 0 : 1);
