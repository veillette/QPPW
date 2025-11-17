/**
 * Exact implementation based on richford/dvr_py SincDVR
 * Using their exact grid and formula
 */

const HBAR = 1.054571817e-34;
const ELECTRON_MASS = 9.1093837015e-31;
const JOULES_TO_EV = 6.241509074e18;

function solveDVR_ReferenceMethod(potential, mass, numStates, x0, L, npts) {
  // Grid spacing: a = L / npts (NOT L/(npts-1))
  const a = L / npts;

  // Generate grid: x_i = x0 + i*a - L/2 + a/2
  // This creates cell-centered points that DON'T include boundaries
  const xGrid = [];
  for (let i = 0; i < npts; i++) {
    const x = x0 + i * a - L / 2 + a / 2;
    xGrid.push(x);
  }

  console.log(`Grid spacing a = ${(a * 1e12).toFixed(2)} pm`);
  console.log(`First point: ${(xGrid[0] * 1e9).toFixed(4)} nm`);
  console.log(`Last point: ${(xGrid[npts-1] * 1e9).toFixed(4)} nm`);
  console.log(`(Boundaries would be at ±${(L/2 * 1e9).toFixed(4)} nm)\n`);

  // Build kinetic energy matrix using reference formula:
  // T_ii = π²/3/a²
  // T_ij = 2·(-1)^(i-j)/(i-j)²/a²
  // Multiplied by ℏ²/(2m)

  const prefactor = (HBAR * HBAR) / (2 * mass);
  const T = [];

  for (let i = 0; i < npts; i++) {
    T[i] = [];
    for (let j = 0; j < npts; j++) {
      if (i === j) {
        // Diagonal: π²/3/a²
        T[i][j] = prefactor * (Math.PI * Math.PI) / (3 * a * a);
      } else {
        // Off-diagonal: 2·(-1)^(i-j)/(i-j)²/a²
        const diff = i - j;
        const sign = Math.pow(-1, diff);
        T[i][j] = prefactor * 2 * sign / (diff * diff * a * a);
      }
    }
  }

  // Build potential energy matrix (diagonal)
  const V = [];
  for (let i = 0; i < npts; i++) {
    V[i] = new Array(npts).fill(0);
    V[i][i] = potential(xGrid[i]);
  }

  // Hamiltonian H = T + V
  const H = [];
  for (let i = 0; i < npts; i++) {
    H[i] = [];
    for (let j = 0; j < npts; j++) {
      H[i][j] = T[i][j] + V[i][j];
    }
  }

  // Diagonalize
  const eigen = diagonalize(H);
  const sorted = eigen.eigenvalues.sort((a, b) => a - b);

  return {
    energies: sorted.slice(0, numStates),
    xGrid
  };
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

console.log('\n=== DVR using EXACT Reference Implementation ===\n');

const mass = ELECTRON_MASS;
const omega = 1e15;
const springConstant = mass * omega * omega;
const numStates = 5;

// Harmonic oscillator potential
const potential = (x) => 0.5 * springConstant * x * x;

// Grid parameters
const x0 = 0;  // Center
const L = 10e-9;  // Total length (same as original test: -5nm to +5nm)
const npts = 200;

const result = solveDVR_ReferenceMethod(potential, mass, numStates, x0, L, npts);

// Analytical
const E_analytical = [];
for (let n = 0; n < numStates; n++) {
  E_analytical.push(HBAR * omega * (n + 0.5));
}

console.log('Energy Levels:');
console.log('State | Numerical (eV) | Analytical (eV) | Error (%)');
console.log('------|----------------|-----------------|----------');

let allPassed = true;
for (let n = 0; n < numStates; n++) {
  const E_num = result.energies[n] * JOULES_TO_EV;
  const E_ana = E_analytical[n] * JOULES_TO_EV;
  const error = Math.abs((result.energies[n] - E_analytical[n]) / E_analytical[n]) * 100;

  const passed = error < 1.0;
  allPassed = allPassed && passed;

  const status = passed ? '✓' : '✗';
  console.log(`E_${n}   ${status} | ${E_num.toFixed(6).padStart(13)} | ${E_ana.toFixed(6).padStart(15)} | ${error.toFixed(4).padStart(8)}`);
}

console.log('\n' + (allPassed ? '✓✓✓ PERFECT! ✓✓✓' : '✗ Still not working'));
process.exit(allPassed ? 0 : 1);
