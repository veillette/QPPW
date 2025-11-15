/**
 * Test SineDVR - the CORRECT DVR for finite domains with ψ=0 at boundaries
 * Formula from richford/dvr_py
 */

const HBAR = 1.054571817e-34;
const ELECTRON_MASS = 9.1093837015e-31;
const JOULES_TO_EV = 6.241509074e18;

function solveSineDVR(potential, mass, numStates, xMin, xMax, npts) {
  const L = xMax - xMin;
  const m = npts + 1;  // Key parameter for SineDVR

  // Generate grid points
  const dx = L / m;
  const xGrid = [];
  for (let i = 1; i <= npts; i++) {  // Note: starts at 1, not 0
    xGrid.push(xMin + i * dx);
  }

  console.log(`SineDVR with m = npts + 1 = ${m}`);
  console.log(`Grid spacing: ${(dx * 1e12).toFixed(2)} pm`);
  console.log(`First point: ${(xGrid[0] * 1e9).toFixed(4)} nm`);
  console.log(`Last point: ${(xGrid[npts-1] * 1e9).toFixed(4)} nm\n`);

  // Build kinetic energy matrix using SineDVR formula:
  // For i ≠ j:
  //   T_ij = (-1)^(i-j) · [1/sin²(π(i-j)/(2m)) - 1/sin²(π(i+j)/(2m))]
  // For i = j:
  //   T_ii = (2m² + 1)/3 - 1/sin²(πi/m)
  // Multiply by: π²/(2L²) · ℏ²/(2m_particle)

  const prefactor = (HBAR * HBAR) / (2 * mass) * (Math.PI * Math.PI) / (2 * L * L);
  const T = [];

  for (let i = 0; i < npts; i++) {
    T[i] = [];
    for (let j = 0; j < npts; j++) {
      const i_idx = i + 1;  // SineDVR uses 1-based indexing
      const j_idx = j + 1;

      if (i === j) {
        // Diagonal
        const sin_term = Math.sin((Math.PI * i_idx) / m);
        T[i][j] = prefactor * ((2 * m * m + 1) / 3 - 1 / (sin_term * sin_term));
      } else {
        // Off-diagonal
        const diff = i_idx - j_idx;
        const sum = i_idx + j_idx;
        const sign = Math.pow(-1, diff);

        const sin_diff = Math.sin((Math.PI * Math.abs(diff)) / (2 * m));
        const sin_sum = Math.sin((Math.PI * sum) / (2 * m));

        T[i][j] = prefactor * sign * (1 / (sin_diff * sin_diff) - 1 / (sin_sum * sin_sum));
      }
    }
  }

  // Build potential energy matrix
  const V = [];
  for (let i = 0; i < npts; i++) {
    V[i] = new Array(npts).fill(0);
    V[i][i] = potential(xGrid[i]);
  }

  // Hamiltonian
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
    energies: sorted.slice(0, numStates)
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

console.log('\n=== Testing SineDVR (Finite Domain with ψ=0 at boundaries) ===\n');

const mass = ELECTRON_MASS;
const omega = 1e15;
const springConstant = mass * omega * omega;
const numStates = 5;

const potential = (x) => 0.5 * springConstant * x * x;

const xMin = -5e-9;
const xMax = 5e-9;
const npts = 200;

const result = solveSineDVR(potential, mass, numStates, xMin, xMax, npts);

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

console.log('\n' + (allPassed ? '🎉 SUCCESS! SINEDVR WORKS! 🎉' : '✗ SineDVR also fails'));
process.exit(allPassed ? 0 : 1);
