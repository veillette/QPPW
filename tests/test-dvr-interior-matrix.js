/**
 * Test DVR with interior matrix extraction (like spectral method)
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

function extractInteriorMatrix(H) {
  const N = H.length;
  const H_interior = [];

  for (let i = 1; i < N - 1; i++) {
    H_interior[i - 1] = [];
    for (let j = 1; j < N - 1; j++) {
      H_interior[i - 1].push(H[i][j]);
    }
  }

  return H_interior;
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

function testDVR_InteriorMatrix() {
  console.log('\n=== DVR with Interior Matrix (Like Spectral Method) ===\n');

  const mass = ELECTRON_MASS;
  const omega = 1e15;
  const springConstant = mass * omega * omega;
  const numStates = 5;

  const xMin = -5e-9;
  const xMax = 5e-9;
  const numPoints = 200;

  // Grid INCLUDING boundaries (current implementation)
  const dx = (xMax - xMin) / (numPoints - 1);

  const xGrid = [];
  for (let i = 0; i < numPoints; i++) {
    xGrid.push(xMin + i * dx);
  }

  console.log(`Grid spacing: dx = ${dx.toExponential(4)} m`);
  console.log(`Total points (with boundaries): ${numPoints}`);
  console.log(`Interior points (without boundaries): ${numPoints - 2}`);

  // Build full matrices
  const V = [];
  for (let i = 0; i < numPoints; i++) {
    V[i] = new Array(numPoints).fill(0);
    V[i][i] = 0.5 * springConstant * xGrid[i] * xGrid[i];
  }

  const T = createKineticEnergyMatrix(numPoints, dx, mass);

  const H_full = [];
  for (let i = 0; i < numPoints; i++) {
    H_full[i] = [];
    for (let j = 0; j < numPoints; j++) {
      H_full[i][j] = T[i][j] + V[i][j];
    }
  }

  // Extract interior matrix (remove boundary points)
  const H = extractInteriorMatrix(H_full);
  console.log(`Hamiltonian size after boundary removal: ${H.length} x ${H.length}`);

  const eigen = diagonalize(H);

  const sorted = eigen.eigenvalues
    .map((e, i) => ({ energy: e, index: i }))
    .sort((a, b) => a.energy - b.energy);

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
    const E_num = sorted[n].energy * JOULES_TO_EV;
    const E_ana = E_analytical[n] * JOULES_TO_EV;
    const error = Math.abs((sorted[n].energy - E_analytical[n]) / E_analytical[n]) * 100;

    const passed = error < 1.0;
    allPassed = allPassed && passed;

    const status = passed ? '✓' : '✗';
    console.log(`E_${n}   ${status} | ${E_num.toFixed(6).padStart(13)} | ${E_ana.toFixed(6).padStart(15)} | ${error.toFixed(4).padStart(8)}`);
  }

  console.log('\n' + (allPassed ? '✓✓✓ ALL TESTS PASSED! ✓✓✓' : '✗ Some tests failed'));
  return allPassed;
}

const result = testDVR_InteriorMatrix();
process.exit(result ? 0 : 1);
