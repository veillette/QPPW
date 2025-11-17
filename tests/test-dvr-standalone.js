/**
 * Standalone DVR test script
 * This directly implements the test without relying on bundled code
 */

// Constants
const HBAR = 1.054571817e-34; // Planck's constant / 2π (J·s)
const ELECTRON_MASS = 9.1093837015e-31; // kg
const JOULES_TO_EV = 6.241509074e18; // Conversion factor

/**
 * Create kinetic energy matrix using Colbert-Miller formula
 */
function createKineticEnergyMatrix(N, dx, mass) {
  const prefactor = (HBAR * HBAR) / (2 * mass * dx * dx);

  const T = [];
  for (let i = 0; i < N; i++) {
    T[i] = [];
    for (let j = 0; j < N; j++) {
      if (i === j) {
        // Diagonal elements
        T[i][j] = prefactor * (Math.PI * Math.PI) / 3;
      } else {
        // Off-diagonal elements
        const diff = i - j;
        const sign = Math.pow(-1, diff);
        T[i][j] = prefactor * (2 * sign) / (diff * diff);
      }
    }
  }

  return T;
}

/**
 * Diagonalize a symmetric matrix using Jacobi eigenvalue algorithm
 */
function diagonalize(matrix) {
  const N = matrix.length;

  // Copy matrix (we'll modify it)
  const A = matrix.map(row => [...row]);

  // Initialize eigenvectors as identity matrix
  const V = [];
  for (let i = 0; i < N; i++) {
    V[i] = new Array(N).fill(0);
    V[i][i] = 1;
  }

  // Jacobi iteration
  const maxIterations = 50 * N * N;
  const tolerance = 1e-12;

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
    const App = c * c * A[p][p] - 2 * s * c * A[p][q] + s * s * A[q][q];
    const Aqq = s * s * A[p][p] + 2 * s * c * A[p][q] + c * c * A[q][q];
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
 * Solve using DVR method
 */
function solveDVR(potential, mass, numStates, gridConfig) {
  const { xMin, xMax, numPoints } = gridConfig;
  const dx = (xMax - xMin) / (numPoints - 1);

  // Generate grid
  const xGrid = [];
  for (let i = 0; i < numPoints; i++) {
    xGrid.push(xMin + i * dx);
  }

  const N = numPoints;

  // Construct potential energy matrix V (diagonal)
  const potentialEnergyMatrix = [];
  for (let i = 0; i < N; i++) {
    potentialEnergyMatrix[i] = new Array(N).fill(0);
    potentialEnergyMatrix[i][i] = potential(xGrid[i]);
  }

  // Construct kinetic energy matrix T
  const kineticEnergyMatrix = createKineticEnergyMatrix(N, dx, mass);

  // Construct Hamiltonian H = T + V
  const hamiltonianMatrix = [];
  for (let i = 0; i < N; i++) {
    hamiltonianMatrix[i] = [];
    for (let j = 0; j < N; j++) {
      hamiltonianMatrix[i][j] = kineticEnergyMatrix[i][j] + potentialEnergyMatrix[i][j];
    }
  }

  // Diagonalize
  const eigen = diagonalize(hamiltonianMatrix);

  // Sort by energy
  const sortedIndices = eigen.eigenvalues
    .map((e, i) => ({ energy: e, index: i }))
    .sort((a, b) => a.energy - b.energy)
    .map(item => item.index);

  // Extract bound states
  const energies = [];
  const V_boundary = Math.max(potential(xMin), potential(xMax));

  for (let i = 0; i < Math.min(numStates, N); i++) {
    const idx = sortedIndices[i];
    const energy = eigen.eigenvalues[idx];

    if (energy < V_boundary) {
      energies.push(energy);
    }
  }

  return { energies, xGrid };
}

/**
 * Analytical harmonic oscillator solution
 */
function solveHarmonicOscillator(springConstant, mass, numStates) {
  const omega = Math.sqrt(springConstant / mass);
  const energies = [];

  for (let n = 0; n < numStates; n++) {
    // E_n = ℏω(n + 1/2)
    energies.push(HBAR * omega * (n + 0.5));
  }

  return { energies };
}

/**
 * Test DVR with harmonic oscillator
 */
function testDVRHarmonicOscillator() {
  console.log('\n=== DVR - Harmonic Oscillator Test ===\n');

  const mass = ELECTRON_MASS;
  const omega = 1e15;
  const springConstant = mass * omega * omega;
  const numStates = 5;

  const gridConfig = {
    xMin: -5e-9,
    xMax: 5e-9,
    numPoints: 200
  };

  const potential = (x) => 0.5 * springConstant * x * x;

  // Get numerical solution
  const numerical = solveDVR(potential, mass, numStates, gridConfig);

  // Get analytical solution
  const analytical = solveHarmonicOscillator(springConstant, mass, numStates);

  // Compare
  console.log('State | Numerical (eV) | Analytical (eV) | Error (%)');
  console.log('------|----------------|-----------------|----------');

  let allPassed = true;
  let maxError = 0;

  for (let n = 0; n < numStates; n++) {
    const E_num = numerical.energies[n] * JOULES_TO_EV;
    const E_ana = analytical.energies[n] * JOULES_TO_EV;
    const error = Math.abs((numerical.energies[n] - analytical.energies[n]) / analytical.energies[n]) * 100;

    const passed = error < 1.0;
    allPassed = allPassed && passed;
    maxError = Math.max(maxError, error);

    const status = passed ? '✓' : '✗';
    console.log(`E_${n}   ${status} | ${E_num.toFixed(6).padStart(13)} | ${E_ana.toFixed(6).padStart(15)} | ${error.toFixed(4).padStart(8)}`);
  }

  console.log('\nMaximum error: ' + maxError.toFixed(4) + '%');
  console.log('Status: ' + (allPassed ? '✓ PASSED' : '✗ FAILED'));

  return allPassed;
}

// Run the test
console.log('========================================');
console.log('DVR Method Accuracy Test');
console.log('========================================');
console.log('Tolerance: 1% error from analytical solution\n');

const passed = testDVRHarmonicOscillator();

console.log('\n========================================');
if (passed) {
  console.log('✓ DVR test PASSED');
} else {
  console.log('✗ DVR test FAILED - Debug needed');
}
console.log('========================================\n');

process.exit(passed ? 0 : 1);
