/**
 * Debug the second derivative matrix to understand the infinite square well issue
 */

function chebyshevPoints(N) {
  const points = [];
  for (let j = 0; j < N; j++) {
    points.push(-Math.cos((Math.PI * j) / (N - 1)));
  }
  return points;
}

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

// Small test case
const N = 11;
const xMin = -0.5e-9;
const xMax = 0.5e-9;
const wellWidth = 1e-9;

console.log('========================================');
console.log('Debugging Spectral Method Second Derivative');
console.log(`N = ${N} points`);
console.log(`Domain: [${xMin}, ${xMax}]`);
console.log('========================================\n');

// Build the second derivative matrix
const D = chebyshevDifferentiationMatrix(N);
const D2 = matrixMultiply(D, D);

// Apply scaling
const scale = 2 / (xMax - xMin);
console.log(`Domain scaling factor: ${scale.toExponential(4)}`);
console.log(`Scale squared: ${(scale * scale).toExponential(4)}\n`);

for (let i = 0; i < N; i++) {
  for (let j = 0; j < N; j++) {
    D2[i][j] *= scale * scale;
  }
}

// Extract interior
const D2_interior = extractInterior(D2);

console.log('Interior second derivative matrix (first 5x5):');
for (let i = 0; i < Math.min(5, D2_interior.length); i++) {
  const row = D2_interior[i].slice(0, 5).map(v => v.toExponential(2)).join('  ');
  console.log(row);
}
console.log();

// Check diagonal
console.log('Diagonal elements of interior D²:');
const diag = D2_interior.map((row, i) => row[i]);
console.log(diag.map(v => v.toExponential(4)).join(', '));
console.log();

// Check if matrix is symmetric
let symmetric = true;
for (let i = 0; i < D2_interior.length && symmetric; i++) {
  for (let j = i + 1; j < D2_interior.length && symmetric; j++) {
    if (Math.abs(D2_interior[i][j] - D2_interior[j][i]) > 1e-10) {
      symmetric = false;
      console.log(`Not symmetric at [${i},${j}]: ${D2_interior[i][j]} vs ${D2_interior[j][i]}`);
    }
  }
}
console.log(`Matrix is ${symmetric ? '' : 'NOT '}symmetric\n`);

// Test analytical formula for infinite square well
console.log('Analytical eigenvalues for particle in a box:');
const HBAR = 1.054571817e-34;
const ELECTRON_MASS = 9.1093837015e-31;
const joulesToEV = (j) => j / 1.602176634e-19;

for (let n = 1; n <= 5; n++) {
  const E_n = (n * n * Math.PI * Math.PI * HBAR * HBAR) / (2 * ELECTRON_MASS * wellWidth * wellWidth);
  console.log(`E_${n} = ${joulesToEV(E_n).toFixed(6)} eV`);
}
console.log();

// The eigenvalue of d²/dx² should be -(nπ/L)²
console.log('Expected eigenvalues of d²/dx² operator (should be negative):');
for (let n = 1; n <= 5; n++) {
  const lambda = -((n * Math.PI / wellWidth) ** 2);
  console.log(`λ_${n} = ${lambda.toExponential(4)}`);
}
console.log();

// When multiplied by -(ℏ²/2m), we get positive energy
const kineticPrefactor = -(HBAR * HBAR) / (2 * ELECTRON_MASS);
console.log(`Kinetic prefactor: -(ℏ²/2m) = ${kineticPrefactor.toExponential(4)}`);
console.log();

console.log('Expected Hamiltonian eigenvalues (should be positive):');
for (let n = 1; n <= 5; n++) {
  const lambda_d2 = -((n * Math.PI / wellWidth) ** 2);
  const E_n = kineticPrefactor * lambda_d2;
  console.log(`E_${n} = ${kineticPrefactor.toExponential(4)} * ${lambda_d2.toExponential(4)} = ${joulesToEV(E_n).toFixed(6)} eV`);
}

console.log('\n========================================');
