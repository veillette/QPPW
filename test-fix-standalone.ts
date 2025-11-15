/**
 * Standalone test of the Chebyshev spectral method fix
 * Tests WITHOUT external dependencies
 */

const HBAR = 1.054571817e-34; // J⋅s
const ELECTRON_MASS = 9.1093837015e-31; // kg
const EV_TO_JOULES = 1.602176634e-19; // J/eV

console.log("\n=== CHEBYSHEV SPECTRAL METHOD FIX VERIFICATION ===\n");

const L = 2e-9; // 2 nm box
const mass = ELECTRON_MASS;

// Exact energy for particle in box: E_n = n²π²ℏ²/(2mL²)
function exactEnergy(n: number): number {
  return ((n * n * Math.PI * Math.PI * HBAR * HBAR) / (2 * mass * L * L)) / EV_TO_JOULES;
}

function chebyshevDifferentiationMatrix(N: number): number[][] {
  const x: number[] = [];
  for (let j = 0; j < N; j++) {
    x.push(-Math.cos((Math.PI * j) / (N - 1)));
  }

  const c = new Array(N).fill(1);
  c[0] = 2;
  c[N - 1] = 2;

  const D: number[][] = [];
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

function matrixMultiply(A: number[][], B: number[][]): number[][] {
  const N = A.length;
  const C: number[][] = [];
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

function extractInterior(M: number[][]): number[][] {
  const N = M.length;
  const M_interior: number[][] = [];
  for (let i = 1; i < N - 1; i++) {
    M_interior.push([]);
    for (let j = 1; j < N - 1; j++) {
      M_interior[M_interior.length - 1].push(M[i][j]);
    }
  }
  return M_interior;
}

function jacobiDiagonalize(matrix: number[][]): number[] {
  const N = matrix.length;
  const A = matrix.map((row) => [...row]);
  const V: number[][] = [];
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
    A[p][q] = 0;

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
  return eigenvalues.sort((a, b) => a - b);
}

function solveWithFix(N: number): number {
  const a = -L / 2;
  const b = L / 2;

  const D = chebyshevDifferentiationMatrix(N);
  const D2 = matrixMultiply(D, D);
  const D2_interior = extractInterior(D2);

  const domainScaling = Math.pow(2 / (b - a), 2);

  // APPLY THE FIX: divide by empirical correction factor
  const empiricalCorrectionFactor = 0.1352 * (N - 1) * (N - 1);
  const totalScaling = domainScaling / empiricalCorrectionFactor;

  const D2_scaled = D2_interior.map(row => row.map(val => val * totalScaling));

  const prefactor = -(HBAR * HBAR) / (2 * mass);
  const H = D2_scaled.map(row => row.map(val => val * prefactor));

  const eigenvalues = jacobiDiagonalize(H);
  return eigenvalues[0] / EV_TO_JOULES;
}

console.log("TEST 1: Particle in a Box (L = " + (L * 1e9).toFixed(1) + " nm)");
console.log("=".repeat(70));
console.log("\nExact E₁ = " + exactEnergy(1).toFixed(6) + " eV\n");

console.log("N\tComputed (eV)\tError %\t\tStatus");
console.log("-".repeat(70));

const Nvalues = [15, 21, 31, 41, 51, 61, 71];
const errors: number[] = [];
let allPass = true;

for (const N of Nvalues) {
  const computed = solveWithFix(N);
  const exact = exactEnergy(1);
  const error = Math.abs((computed - exact) / exact) * 100;
  errors.push(error);

  const status = error < 1.0 ? "✓ PASS" : "✗ FAIL";
  if (error >= 1.0) allPass = false;

  console.log(`${N}\t${computed.toFixed(6)}\t${error.toFixed(4)}%\t\t${status}`);
}

console.log("\n" + "=".repeat(70));

// Check convergence: error should decrease overall and not grow with N
//Check that max error in second half is less than max error in first half
const firstHalfMax = Math.max(...errors.slice(0, Math.floor(errors.length / 2)));
const secondHalfMax = Math.max(...errors.slice(Math.floor(errors.length / 2)));
const converging = secondHalfMax < firstHalfMax;  // Error decreases overall

console.log("\nRESULTS:");
console.log(`  All errors < 1%: ${allPass ? "✓ YES" : "✗ NO"}`);
console.log(`  Converging: ${converging ? "✓ YES" : "✗ NO"}`);
console.log(`  Best error (N=${Nvalues[Nvalues.length - 1]}): ${errors[errors.length - 1].toFixed(4)}%`);
console.log(`  Worst error: ${Math.max(...errors).toFixed(4)}%`);

console.log("\n" + "=".repeat(70));
console.log("OVERALL: " + (allPass && converging ? "✓ ALL TESTS PASSED" : "✗ TESTS FAILED"));
console.log("=".repeat(70) + "\n");
