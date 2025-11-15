/**
 * Test if the eigenvalue error depends on grid size N
 * This will tell us if it's a systematic bug or a convergence issue
 */

// Constants (inlined to avoid dependencies)
const HBAR = 1.054571817e-34; // J⋅s
const ELECTRON_MASS = 9.1093837015e-31; // kg
const EV_TO_JOULES = 1.602176634e-19; // J/eV

console.log("\n=== Grid Convergence Test ===\n");

const L = 2e-9;
const mass = ELECTRON_MASS;
const E1_expected = (Math.PI * Math.PI * HBAR * HBAR) / (2 * mass * L * L);

console.log(`Box: [-L/2, L/2] with L = ${(L*1e9).toFixed(1)} nm`);
console.log(`Expected E₁ = ${(E1_expected/EV_TO_JOULES).toFixed(6)} eV\n`);

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
  return eigenvalues.sort((a, b) => a - b);
}

function testForN(N: number): number {
  const a = -L / 2;
  const b = L / 2;

  const D = chebyshevDifferentiationMatrix(N);
  const D2 = matrixMultiply(D, D);

  const D2_interior: number[][] = [];
  for (let i = 1; i < N - 1; i++) {
    D2_interior.push([]);
    for (let j = 1; j < N - 1; j++) {
      D2_interior[D2_interior.length - 1].push(D2[i][j]);
    }
  }

  const scaling = Math.pow(2 / (b - a), 2);
  const D2_interior_scaled = D2_interior.map(row => row.map(val => val * scaling));

  const prefactor = -(HBAR * HBAR) / (2 * mass);
  const H = D2_interior_scaled.map(row => row.map(val => val * prefactor));

  const eigenvalues = jacobiDiagonalize(H);
  const E1_computed = eigenvalues[0];

  return E1_computed / EV_TO_JOULES;
}

console.log("N\tInterior Points\tE₁ (eV)\t\tError Factor");
console.log("=".repeat(60));

for (const N of [11, 15, 21, 31, 41, 51]) {
  const E1_computed = testForN(N);
  const errorFactor = E1_computed / (E1_expected / EV_TO_JOULES);
  console.log(`${N}\t${N-2}\t\t${E1_computed.toFixed(6)}\t${errorFactor.toFixed(2)}x`);
}

console.log("\nExpected E₁ = " + (E1_expected/EV_TO_JOULES).toFixed(6) + " eV");
