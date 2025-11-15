/**
 * Check the actual kinetic energy matrix elements
 */

const HBAR = 1.054571817e-34;
const ELECTRON_MASS = 9.1093837015e-31;
const JOULES_TO_EV = 6.241509074e18;

const mass = ELECTRON_MASS;
const xMin = -5e-9;
const xMax = 5e-9;
const numPoints = 10;  // Small for easy inspection

const dx = (xMax - xMin) / (numPoints - 1);
const prefactor = (HBAR * HBAR) / (2 * mass * dx * dx);

console.log('\n=== Kinetic Energy Matrix Elements ===\n');
console.log(`dx = ${(dx * 1e9).toFixed(4)} nm`);
console.log(`prefactor = ℏ²/(2m·dx²) = ${(prefactor * JOULES_TO_EV).toExponential(6)} eV\n`);

// Build matrix
const T = [];
for (let i = 0; i < numPoints; i++) {
  T[i] = [];
  for (let j = 0; j < numPoints; j++) {
    if (i === j) {
      T[i][j] = prefactor * (Math.PI * Math.PI) / 3;
    } else {
      const diff = i - j;
      const sign = Math.pow(-1, diff);
      T[i][j] = prefactor * (2 * sign) / (diff * diff);
    }
  }
}

// Print diagonal elements
console.log('Diagonal elements T_ii = prefactor · π²/3:');
for (let i = 0; i < Math.min(5, numPoints); i++) {
  console.log(`  T[${i}][${i}] = ${(T[i][i] * JOULES_TO_EV).toExponential(6)} eV`);
}

// Print some off-diagonal elements
console.log('\nOff-diagonal elements T_ij = prefactor · 2·(-1)^(i-j)/(i-j)²:');
console.log(`  T[0][1] (i=0, j=1, diff=-1) = ${(T[0][1] * JOULES_TO_EV).toExponential(6)} eV`);
console.log(`  T[0][2] (i=0, j=2, diff=-2) = ${(T[0][2] * JOULES_TO_EV).toExponential(6)} eV`);
console.log(`  T[1][0] (i=1, j=0, diff=+1) = ${(T[1][0] * JOULES_TO_EV).toExponential(6)} eV`);
console.log(`  T[1][2] (i=1, j=2, diff=-1) = ${(T[1][2] * JOULES_TO_EV).toExponential(6)} eV`);

// Check ratio of off-diagonal to diagonal
const ratio_01 = Math.abs(T[0][1] / T[0][0]);
const ratio_02 = Math.abs(T[0][2] / T[0][0]);
console.log(`\nRatio |T[0][1]| / T[0][0] = ${ratio_01.toFixed(4)}`);
console.log(`Ratio |T[0][2]| / T[0][0] = ${ratio_02.toFixed(4)}`);

// Print full matrix (in eV, scientific notation)
console.log('\nFull T matrix (eV, first 5x5):');
for (let i = 0; i < Math.min(5, numPoints); i++) {
  let row = '  ';
  for (let j = 0; j < Math.min(5, numPoints); j++) {
    const val = T[i][j] * JOULES_TO_EV;
    row += val.toExponential(3).padStart(12) + '  ';
  }
  console.log(row);
}

// Check if matrix is symmetric
let isSymmetric = true;
for (let i = 0; i < numPoints; i++) {
  for (let j = i + 1; j < numPoints; j++) {
    if (Math.abs(T[i][j] - T[j][i]) > 1e-15) {
      isSymmetric = false;
      console.log(`\nNOT SYMMETRIC: T[${i}][${j}] = ${T[i][j]}, T[${j}][${i}] = ${T[j][i]}`);
    }
  }
}
if (isSymmetric) {
  console.log('\n✓ Matrix is symmetric');
}
