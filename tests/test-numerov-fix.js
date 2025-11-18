#!/usr/bin/env node
/**
 * Standalone test for verifying Numerov solver fixes
 */

// Set up browser globals
globalThis.self = globalThis;
globalThis.window = globalThis;
globalThis.location = {
  href: 'http://localhost',
  search: '',
  hash: '',
  pathname: '/',
  host: 'localhost',
  hostname: 'localhost',
  port: '',
  protocol: 'http:',
  origin: 'http://localhost'
};

try {
  Object.defineProperty(globalThis, 'navigator', {
    value: { userAgent: 'node', platform: 'node', language: 'en-US' },
    writable: true,
    configurable: true
  });
} catch {}

globalThis.document = {
  createElement: () => ({ getContext: () => null, style: {}, setAttribute: () => {}, appendChild: () => {} }),
  body: { appendChild: () => {} },
  documentElement: { style: {} },
  createElementNS: () => ({ setAttribute: () => {} }),
  createTextNode: () => ({}),
  querySelector: () => null,
  querySelectorAll: () => []
};

globalThis.phet = globalThis.phet || {};
globalThis.phet.chipper = globalThis.phet.chipper || {};
globalThis.phet.joist = globalThis.phet.joist || {};
globalThis.phet.chipper.packageObject = { name: 'scenerystack' };
globalThis.phet.chipper.queryParameters = {};

globalThis.requestAnimationFrame = (cb) => setTimeout(cb, 16);
globalThis.cancelAnimationFrame = (id) => clearTimeout(id);

async function runTests() {
  console.log('╔════════════════════════════════════════════════════════════╗');
  console.log('║     Testing Numerov Solver Fixes                          ║');
  console.log('╚════════════════════════════════════════════════════════════╝\n');

  // Import modules
  const { solveMatrixNumerov } = await import('../src/common/model/MatrixNumerovSolver.js');
  const { integrateNumerovFromCenter, normalizeWavefunction } = await import('../src/common/model/NumerovSolver.js');
  const QuantumConstants = (await import('../src/common/model/QuantumConstants.js')).default;

  const HBAR = QuantumConstants.HBAR;
  const ELECTRON_MASS = QuantumConstants.ELECTRON_MASS;
  const EV_TO_JOULES = QuantumConstants.EV_TO_JOULES;

  let passed = 0;
  let failed = 0;

  // Test 1: MatrixNumerov - Harmonic Oscillator
  console.log('=== Test 1: MatrixNumerov - Harmonic Oscillator ===');
  {
    const mass = ELECTRON_MASS;
    const omega = 1e15;
    const springConstant = mass * omega * omega;

    const gridConfig = {
      xMin: -5e-9,
      xMax: 5e-9,
      numPoints: 128
    };

    const potential = (x) => 0.5 * springConstant * x * x;

    const result = solveMatrixNumerov(potential, mass, 5, gridConfig);

    // Analytical energies: E_n = ℏω(n + 1/2)
    const analytical = [];
    for (let n = 0; n < 5; n++) {
      analytical.push(HBAR * omega * (n + 0.5));
    }

    // Debug: check boundary potential
    const V_boundary = potential(gridConfig.xMin);
    console.log(`V_boundary = ${V_boundary / EV_TO_JOULES} eV`);
    console.log(`Expected E_0 = ${(HBAR * omega * 0.5) / EV_TO_JOULES} eV`);

    console.log(`Found ${result.energies.length} energy levels:`);
    let allPass = true;
    if (result.energies.length === 0) {
      console.log('  No energy levels found!');
      allPass = false;
    }
    for (let n = 0; n < Math.min(5, result.energies.length); n++) {
      const numEV = result.energies[n] / EV_TO_JOULES;
      const anaEV = analytical[n] / EV_TO_JOULES;
      const error = Math.abs((result.energies[n] - analytical[n]) / analytical[n]) * 100;
      const pass = error < 1.0; // 1% tolerance
      if (!pass) allPass = false;
      console.log(`  E_${n}: Numerical=${numEV.toFixed(6)} eV, Analytical=${anaEV.toFixed(6)} eV, Error=${error.toFixed(4)}% ${pass ? '✓' : '✗'}`);
    }

    if (allPass) {
      console.log('Result: PASSED ✓\n');
      passed++;
    } else {
      console.log('Result: FAILED ✗\n');
      failed++;
    }
  }

  // Test 2: MatrixNumerov - Finite Square Well
  console.log('=== Test 2: MatrixNumerov - Finite Square Well ===');
  {
    const mass = ELECTRON_MASS;
    const wellWidth = 1e-9;
    const wellDepth = 5 * EV_TO_JOULES;

    const gridConfig = {
      xMin: -2e-9,
      xMax: 2e-9,
      numPoints: 128
    };

    const potential = (x) => {
      const halfWidth = wellWidth / 2;
      return (x >= -halfWidth && x <= halfWidth) ? 0 : wellDepth;
    };

    const result = solveMatrixNumerov(potential, mass, 5, gridConfig);

    // Check that we get positive, increasing energies
    console.log(`Found ${result.energies.length} energy levels:`);
    let allPass = true;
    let prevE = -Infinity;
    for (let n = 0; n < Math.min(5, result.energies.length); n++) {
      const eEV = result.energies[n] / EV_TO_JOULES;
      const isIncreasing = result.energies[n] > prevE;
      const isBound = result.energies[n] < wellDepth;
      const pass = isIncreasing && isBound;
      if (!pass) allPass = false;
      console.log(`  E_${n}: ${eEV.toFixed(6)} eV (bound=${isBound}, increasing=${isIncreasing}) ${pass ? '✓' : '✗'}`);
      prevE = result.energies[n];
    }

    if (result.energies.length === 0) {
      console.log('Result: FAILED ✗ (No bound states found)\n');
      failed++;
    } else if (allPass) {
      console.log('Result: PASSED ✓ (Energies are positive, increasing, and bound)\n');
      passed++;
    } else {
      console.log('Result: FAILED ✗\n');
      failed++;
    }
  }

  // Test 3: integrateNumerovFromCenter - Symmetric state
  console.log('=== Test 3: NumerovSolver - Symmetric State Integration ===');
  {
    const mass = ELECTRON_MASS;
    const wellDepth = 5 * EV_TO_JOULES;
    const wellWidth = 1e-9;

    const numPoints = 201;
    const xMin = -2e-9;
    const xMax = 2e-9;
    const dx = (xMax - xMin) / (numPoints - 1);

    const xGrid = [];
    for (let i = 0; i < numPoints; i++) {
      xGrid.push(xMin + i * dx);
    }

    const V = xGrid.map(x => {
      const halfWidth = wellWidth / 2;
      return (x >= -halfWidth && x <= halfWidth) ? 0 : wellDepth;
    });

    // Test with ground state energy estimate
    const E = 0.5 * EV_TO_JOULES;

    const psi = integrateNumerovFromCenter(E, V, xGrid, dx, mass, 'symmetric');
    const normalizedPsi = normalizeWavefunction(psi, dx);

    // Check symmetry: ψ(-x) = ψ(x)
    const centerIdx = Math.floor(numPoints / 2);
    let symmetryError = 0;
    const checkPoints = 20;
    for (let i = 1; i <= checkPoints; i++) {
      const leftIdx = centerIdx - i;
      const rightIdx = centerIdx + i;
      if (leftIdx >= 0 && rightIdx < numPoints) {
        symmetryError += Math.abs(normalizedPsi[leftIdx] - normalizedPsi[rightIdx]);
      }
    }
    symmetryError /= checkPoints;

    const pass = symmetryError < 0.01; // Small symmetry error
    console.log(`  Symmetry error: ${symmetryError.toFixed(6)}`);
    console.log(`  Center value: ${normalizedPsi[centerIdx].toFixed(6)}`);

    if (pass) {
      console.log('Result: PASSED ✓ (Wavefunction is symmetric)\n');
      passed++;
    } else {
      console.log('Result: FAILED ✗\n');
      failed++;
    }
  }

  // Test 4: integrateNumerovFromCenter - Antisymmetric state
  console.log('=== Test 4: NumerovSolver - Antisymmetric State Integration ===');
  {
    const mass = ELECTRON_MASS;
    const wellDepth = 5 * EV_TO_JOULES;
    const wellWidth = 1e-9;

    const numPoints = 201;
    const xMin = -2e-9;
    const xMax = 2e-9;
    const dx = (xMax - xMin) / (numPoints - 1);

    const xGrid = [];
    for (let i = 0; i < numPoints; i++) {
      xGrid.push(xMin + i * dx);
    }

    const V = xGrid.map(x => {
      const halfWidth = wellWidth / 2;
      return (x >= -halfWidth && x <= halfWidth) ? 0 : wellDepth;
    });

    // Test with first excited state energy estimate
    const E = 1.5 * EV_TO_JOULES;

    const psi = integrateNumerovFromCenter(E, V, xGrid, dx, mass, 'antisymmetric');
    const normalizedPsi = normalizeWavefunction(psi, dx);

    // Check antisymmetry: ψ(-x) = -ψ(x)
    const centerIdx = Math.floor(numPoints / 2);
    let antisymmetryError = 0;
    const checkPoints = 20;
    for (let i = 1; i <= checkPoints; i++) {
      const leftIdx = centerIdx - i;
      const rightIdx = centerIdx + i;
      if (leftIdx >= 0 && rightIdx < numPoints) {
        antisymmetryError += Math.abs(normalizedPsi[leftIdx] + normalizedPsi[rightIdx]);
      }
    }
    antisymmetryError /= checkPoints;

    const centerNearZero = Math.abs(normalizedPsi[centerIdx]) < 0.01;
    const pass = antisymmetryError < 0.01 && centerNearZero;

    console.log(`  Antisymmetry error: ${antisymmetryError.toFixed(6)}`);
    console.log(`  Center value: ${normalizedPsi[centerIdx].toFixed(6)}`);

    if (pass) {
      console.log('Result: PASSED ✓ (Wavefunction is antisymmetric)\n');
      passed++;
    } else {
      console.log('Result: FAILED ✗\n');
      failed++;
    }
  }

  // Summary
  console.log('╔════════════════════════════════════════════════════════════╗');
  console.log('║     Test Summary                                          ║');
  console.log('╚════════════════════════════════════════════════════════════╝');
  console.log(`Total: ${passed + failed}, Passed: ${passed}, Failed: ${failed}`);

  if (failed === 0) {
    console.log('\n✓ All tests passed! Numerov solver fixes are working correctly.');
  } else {
    console.log('\n✗ Some tests failed. Please review the results above.');
    process.exit(1);
  }
}

runTests().catch(err => {
  console.error('Test error:', err);
  process.exit(1);
});
