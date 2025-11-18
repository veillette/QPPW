#!/usr/bin/env node
/**
 * Debug test for MatrixNumerov solver
 */

// Set up browser globals
globalThis.self = globalThis;
globalThis.window = globalThis;
globalThis.location = {
  href: 'http://localhost', search: '', hash: '', pathname: '/',
  host: 'localhost', hostname: 'localhost', port: '', protocol: 'http:', origin: 'http://localhost'
};
try {
  Object.defineProperty(globalThis, 'navigator', {
    value: { userAgent: 'node', platform: 'node', language: 'en-US' },
    writable: true, configurable: true
  });
} catch {
  // Ignore error if navigator is already defined
}
globalThis.document = {
  createElement: () => ({ getContext: () => null, style: {}, setAttribute: () => {}, appendChild: () => {} }),
  body: { appendChild: () => {} }, documentElement: { style: {} },
  createElementNS: () => ({ setAttribute: () => {} }), createTextNode: () => ({}),
  querySelector: () => null, querySelectorAll: () => []
};
globalThis.phet = globalThis.phet || {};
globalThis.phet.chipper = globalThis.phet.chipper || { packageObject: { name: 'scenerystack' }, queryParameters: {} };
globalThis.phet.joist = globalThis.phet.joist || {};
globalThis.requestAnimationFrame = (cb) => setTimeout(cb, 16);
globalThis.cancelAnimationFrame = (id) => clearTimeout(id);

async function debug() {
  console.log('=== Debug: MatrixNumerov Eigenvalue Check ===\n');

  const { DotMatrix, diagonalize, matrixToArray } = await import('../src/common/model/LinearAlgebraUtils.js');
  const QuantumConstants = (await import('../src/common/model/QuantumConstants.js')).default;

  const HBAR = QuantumConstants.HBAR;
  const ELECTRON_MASS = QuantumConstants.ELECTRON_MASS;
  const EV_TO_JOULES = QuantumConstants.EV_TO_JOULES;

  // Small test case - harmonic oscillator
  const mass = ELECTRON_MASS;
  const omega = 1e15;
  const springConstant = mass * omega * omega;

  const numPoints = 32; // Small for debugging
  const xMin = -5e-9;
  const xMax = 5e-9;
  const dx = (xMax - xMin) / (numPoints - 1);

  // Generate grid
  const xGrid = [];
  for (let i = 0; i < numPoints; i++) {
    xGrid.push(xMin + i * dx);
  }

  // Potential
  const V = xGrid.map(x => 0.5 * springConstant * x * x);

  const N = numPoints;
  const coeff = (2 * mass) / (HBAR * HBAR);
  const h2 = dx * dx;
  // FIXED: Use separate factors for H and S matrices
  const hFactor = (h2 / 12) * coeff;  // For H matrix potential terms
  const sFactor = coeff / 12;          // For S matrix overlap terms (no h²!)

  // Construct matrices
  const H = new DotMatrix(N, N);
  const S = new DotMatrix(N, N);

  for (let i = 0; i < N; i++) {
    H.set(i, i, (2 / h2) - (10 * hFactor * V[i]));
    S.set(i, i, 10 * sFactor);

    if (i > 0) {
      const avgV = (V[i - 1] + V[i]) / 2;
      const offDiagValue = -(1 / h2) - (hFactor * avgV);
      H.set(i, i - 1, offDiagValue);
      H.set(i - 1, i, offDiagValue);
      S.set(i, i - 1, sFactor);
      S.set(i - 1, i, sFactor);
    }
  }

  console.log(`Grid: ${N} points, dx = ${dx.toExponential(3)} m`);
  console.log(`V_boundary = ${(V[0] / EV_TO_JOULES).toFixed(3)} eV`);
  console.log(`Expected E_0 = ${(HBAR * omega * 0.5 / EV_TO_JOULES).toFixed(6)} eV`);
  console.log('');

  // Check some matrix values
  console.log('H matrix samples:');
  console.log(`  H[0,0] = ${H.get(0, 0).toExponential(3)}`);
  console.log(`  H[N/2,N/2] = ${H.get(Math.floor(N/2), Math.floor(N/2)).toExponential(3)}`);
  console.log(`  H[1,0] = ${H.get(1, 0).toExponential(3)}`);
  console.log('');

  console.log('S matrix samples:');
  console.log(`  S[0,0] = ${S.get(0, 0).toExponential(3)}`);
  console.log(`  S[1,0] = ${S.get(1, 0).toExponential(3)}`);
  console.log('');

  // Solve generalized eigenvalue problem
  try {
    const Sinv = S.inverse();
    const SinvH = Sinv.times(H);

    console.log('Sinv*H matrix samples:');
    console.log(`  (S^-1*H)[0,0] = ${SinvH.get(0, 0).toExponential(3)}`);
    console.log(`  (S^-1*H)[N/2,N/2] = ${SinvH.get(Math.floor(N/2), Math.floor(N/2)).toExponential(3)}`);
    console.log('');

    const SinvHArray = matrixToArray(SinvH);
    const eigen = diagonalize(SinvHArray);

    // Sort eigenvalues
    const sortedEigs = [...eigen.eigenvalues].sort((a, b) => a - b);

    console.log('First 10 eigenvalues:');
    for (let i = 0; i < Math.min(10, sortedEigs.length); i++) {
      const eEV = sortedEigs[i] / EV_TO_JOULES;
      console.log(`  λ_${i} = ${sortedEigs[i].toExponential(4)} J = ${eEV.toFixed(6)} eV`);
    }
    console.log('');

    // Check boundary condition
    const V_boundary = Math.max(V[0], V[N - 1]);
    const boundStates = sortedEigs.filter(e => e < V_boundary && isFinite(e));
    console.log(`Bound states (E < ${(V_boundary/EV_TO_JOULES).toFixed(2)} eV): ${boundStates.length}`);

  } catch (err) {
    console.error('Error:', err.message);
  }
}

debug().catch(err => console.error('Debug error:', err));
