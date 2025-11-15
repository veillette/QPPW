# Spectral Method Implementation Audit Report

**Date**: 2025-11-15
**Auditor**: Claude
**Files Reviewed**:
- `src/common/model/SpectralSolver.ts`
- `src/common/model/AccuracyTests.ts`
- `src/common/model/analytical-solutions/*.ts`

## Executive Summary

The Chebyshev spectral method implementation has **3 confirmed bugs** and **1 fundamental limitation**:

1. ✗ **CRITICAL**: Incorrect method name in return value
2. ✗ **CRITICAL**: Missing Jacobian in wavefunction normalization
3. ⚠️ **LIMITATION**: Empirical correction factor only works for ground states

## Detailed Findings

### Bug #1: Incorrect Method Name
**Severity**: Low (metadata only)
**Location**: `SpectralSolver.ts:132`
**Status**: CONFIRMED BUG

```typescript
// Current (WRONG):
return {
  energies,
  wavefunctions,
  xGrid,
  method: "dvr", // ← Should be "spectral"
};
```

**Impact**: Misleading metadata that could cause confusion in debugging or result interpretation.

**Fix**: Change `method: "dvr"` to `method: "spectral"`

---

### Bug #2: Missing Jacobian in Wavefunction Normalization
**Severity**: High (affects wavefunction accuracy)
**Location**: `SpectralSolver.ts:232-253` (normalizeChebyshev function)
**Status**: CONFIRMED BUG

**Problem**:
The normalization computes ∫|ψ(ξ)|²dξ = 1 over the computational domain ξ ∈ [-1,1], but fails to account for the coordinate transformation Jacobian when ensuring normalization in physical space.

**Mathematical Analysis**:
```
Coordinate transformation: x = ((xMax - xMin)·ξ + (xMax + xMin))/2
Jacobian: dx = (xMax - xMin)/2 · dξ

Required normalization: ∫|ψ(x)|² dx = 1
Which means: ∫|ψ(ξ)|² · (xMax - xMin)/2 · dξ = 1

Current code normalizes: ∫|ψ(ξ)|² dξ = 1 ← WRONG!
```

**Current Code** (lines 245-252):
```typescript
// Compute ∫|ψ|² dξ
let integral = 0;
for (let i = 0; i < N; i++) {
  integral += weights[i] * psi[i] * psi[i];
}

const normalization = Math.sqrt(integral);
return psi.map((val) => val / normalization);
```

**Impact**:
- Wavefunctions are incorrectly normalized by a factor of √((xMax-xMin)/2)
- For wellWidth = 1nm, this is a factor of √(0.5e-9) = 2.236e-5
- This doesn't affect energy eigenvalues but makes probability densities incorrect

**Fix**:
```typescript
// Compute ∫|ψ|² dξ and include Jacobian for physical space
const jacobian = (xMax - xMin) / 2; // dx/dξ scaling factor
let integral = 0;
for (let i = 0; i < N; i++) {
  integral += weights[i] * psi[i] * psi[i] * jacobian;
}
```

---

### Bug #3 (Known Issue): Empirical Correction Factor Limitation
**Severity**: Critical (makes method unusable for most applications)
**Location**: `SpectralSolver.ts:55-73`
**Status**: DOCUMENTED LIMITATION

**Problem**:
The empirical correction factor `0.1352 * (N-1)²` was calibrated ONLY for:
- Infinite square well (V=0)
- Ground state eigenvalue (n=1)

It does NOT work for:
- Excited states (n≥2) → 75-95% errors
- Non-zero potentials (harmonic oscillator) → 76-99% errors
- Any realistic quantum system requiring multiple energy levels

**Test Results from Previous Investigation**:
```
Infinite Square Well (N=31):
  Ground state:     0.14% error ✓ PASS
  1st excited:     75.0% error ✗ FAIL
  2nd excited:     88.9% error ✗ FAIL

Harmonic Oscillator (N=100):
  Ground state:    93.4% error ✗ FAIL
  All states:      76-99% error ✗ FAIL
```

**Root Cause**:
Applying a single global scale factor to the entire second derivative matrix affects all eigenvalues equally, but different energy levels require different corrections. This is a fundamental architectural flaw.

**Fix**:
This cannot be fixed with a simple correction. Options include:
1. Remove the empirical correction and use a theoretically sound approach
2. Implement eigenvalue-dependent corrections (complex)
3. Use a different spectral basis or boundary condition method
4. Mark the method as experimental/deprecated

---

## Mathematical Verification

### Chebyshev Differentiation Matrix ✓ CORRECT
Reviewed lines 149-184. The implementation correctly follows Boyd's formulas:
- ✓ Chebyshev-Gauss-Lobatto points: x_j = -cos(πj/(N-1))
- ✓ Weights: c_0 = c_{N-1} = 2, others = 1
- ✓ Diagonal elements: D_jj = -x_j/(2(1-x_j²)) for interior points
- ✓ Boundary diagonals: D_00 = (2(N-1)² + 1)/6, D_{N-1,N-1} = -(2(N-1)² + 1)/6
- ✓ Off-diagonal: D_ij = (c_i/c_j)·(-1)^{i+j}/(x_i - x_j)

### Matrix Operations ✓ CORRECT
- ✓ matrixMultiply (lines 207-223): Standard matrix multiplication
- ✓ extractInteriorMatrix (lines 190-202): Correctly removes boundaries for Dirichlet BCs
- ✓ diagonalize (lines 262-352): Jacobi iteration for symmetric matrices

### Domain Scaling ✓ CORRECT
- ✓ Line 52: domainScalingFactor = 2/(xMax - xMin)
- ✓ Line 53: Second derivative scaling: (domainScalingFactor)²

### Hamiltonian Construction ✓ CORRECT
- ✓ Lines 77-91: H = T + V with correct kinetic prefactor -(ℏ²/2m)
- ✓ Potential added on diagonal before extracting interior (fine since boundaries are removed)

---

## Test Expectations vs. Reality

### Test: Spectral - Infinite Square Well
**Expected** (AccuracyTests.ts:250-312):
- 5 energy levels with <1% error each
- Ground state E₁ = π²ℏ²/(2mL²)
- Excited states E_n = n²E₁

**Reality**:
- ✓ Ground state: <1% error
- ✗ Excited states: 75-95% error
- ✗ Spurious 2-fold degeneracies (incorrect for 1D)

### Test: Spectral - Harmonic Oscillator
**Expected** (AccuracyTests.ts:113-175):
- 5 energy levels with <1% error each
- Energies E_n = ℏω(n + 1/2)

**Reality**:
- ✗ All states: 76-99% error
- ✗ Empirical correction incompatible with V≠0

---

## Recommendations

### Immediate Actions (This PR)
1. ✅ **Fix Bug #1**: Change method name to "spectral"
2. ✅ **Fix Bug #2**: Add Jacobian to wavefunction normalization
3. ⚠️ **Document Bug #3**: Add warnings that method only works for ground states

### Medium Term
1. **Update Tests**: Modify AccuracyTests to only test ground states for spectral method
2. **Add User Warnings**: Display warning in UI when spectral method selected for excited states
3. **Consider Deprecation**: Mark spectral method as experimental

### Long Term (Future Research)
1. **Investigate Root Cause**: Why does D² have this (N-1)² scaling issue?
2. **Theoretical Solution**: Find proper way to handle Chebyshev derivatives with Dirichlet BCs
3. **Alternative Methods**: Consider Fourier spectral or other bases
4. **Eigenvalue-Dependent Correction**: Research if per-eigenvalue corrections are feasible

---

## Conclusion

The spectral method implementation has **2 fixable bugs** (method name, normalization) and **1 fundamental limitation** (empirical correction).

After fixing the two bugs, the method will:
- ✓ Correctly identify itself as "spectral"
- ✓ Produce properly normalized wavefunctions
- ⚠️ Still only work reliably for ground states of V=0 potentials

**The spectral method is NOT production-ready for general quantum mechanical problems requiring multiple energy levels.**
