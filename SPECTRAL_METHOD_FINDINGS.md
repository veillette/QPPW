# Spectral Method Investigation Findings

## Summary

The Chebyshev spectral method implementation has **fundamental issues** with excited states that were not caught in the original testing.

## Test Results

### What Works ✓
- **Infinite Square Well (V=0) - Ground State Only**: <1% error for N=15 to N=250
- Empirical correction factor formula: `0.1352 * (N-1)²` is verified for this specific case

### What Doesn't Work ✗
- **All Excited States**: 70-95% error for n≥2 (even for infinite square well)
- **Harmonic Oscillator - All States**: 76-99% error for ground state and all excited states
- **Any Non-Zero Potential**: The empirical correction factor only works for V=0
- **Spurious Degeneracies**: All excited states appear in degenerate pairs (incorrect for 1D)

### Critical Limitation

**The empirical correction factor ONLY works for potentials with V=0 (infinite square well).**
It does NOT work for:
- Harmonic oscillator (V = ½kx²)
- Any other potential with V≠0

This makes the spectral method essentially unusable for realistic quantum mechanical problems.

## Root Cause

1. **Original Testing Was Incomplete**: The commit (0bc976d) that added the empirical correction factor **only tested ground states**:
   ```typescript
   console.log("Testing ground state (n=1) energy\n");
   ```
   No verification of excited states was performed.

2. **Empirical Correction Only Works for Ground State**: The factor `0.1352*(N-1)²` correctly scales the ground state eigenvalue, but does NOT work for higher eigenvalues.

3. **Fundamental Scaling Issue**: Different eigenvalues are being scaled by different amounts, but we're applying a single global correction factor to the entire second derivative matrix.

## Evidence

### Infinite Square Well Test (N=31)
```
Spectral energies:
  E_0: 0.376573 eV  ← ground state: 0.14% error ✓
  E_1: 0.380747 eV  ← should be 1.504 eV (74% error) ✗
  E_2: 0.380747 eV  ← should be 3.384 eV (88% error) ✗
  E_3: 0.393646 eV  ← should be 6.016 eV (93% error) ✗
  E_4: 0.393646 eV  ← should be 9.401 eV (95% error) ✗

Analytical energies:
  E_1: 0.376030 eV
  E_2: 1.504121 eV  (4× E_1)
  E_3: 3.384272 eV  (9× E_1)
  E_4: 6.016483 eV  (16× E_1)
  E_5: 9.400754 eV  (25× E_1)
```

Note the spurious 2-fold degeneracies in the spectral results (E_1=E_2, E_3=E_4), which should not occur in 1D.

## Impact on AccuracyTests

The `AccuracyTests.ts` file was added in commit 02ced19 and includes tests for:
- DVR - Harmonic Oscillator ← **FAILING** (separate DVR bug)
- Spectral - Harmonic Oscillator ← **FAILING** (95-99% error all states)
- DVR - Infinite Square Well ← **FAILING** (separate DVR bug)
- Spectral - Infinite Square Well ← **FAILING** (only ground state works)

**None of these tests were actually passing when committed.**

## Bugs Fixed in This Session

1. **Bound State Filter Bug**: The code had:
   ```typescript
   if (energy < V_boundary)
   ```
   For V=0 potentials, this rejected all positive energy states. Fixed by removing this filter for spectral method with Dirichlet BCs.

2. **Updated Documentation**: Improved comments explaining the empirical correction factor and its verified range.

## Recommendations

### Short Term (This PR)
1. **Update AccuracyTests** to:
   - Mark spectral tests as "known issues" or skip them
   - Only test ground states for spectral method
   - Focus on DVR method accuracy

2. **Add Clear Documentation** warning users that:
   - Spectral method only reliable for ground states
   - Excited states have significant errors
   - Method should be considered experimental

3. **Add User Task Tests**: Test that spectral method works for different numbers of points **for ground states only**.

### Long Term (Future Work)
1. **Research Eigenvalue-Dependent Scaling**: Investigate why different eigenvalues scale differently
2. **Consider Alternative Approaches**:
   - Different boundary condition handling
   - Different spectral basis (Fourier, Legendre, etc.)
   - Eigenvalue-specific correction factors
3. **Fix or Remove Spectral Method**: Either fix properly or mark as not ready for production use

## Conclusion

The spectral method implementation is **not production-ready**. The empirical correction factor approach is fundamentally flawed because it applies a single global scale factor when different eigenvalues require different corrections.

The ground state accuracy (<1% error) is excellent, but this doesn't extend to excited states, making the method unsuitable for most quantum mechanical applications which require multiple energy levels.
