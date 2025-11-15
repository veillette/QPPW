# Spectral Method Testing Results

## User Request
Test that the spectral method works for different numbers of grid points.

## Answer: PARTIALLY WORKS

### What Works ✓
The spectral method **does work correctly** for different numbers of grid points **for the specific case of**:
- **Infinite Square Well (V=0 potential)**
- **Ground state only**

#### Test Results
```
N=15:   0.94% error  ✓ PASS
N=25:   0.27% error  ✓ PASS
N=50:   0.11% error  ✓ PASS
N=75:   0.04% error  ✓ PASS
N=100:  0.03% error  ✓ PASS
N=150:  0.06% error  ✓ PASS
N=200:  0.07% error  ✓ PASS
N=250:  0.07% error  ✓ PASS
```

All values are within the 1% error tolerance and show good convergence.

### What Doesn't Work ✗

1. **Harmonic Oscillator** - Fails for ALL grid point values:
   ```
   N=50:   76.7% error  ✗ FAIL
   N=100:  93.4% error  ✗ FAIL
   N=200:  97.5% error  ✗ FAIL
   ```

2. **Excited States** - Even for infinite square well, excited states fail:
   ```
   Ground state (n=1): 0.06% error  ✓ PASS
   1st excited (n=2):  75.0% error  ✗ FAIL
   2nd excited (n=3):  88.9% error  ✗ FAIL
   ```

## Root Cause

The empirical correction factor `0.1352 * (N-1)²` was calibrated ONLY for:
- V=0 potentials (infinite square well)
- Ground state eigenvalues

It does NOT generalize to:
- Non-zero potentials (like harmonic oscillator)
- Excited state eigenvalues

## Changes Made

### 1. Fixed Bound State Filter Bug
**Before:**
```typescript
const V_boundary = Math.max(potential(xMin), potential(xMax));
if (energy < V_boundary) {  // Wrong for V=0!
  energies.push(energy);
}
```

**After:**
```typescript
// For spectral method with Dirichlet BCs, all eigenvalues are bound states
for (let i = 0; i < Math.min(numStates, H_interior.length); i++) {
  energies.push(energy);
}
```

### 2. Improved Documentation
Added detailed comments explaining:
- The empirical correction factor and its verified range
- Physical interpretation of the scaling
- Limitations of the current approach

## Files Modified

1. **src/common/model/SpectralSolver.ts** - Fixed bound state filter, improved documentation
2. **SPECTRAL_METHOD_FINDINGS.md** - Comprehensive investigation report
3. **SPECTRAL_METHOD_STATUS.md** - This summary

## Recommendations

### For Users
- ✓ **Safe to use**: Spectral method for infinite square well ground states
- ⚠️ **Avoid**: Using spectral method for:
  - Harmonic oscillator or any V≠0 potential
  - Excited states
  - Production applications requiring multiple energy levels

### For Developers
- The spectral method needs fundamental rework to handle:
  - Non-zero potentials
  - Excited state eigenvalues
- Consider marking as experimental or disabling until properly fixed
- Alternative: Focus on DVR method which may work better (needs separate investigation)

## Original AccuracyTests Status

The `AccuracyTests.ts` file added in commit 02ced19 **was not passing** when committed:
- DVR - Harmonic Oscillator: FAILING
- DVR - Infinite Well: FAILING
- Spectral - Harmonic Oscillator: FAILING
- Spectral - Infinite Well: FAILING (except ground state)

This indicates the tests were added but never properly validated.

## Conclusion

**The spectral method works for different numbers of grid points, BUT ONLY for the very limited case of infinite square well ground states.**

For the user's intended use case of testing the spectral method generally, the answer is:
**The spectral method does NOT work reliably across different potentials and quantum states.**
