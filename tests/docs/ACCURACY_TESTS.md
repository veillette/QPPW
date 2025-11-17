# DVR and Spectral Method Accuracy Tests

This document describes the accuracy validation tests for the DVR (Discrete Variable Representation) and Spectral methods used in the QPPW quantum physics simulation.

## Overview

The accuracy tests verify that the numerical methods (DVR and Spectral) produce results within **1% error** of known analytical solutions. This ensures the numerical solvers are working correctly and provide reliable results.

## Test Cases

### 1. Harmonic Oscillator

The quantum harmonic oscillator has an exact analytical solution:

```
E_n = ℏω(n + 1/2)
```

where:
- `n` = quantum number (0, 1, 2, ...)
- `ℏ` = reduced Planck constant
- `ω` = angular frequency

**Test Parameters:**
- Mass: Electron mass (9.109 × 10⁻³¹ kg)
- Angular frequency: 10¹⁵ rad/s
- Grid: -5 nm to +5 nm with 200 points
- States tested: Ground state + 4 excited states

### 2. Infinite Square Well

The infinite square well (particle in a box) has an exact analytical solution:

```
E_n = (n²π²ℏ²) / (2mL²)
```

where:
- `n` = quantum number (1, 2, 3, ...)
- `m` = particle mass
- `L` = well width

**Test Parameters:**
- Mass: Electron mass (9.109 × 10⁻³¹ kg)
- Well width: 1 nm
- Grid: Various configurations for DVR vs Spectral
- States tested: Ground state + 4 excited states

## Running the Tests

### Option 1: Browser-based Testing (Recommended)

1. Build the project:
   ```bash
   npm run build
   ```

2. Open `../accuracy-tests.html` in a web browser

3. Click "Run Full Tests" to run all tests or "Run Quick Check" for a fast validation

The results will be displayed in the browser console with color-coded pass/fail indicators.

### Option 2: Programmatic Testing

Import and run the tests in your code:

```typescript
import { runAccuracyTests, runQuickAccuracyCheck } from './common/model/AccuracyTests.js';

// Run all tests
runAccuracyTests();

// Or run quick check
runQuickAccuracyCheck();
```

### Option 3: Development Server

1. Start the development server:
   ```bash
   npm start
   ```

2. Open the browser console and run:
   ```javascript
   import('./common/model/AccuracyTests.js').then(m => m.runAccuracyTests());
   ```

## Success Criteria

Each test is considered **PASSED** if the numerical result is within **1% error** of the analytical solution.

For example, if the analytical ground state energy is 0.412967 eV:
- ✓ Numerical result of 0.412500 eV is acceptable (0.11% error)
- ✓ Numerical result of 0.408000 eV is acceptable (0.90% error)
- ✗ Numerical result of 0.408000 eV would fail if error > 1%

## Test Output

The test output includes:

1. **Individual test results** for each method and potential:
   - Test name (e.g., "DVR - Harmonic Oscillator")
   - Pass/fail status
   - Maximum error percentage
   - Detailed energy level comparisons

2. **Summary statistics**:
   - Total tests run
   - Number passed
   - Number failed
   - Overall pass/fail status

### Example Output

```
========================================
DVR and Spectral Method Accuracy Tests
========================================
Tolerance: 1% error from analytical solutions

=== DVR - Harmonic Oscillator ===
Status: ✓ PASSED
Maximum error: 0.0234%

Testing 5 energy levels:
  E_0: ✓ PASS - Numerical: 0.412967 eV, Analytical: 0.412967 eV, Error: 0.0001%
  E_1: ✓ PASS - Numerical: 1.238902 eV, Analytical: 1.238902 eV, Error: 0.0002%
  E_2: ✓ PASS - Numerical: 2.064836 eV, Analytical: 2.064837 eV, Error: 0.0003%
  E_3: ✓ PASS - Numerical: 2.890771 eV, Analytical: 2.890771 eV, Error: 0.0004%
  E_4: ✓ PASS - Numerical: 3.716705 eV, Analytical: 3.716706 eV, Error: 0.0005%

...

========================================
Test Summary
========================================
Total tests: 4
Passed: 4
Failed: 0

✓ All tests passed!
Both DVR and Spectral methods produce results within 1% of analytical solutions.
========================================
```

## Implementation Details

### Test File Location

- **Source:** `src/common/model/AccuracyTests.ts`
- **Compiled:** `dist/common/model/AccuracyTests.js`
- **Test Runner:** `../accuracy-tests.html`

### Methods Tested

1. **DVR (Discrete Variable Representation)**
   - File: `src/common/model/DVRSolver.ts`
   - Uses Colbert-Miller formula for kinetic energy
   - Matrix diagonalization approach

2. **Spectral (Chebyshev)**
   - File: `src/common/model/SpectralSolver.ts`
   - Uses Chebyshev polynomial expansion
   - Provides spectral (exponential) convergence

### Analytical Solutions

The analytical solutions are implemented in:
- `src/common/model/analytical-solutions/harmonic-oscillator.ts`
- `src/common/model/analytical-solutions/infinite-square-well.ts`

## Troubleshooting

### Tests Fail to Run

1. Ensure the project is built: `npm run build`
2. Check that all dependencies are installed: `npm install`
3. Verify TypeScript compilation succeeds: `npm run check`

### Tests Show High Error

If tests show errors exceeding 1%:

1. **Check grid resolution:** Increase `numPoints` in grid configuration
2. **Check grid extent:** Ensure grid covers the region where wavefunction is significant
3. **Check potential definition:** Verify potential function is correctly implemented
4. **Check numerical precision:** Ensure calculations use sufficient precision

### Browser Console Errors

If the browser shows module loading errors:
1. Rebuild the project: `npm run build`
2. Clear browser cache
3. Check browser console for specific error messages

## Adding New Tests

To add new test cases:

1. Add a new test function in `AccuracyTests.ts`:
   ```typescript
   function testDVRNewPotential(): TestResult {
     // Implementation
   }
   ```

2. Call it from `runAccuracyTests()`:
   ```typescript
   results.push(testDVRNewPotential());
   ```

3. Ensure an analytical solution exists in `analytical-solutions/`

4. Rebuild and run tests to verify

## References

- **DVR Method:** Colbert & Miller, J. Chem. Phys. 96, 1982 (1992)
- **Spectral Method:** Boyd, "Chebyshev and Fourier Spectral Methods", 2nd ed. (2001)
- **Quantum Mechanics:** Griffiths, "Introduction to Quantum Mechanics", 3rd ed.

## License

MIT License - Same as the main QPPW project
