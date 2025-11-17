# Tests Directory

This directory contains comprehensive test suites for the QPPW quantum physics simulation.

## Structure

The test directory has been simplified to include only essential files:

- **AccuracyTests.ts** - Comprehensive test suite in `src/common/model/`
- **accuracy-tests.html** - Browser-based test runner with visual interface
- **run-tests.js** - Simple command-line test runner
- **docs/** - Test documentation

## Test Coverage

The comprehensive test suite (`AccuracyTests.ts`) validates all numerical methods (DVR, Spectral, Matrix Numerov, and FGH) against analytical solutions across multiple potential types and grid sizes:

### 1. Harmonic Oscillator
- Tests across grid sizes: 100, 150, 200, 256 points
- All 4 methods tested
- 1% error tolerance

### 2. Finite Square Wells
- **4 configurations tested:**
  - Shallow well: 1nm width, 10eV depth
  - Deep well: 1nm width, 50eV depth
  - Wide shallow well: 2nm width, 10eV depth
  - Narrow medium well: 0.5nm width, 30eV depth
- Grid sizes: 150, 200 points
- DVR, Matrix Numerov, and FGH methods
- 2% error tolerance

### 3. 3D Coulomb Potential (Hydrogen Atom)
- Tests the radial Schrödinger equation for s-waves
- Grid sizes: 150, 200, 256 points
- DVR, Matrix Numerov, and FGH methods
- 3% error tolerance

### 4. Double Square Wells
- **3 configurations tested:**
  - Symmetric, low barrier: 0.5nm wells, 0.3nm barrier, 30eV depth, 10eV barrier
  - Symmetric, medium barrier: 0.5nm wells, 0.5nm barrier, 40eV depth, 20eV barrier
  - Wide wells, low barrier: 0.6nm wells, 0.4nm barrier, 35eV depth, 5eV barrier
- Grid sizes: 200, 256 points
- Methods compared against DVR as reference
- 5% error tolerance

## Running Tests

### Browser-Based Tests (Recommended)

1. Build the project from the root directory:
   ```bash
   npm run build
   ```

2. Open `tests/accuracy-tests.html` in a web browser

3. Click "Run Full Tests" or "Run Quick Check"

### Command-Line Tests

```bash
# From the root directory
node tests/run-tests.js
```

### Development Server

```bash
# Start the dev server from root
npm start

# Then in browser console:
import('./common/model/AccuracyTests.js').then(m => m.runAccuracyTests());
```

## Test Functions

The AccuracyTests module exports two main functions:

- **`runAccuracyTests()`** - Runs the full comprehensive test suite
- **`runQuickAccuracyCheck()`** - Runs a quick subset of tests for rapid validation

## Understanding Test Results

Tests report:
- **Status**: ✓ PASSED or ✗ FAILED
- **Maximum error**: The largest percentage error found
- **Individual energy levels**: Comparison of numerical vs analytical (or reference) values
- **Summary**: Total tests run, passed, and failed

### Error Tolerances

Different potentials have different error tolerances based on their complexity:
- Harmonic Oscillator: 1%
- Finite Square Wells: 2%
- 3D Coulomb: 3%
- Double Square Wells: 5% (compared to DVR reference)

## Adding New Tests

To add new test cases:

1. Open `src/common/model/AccuracyTests.ts`
2. Create a new test function following the pattern of existing tests
3. Add your test to the `runAccuracyTests()` function
4. Rebuild the project: `npm run build`

See `docs/ACCURACY_TESTS.md` for more details on test design.

## Notes

- These tests ensure that numerical methods produce accurate results within acceptable error tolerances
- Tests validate correctness across different grid sizes to ensure robustness
- The comprehensive test suite now covers all major quantum potential types
- All debug and redundant test files have been removed to simplify maintenance
