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

### Method 1: Browser-Based Tests (Recommended)

This is the easiest way to run tests with a visual interface.

1. **Build the project** from the root directory:
   ```bash
   npm run build
   ```

2. **Open the test runner** in a web browser:
   ```bash
   # Open tests/accuracy-tests.html in your browser
   # On macOS:
   open tests/accuracy-tests.html
   # On Linux:
   xdg-open tests/accuracy-tests.html
   # On Windows:
   start tests/accuracy-tests.html
   ```

3. **Run tests** by clicking one of the buttons:
   - **"Run Full Tests"** - Runs comprehensive test suite (~60+ tests)
   - **"Run Quick Check"** - Runs a subset for rapid validation

### Method 2: Command-Line Tests

Run tests directly from the command line.

```bash
# From the root directory
node tests/run-tests.js
```

This will execute the comprehensive test suite and display results in the terminal, including:
- Pass/fail status for each test
- Numerical vs analytical energy comparisons
- Error percentages
- **Execution times for each method**
- **Performance summary with timing statistics**

### Method 3: Development Server + Browser Console

Useful during active development.

1. **Start the development server** from the root directory:
   ```bash
   npm start
   ```

2. **Open your browser** and navigate to the local server URL (usually `http://localhost:8080`)

3. **Open browser console** (F12 or Cmd+Option+I) and run:
   ```javascript
   // Run full comprehensive test suite
   import('./common/model/AccuracyTests.js').then(m => m.runAccuracyTests());

   // Or run quick check only
   import('./common/model/AccuracyTests.js').then(m => m.runQuickAccuracyCheck());
   ```

### Method 4: Direct TypeScript Execution (For Development)

If you have `tsx` installed, you can run tests directly on the TypeScript source:

```bash
# Install tsx if needed
npm install -g tsx

# Run tests directly
npx tsx src/common/model/AccuracyTests.ts
```

## Test Functions

The AccuracyTests module exports two main functions:

- **`runAccuracyTests()`** - Runs the full comprehensive test suite
- **`runQuickAccuracyCheck()`** - Runs a quick subset of tests for rapid validation

## Understanding Test Results

### Individual Test Output

Each test reports:
- **Status**: ✓ PASSED or ✗ FAILED
- **Maximum error**: The largest percentage error found across all energy levels
- **Execution time**: Time taken to solve in milliseconds
- **Individual energy levels**: Detailed comparison showing:
  - Numerical result (in eV)
  - Analytical/reference result (in eV)
  - Percentage error

Example output:
```
=== DVR - Harmonic Oscillator (N=150) ===
Status: ✓ PASSED
Maximum error: 0.0234%

Testing 5 energy levels with 150 grid points:
Execution time: 12.45 ms
  E_0: ✓ Num: 0.413567 eV, Ana: 0.413470 eV, Err: 0.0234%
  E_1: ✓ Num: 1.240567 eV, Ana: 1.240410 eV, Err: 0.0127%
  ...
```

### Performance Summary

At the end of the test run, a performance summary shows timing statistics for each method:
```
--- Performance Summary ---
DVR:
  Average: 15.34 ms | Min: 8.21 ms | Max: 45.67 ms | Total: 245.44 ms
MatrixNumerov:
  Average: 12.89 ms | Min: 7.45 ms | Max: 38.23 ms | Total: 206.24 ms
Spectral:
  Average: 18.76 ms | Min: 10.34 ms | Max: 52.11 ms | Total: 300.16 ms
FGH:
  Average: 22.14 ms | Min: 11.56 ms | Max: 68.92 ms | Total: 354.24 ms
```

This allows you to:
- Compare performance across different methods
- Identify which method is fastest for your use case
- Track performance regressions over time

### Error Tolerances

Different potentials have different error tolerances based on their complexity:
- **Harmonic Oscillator**: 1% - Very well-behaved, high accuracy expected
- **Finite Square Wells**: 2% - Moderate complexity with discontinuous potential
- **3D Coulomb**: 3% - Singular potential requires careful handling
- **Double Square Wells**: 5% - Compared against DVR reference (no analytical solution)

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
