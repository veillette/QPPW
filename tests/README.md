# Tests Directory

This directory contains test files and debugging utilities for the QPPW quantum physics simulation.

## Contents

### Test Files

- **test-*.js/ts** - Unit and integration tests for various solvers and components
- **debug-*.js/mjs** - Debugging scripts for investigating solver behavior
- **check-*.js** - Validation scripts
- **diagnose-*.js** - Diagnostic tools for troubleshooting
- **verify-*.ts** - Verification scripts for specific features

### Test Runners

- **run-tests.js** - Main test runner
- **run-accuracy-test.mjs** - Accuracy test runner for numerical methods
- **accuracy-tests.html** - Browser-based test runner with visual interface

## Documentation

See the `docs/` subdirectory for detailed documentation:
- **ACCURACY_TESTS.md** - Instructions for running and understanding accuracy tests
- **ACCURACY_VERIFICATION.md** - Verification procedures for numerical solvers

## Running Tests

### Browser-Based Tests (Recommended)

1. Build the project from the root directory:
   ```bash
   npm run build
   ```

2. Open `accuracy-tests.html` in a web browser

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

## Adding New Tests

1. Create a new test file in this directory (e.g., `test-new-feature.js`)
2. Follow the existing test patterns
3. Update this README if needed
4. See `docs/ACCURACY_TESTS.md` for details on adding accuracy tests

## Note

These tests are primarily for development and validation purposes. They ensure that numerical methods (DVR, Spectral, Numerov, FGH) produce accurate results within acceptable error tolerances.
