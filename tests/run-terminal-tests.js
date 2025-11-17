#!/usr/bin/env node
/**
 * Terminal-based test runner for QPPW accuracy tests
 *
 * This script runs the comprehensive test suite directly in the terminal
 * without requiring a browser.
 *
 * Usage:
 *   npx tsx tests/run-terminal-tests.js
 *   npm test
 */

console.log('╔════════════════════════════════════════════════════════════╗');
console.log('║     QPPW Quantum Mechanics Accuracy Test Suite            ║');
console.log('║     Running in Terminal Mode                               ║');
console.log('╚════════════════════════════════════════════════════════════╝\n');

// Note: This file must be run with tsx to handle TypeScript imports
// Run with: npx tsx tests/run-terminal-tests.js

import { runAccuracyTests } from '../src/common/model/AccuracyTests.js';

// Run the full test suite
runAccuracyTests();
