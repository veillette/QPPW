/**
 * Simple test runner for DVR accuracy tests
 * Run with: node test-dvr.js
 */

import { runAccuracyTests } from './dist/common/model/AccuracyTests.js';

console.log('Building project first...\n');

// Run the accuracy tests
runAccuracyTests();
