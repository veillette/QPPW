# Performance Monitoring

This document describes the performance monitoring and memoization features added to the QPPW (Quantum Particle in Potential Wells) simulation.

## Overview

The performance monitoring system helps identify and track slow operations that may impact the user experience. It consists of three main utilities:

1. **Logger** - Simple logging utility with different log levels
2. **PerformanceMonitor** - Tracks execution times and warns about slow operations
3. **MemoizationCache** - Caches expensive calculations to avoid redundant computation

## Components

### Logger (`src/common/utils/Logger.ts`)

A simple logging utility that provides different log levels.

```typescript
import Logger, { LogLevel } from "./common/utils/Logger.js";

// Set log level (default is INFO)
Logger.setLevel(LogLevel.DEBUG);

// Log messages at different levels
Logger.debug("Debug message", { data: "value" });
Logger.info("Info message");
Logger.warn("Warning message");
Logger.error("Error message");
```

**Log Levels:**

- `DEBUG` - Detailed debugging information
- `INFO` - General informational messages
- `WARN` - Warning messages
- `ERROR` - Error messages
- `NONE` - Disable all logging

### PerformanceMonitor (`src/common/utils/PerformanceMonitor.ts`)

Tracks execution times of operations and warns when they exceed the frame budget (16.67ms at 60fps).

```typescript
import PerformanceMonitor from "./common/utils/PerformanceMonitor.js";

// Measure a synchronous operation
const result = PerformanceMonitor.measure("solve-schrodinger", () => {
  return solver.solveNumerical(potential, mass, numStates, grid);
});

// Measure an async operation
const data = await PerformanceMonitor.measureAsync("load-data", async () => {
  return await fetch("/api/data").then((r) => r.json());
});

// Get statistics for a specific operation
const stats = PerformanceMonitor.getStats("solve-schrodinger");
console.log(stats); // { count, mean, max, min, total }

// Get all statistics
const allStats = PerformanceMonitor.getAllStats();

// Log a summary of all metrics
PerformanceMonitor.logSummary();

// Clear statistics
PerformanceMonitor.clearStats("solve-schrodinger"); // Clear specific operation
PerformanceMonitor.clearAllStats(); // Clear all statistics

// Enable/disable monitoring
PerformanceMonitor.setEnabled(false); // Disable
PerformanceMonitor.setEnabled(true); // Enable
```

**Features:**

- Automatically warns when operations exceed 16.67ms (1 frame at 60fps)
- Collects statistics: count, mean, min, max, total time
- Can be enabled/disabled globally
- Supports both sync and async functions

### MemoizationCache (`src/common/utils/MemoizationHelper.ts`)

Caches the results of expensive calculations to avoid redundant computation.

```typescript
import MemoizationCache, { memoize } from "./common/utils/MemoizationHelper.js";

// Using MemoizationCache directly
class AnalyticalSolution {
  private wavefunctionCache = new MemoizationCache<number[]>(100); // Max 100 entries

  getWavefunction(n: number, width: number, depth: number): number[] {
    const key = MemoizationCache.generateKey(n, width, depth);

    return this.wavefunctionCache.getOrCompute(key, () => {
      // Expensive calculation only runs if not cached
      return this.calculateWavefunction(n, width, depth);
    });
  }

  private calculateWavefunction(
    n: number,
    width: number,
    depth: number,
  ): number[] {
    // ... expensive calculation ...
  }
}

// Using the memoize function
const memoizedCalculation = memoize(
  (n: number, width: number, depth: number) =>
    expensiveCalculation(n, width, depth),
  (n, width, depth) => `${n}-${width}-${depth}`, // Custom key generator
  100, // Max cache size
);

const result = memoizedCalculation(1, 2.0, 3.0);
```

**Features:**

- Automatic cache key generation
- Configurable maximum cache size (FIFO eviction)
- Generic type support
- Helper function for functional-style memoization

## Usage in the Codebase

### Quantum Solver Performance Monitoring

The `QuantumBoundStateSolver` has been enhanced with performance monitoring:

```typescript
import { solveQuantumBound } from "./common/model/QuantumBoundStateSolver.js";
import PerformanceMonitor from "./common/utils/PerformanceMonitor.js";

// Solve the quantum bound state problem
// Performance will be automatically monitored
const result = solveQuantumBound(potential, mass, numStates, gridConfig);

// Check performance statistics
PerformanceMonitor.logSummary();
```

**Monitored Operations:**

- `solve-quantum-bound` - Top-level solve function
- `find-eigenstate-{n}` - Finding individual eigenstates
- `find-multiple-eigenstates-{nStates}` - Finding multiple eigenstates
- `calculate-wavefunction` - Wavefunction calculation

### Example: Memoization for Analytical Solutions

For analytical solutions where wavefunctions are calculated repeatedly with the same parameters:

```typescript
class FiniteWellSolution extends AnalyticalSolution {
  private memoizedWavefunctions = new MemoizationCache<number[]>();

  public getWavefunction(n: number): number[] | null {
    const key = MemoizationCache.generateKey(n, this.wellWidth, this.wellDepth);

    if (this.memoizedWavefunctions.has(key)) {
      return this.memoizedWavefunctions.get(key)!;
    }

    const wavefunction = this.calculateWavefunction(n);
    this.memoizedWavefunctions.set(key, wavefunction);
    return wavefunction;
  }
}
```

## Performance Threshold

The default performance threshold is **16.67ms**, which corresponds to one frame at 60fps. Operations exceeding this threshold will trigger a warning:

```
[WARN] Slow operation: solve-schrodinger took 25.42ms (>16.67ms threshold)
```

This helps identify performance bottlenecks that may cause frame drops or UI lag.

## Best Practices

1. **Use performance monitoring for expensive operations** - Solver computations, large data processing, complex calculations

2. **Use memoization for pure functions** - Functions that return the same result for the same inputs

3. **Set appropriate cache sizes** - Balance memory usage vs. cache hit rate

4. **Review performance logs regularly** - Use `PerformanceMonitor.logSummary()` during development

5. **Disable in production if needed** - Use `PerformanceMonitor.setEnabled(false)` to reduce overhead

## Example Output

```
[INFO] === Performance Metrics Summary ===
[INFO] solve-quantum-bound: count=1, mean=45.23ms, min=45.23ms, max=45.23ms, total=45.23ms
[INFO] find-eigenstate-1: count=1, mean=5.12ms, min=5.12ms, max=5.12ms, total=5.12ms
[INFO] find-eigenstate-2: count=1, mean=6.34ms, min=6.34ms, max=6.34ms, total=6.34ms
[INFO] calculate-wavefunction: count=10, mean=3.45ms, min=2.12ms, max=5.67ms, total=34.50ms
[INFO] ===================================
```

## See Also

- [Solver Documentation](./SOLVER_DOCUMENTATION.md)
- [Code Review](../CODE_REVIEW.md)
- [Improvement Suggestions](../IMPROVEMENT_SUGGESTIONS.md)
