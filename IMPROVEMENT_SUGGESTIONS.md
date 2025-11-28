# Actionable Improvement Suggestions

This document provides specific, implementable code improvements for the QPPW codebase.

---

## Quick Wins (Can be implemented immediately)

### 1. Add Automated Tests to CI Pipeline (15 minutes)

**File**: `.github/workflows/deploy.yml`

**Current** (lines 28-33):
```yaml
- name: Install dependencies
  run: npm ci
- name: Type check
  run: npm run check
- name: Lint
  run: npm run lint
```

**Improved**:
```yaml
- name: Install dependencies
  run: npm ci
- name: Run tests
  run: npm test
- name: Type check
  run: npm run check
- name: Lint
  run: npm run lint
- name: Security audit
  run: npm audit --audit-level=high
- name: Check formatting
  run: npm run format:check
```

---

### 2. Create Logger Utility (30 minutes)

**New File**: `src/common/utils/Logger.ts`

```typescript
/**
 * Centralized logging utility with level-based filtering.
 * Prevents console spam in production while enabling debug logs in development.
 */

export enum LogLevel {
  DEBUG = 0,
  INFO = 1,
  WARN = 2,
  ERROR = 3,
  NONE = 4,
}

export class Logger {
  private static level: LogLevel =
    process.env.NODE_ENV === "production" ? LogLevel.WARN : LogLevel.DEBUG;

  /**
   * Set the minimum log level. Messages below this level will be suppressed.
   */
  static setLevel(level: LogLevel): void {
    this.level = level;
  }

  /**
   * Log debug information (development only).
   * Use for verbose debugging information.
   */
  static debug(message: string, context?: Record<string, any>): void {
    if (this.level <= LogLevel.DEBUG) {
      console.log(`[DEBUG] ${message}`, context ?? "");
    }
  }

  /**
   * Log informational messages.
   * Use for significant application events.
   */
  static info(message: string, context?: Record<string, any>): void {
    if (this.level <= LogLevel.INFO) {
      console.info(`[INFO] ${message}`, context ?? "");
    }
  }

  /**
   * Log warning messages.
   * Use for recoverable errors or unexpected conditions.
   */
  static warn(message: string, context?: Record<string, any>): void {
    if (this.level <= LogLevel.WARN) {
      console.warn(`[WARN] ${message}`, context ?? "");
    }
  }

  /**
   * Log error messages.
   * Use for unrecoverable errors.
   */
  static error(
    message: string,
    error?: Error,
    context?: Record<string, any>,
  ): void {
    if (this.level <= LogLevel.ERROR) {
      console.error(`[ERROR] ${message}`, error ?? "", context ?? "");
    }
  }

  /**
   * Create a performance timer for measuring operation duration.
   * Automatically logs if operation exceeds threshold.
   *
   * @param operationName - Name of the operation being timed
   * @param warnThresholdMs - Log warning if operation exceeds this (default: 16.67ms = 1 frame)
   */
  static time(
    operationName: string,
    warnThresholdMs = 16.67,
  ): { end: () => void } {
    const start = performance.now();
    return {
      end: () => {
        const duration = performance.now() - start;
        if (duration > warnThresholdMs) {
          this.warn(`Slow operation: ${operationName}`, {
            duration: `${duration.toFixed(2)}ms`,
          });
        } else {
          this.debug(`${operationName} completed`, {
            duration: `${duration.toFixed(2)}ms`,
          });
        }
      },
    };
  }
}

export default Logger;
```

**Usage Example** - Update `src/common/model/Schrodinger1DSolver.ts`:

```typescript
import Logger from "../utils/Logger.js";

public solveNumerical(
  potential: PotentialFunction,
  mass: number,
  numStates: number,
  gridConfig: GridConfig,
  energyRange?: [number, number],
): BoundStateResult {
  const timer = Logger.time("Numerical solve");

  Logger.debug("Starting numerical solve", {
    mass,
    numStates,
    method: this.numericalMethod,
    gridPoints: gridConfig.numPoints,
  });

  try {
    // ... existing solve logic
    const result = /* ... */;
    timer.end();
    return result;
  } catch (error) {
    Logger.error("Numerical solve failed", error as Error, {
      mass,
      numStates,
      gridConfig,
    });
    throw error;
  }
}
```

---

### 3. Add Input Validation Utility (1 hour)

**New File**: `src/common/utils/Validation.ts`

```typescript
/**
 * Validation utilities for input parameters.
 * Throws descriptive errors for invalid inputs.
 */

export class ValidationError extends Error {
  constructor(
    message: string,
    public readonly field: string,
    public readonly value: any,
  ) {
    super(message);
    this.name = "ValidationError";
  }
}

/**
 * Validates that a value is a finite number.
 */
export function validateFiniteNumber(
  value: number,
  fieldName: string,
): void {
  if (!Number.isFinite(value)) {
    throw new ValidationError(
      `${fieldName} must be a finite number, got ${value}`,
      fieldName,
      value,
    );
  }
}

/**
 * Validates that a value is positive (> 0).
 */
export function validatePositive(value: number, fieldName: string): void {
  validateFiniteNumber(value, fieldName);
  if (value <= 0) {
    throw new ValidationError(
      `${fieldName} must be positive, got ${value}`,
      fieldName,
      value,
    );
  }
}

/**
 * Validates that a value is non-negative (>= 0).
 */
export function validateNonNegative(
  value: number,
  fieldName: string,
): void {
  validateFiniteNumber(value, fieldName);
  if (value < 0) {
    throw new ValidationError(
      `${fieldName} must be non-negative, got ${value}`,
      fieldName,
      value,
    );
  }
}

/**
 * Validates that a value is within a specified range [min, max].
 */
export function validateRange(
  value: number,
  min: number,
  max: number,
  fieldName: string,
): void {
  validateFiniteNumber(value, fieldName);
  if (value < min || value > max) {
    throw new ValidationError(
      `${fieldName} must be between ${min} and ${max}, got ${value}`,
      fieldName,
      value,
    );
  }
}

/**
 * Validates that a value is an integer.
 */
export function validateInteger(value: number, fieldName: string): void {
  validateFiniteNumber(value, fieldName);
  if (!Number.isInteger(value)) {
    throw new ValidationError(
      `${fieldName} must be an integer, got ${value}`,
      fieldName,
      value,
    );
  }
}

/**
 * Validates that an array has the expected length.
 */
export function validateArrayLength<T>(
  array: T[],
  expectedLength: number,
  fieldName: string,
): void {
  if (array.length !== expectedLength) {
    throw new ValidationError(
      `${fieldName} must have length ${expectedLength}, got ${array.length}`,
      fieldName,
      array,
    );
  }
}

/**
 * Validates grid configuration parameters.
 */
export function validateGridConfig(gridConfig: {
  xMin: number;
  xMax: number;
  numPoints: number;
}): void {
  validateFiniteNumber(gridConfig.xMin, "gridConfig.xMin");
  validateFiniteNumber(gridConfig.xMax, "gridConfig.xMax");
  validateInteger(gridConfig.numPoints, "gridConfig.numPoints");
  validateRange(gridConfig.numPoints, 10, 10000, "gridConfig.numPoints");

  if (gridConfig.xMax <= gridConfig.xMin) {
    throw new ValidationError(
      `gridConfig.xMax (${gridConfig.xMax}) must be greater than gridConfig.xMin (${gridConfig.xMin})`,
      "gridConfig",
      gridConfig,
    );
  }
}
```

**Usage** - Update `src/common/model/Schrodinger1DSolver.ts`:

```typescript
import {
  validatePositive,
  validateGridConfig,
  ValidationError,
} from "../utils/Validation.js";

public solveNumerical(
  potential: PotentialFunction,
  mass: number,
  numStates: number,
  gridConfig: GridConfig,
  energyRange?: [number, number],
): BoundStateResult {
  // Validate inputs
  validatePositive(mass, "mass");
  validatePositive(numStates, "numStates");
  validateGridConfig(gridConfig);

  if (energyRange) {
    validateFiniteNumber(energyRange[0], "energyRange[0]");
    validateFiniteNumber(energyRange[1], "energyRange[1]");
    if (energyRange[1] <= energyRange[0]) {
      throw new ValidationError(
        `energyRange[1] must be greater than energyRange[0]`,
        "energyRange",
        energyRange,
      );
    }
  }

  // ... rest of method
}
```

---

### 4. Add Dependabot for Dependency Updates (5 minutes)

**New File**: `.github/dependabot.yml`

```yaml
version: 2
updates:
  # Maintain npm dependencies
  - package-ecosystem: "npm"
    directory: "/"
    schedule:
      interval: "weekly"
      day: "monday"
      time: "09:00"
    open-pull-requests-limit: 5
    reviewers:
      - "veillette"
    labels:
      - "dependencies"
      - "automated"
    # Group minor and patch updates together
    groups:
      development-dependencies:
        patterns:
          - "@typescript/*"
          - "eslint*"
          - "prettier"
          - "typescript"
          - "vite*"
        update-types:
          - "minor"
          - "patch"

  # Maintain GitHub Actions
  - package-ecosystem: "github-actions"
    directory: "/"
    schedule:
      interval: "monthly"
    reviewers:
      - "veillette"
    labels:
      - "dependencies"
      - "github-actions"
```

---

## Medium-Term Improvements (1-2 days each)

### 5. Refactor Code Duplication in Schrodinger1DSolver

**Current Issue**: Lines 194-353 and 573-710 are nearly identical.

**Solution**: Create a configuration-based factory.

**New File**: `src/common/model/PotentialFactory.ts`

```typescript
import {
  AnalyticalSolution,
  InfiniteSquareWellSolution,
  FiniteSquareWellSolution,
  HarmonicOscillatorSolution,
  MorsePotentialSolution,
  PoschlTellerPotentialSolution,
  RosenMorsePotentialSolution,
  EckartPotentialSolution,
  AsymmetricTrianglePotentialSolution,
  Coulomb1DPotentialSolution,
  Coulomb3DPotentialSolution,
  TriangularPotentialSolution,
} from "./analytical-solutions/index.js";

import {
  BasePotential,
  InfiniteSquareWellPotential,
  FiniteSquareWellPotential,
  HarmonicOscillatorPotential,
  MorsePotential,
  PoschlTellerPotential,
  RosenMorsePotential,
  EckartPotential,
  AsymmetricTrianglePotential,
  Coulomb1DPotential,
  Coulomb3DPotential,
  TriangularPotential,
} from "./potentials/index.js";

import { PotentialType, WellParameters } from "./PotentialFunction.js";

/**
 * Configuration for creating analytical solutions and potential classes.
 */
interface PotentialConfig {
  /** Constructor for analytical solution class */
  solutionClass: new (...args: any[]) => AnalyticalSolution;
  /** Constructor for potential class */
  potentialClass: new (...args: any[]) => BasePotential;
  /** Extract constructor arguments from well parameters */
  extractArgs: (params: WellParameters, mass: number) => any[];
  /** Required parameter names */
  requiredParams: (keyof WellParameters)[];
}

/**
 * Registry of all potential type configurations.
 * Single source of truth for creating analytical solutions and potentials.
 */
const POTENTIAL_CONFIGS: Partial<Record<PotentialType, PotentialConfig>> = {
  [PotentialType.INFINITE_WELL]: {
    solutionClass: InfiniteSquareWellSolution,
    potentialClass: InfiniteSquareWellPotential,
    extractArgs: (p, m) => [p.wellWidth!, m],
    requiredParams: ["wellWidth"],
  },

  [PotentialType.FINITE_WELL]: {
    solutionClass: FiniteSquareWellSolution,
    potentialClass: FiniteSquareWellPotential,
    extractArgs: (p, m) => [p.wellWidth!, p.wellDepth!, m],
    requiredParams: ["wellWidth", "wellDepth"],
  },

  [PotentialType.HARMONIC_OSCILLATOR]: {
    solutionClass: HarmonicOscillatorSolution,
    potentialClass: HarmonicOscillatorPotential,
    extractArgs: (p, m) => [p.springConstant!, m],
    requiredParams: ["springConstant"],
  },

  [PotentialType.MORSE]: {
    solutionClass: MorsePotentialSolution,
    potentialClass: MorsePotential,
    extractArgs: (p, m) => [
      p.dissociationEnergy!,
      p.wellWidth!,
      p.equilibriumPosition!,
      m,
    ],
    requiredParams: ["dissociationEnergy", "wellWidth", "equilibriumPosition"],
  },

  [PotentialType.POSCHL_TELLER]: {
    solutionClass: PoschlTellerPotentialSolution,
    potentialClass: PoschlTellerPotential,
    extractArgs: (p, m) => [p.potentialDepth!, p.wellWidth!, m],
    requiredParams: ["potentialDepth", "wellWidth"],
  },

  [PotentialType.ROSEN_MORSE]: {
    solutionClass: RosenMorsePotentialSolution,
    potentialClass: RosenMorsePotential,
    extractArgs: (p, m) => [
      p.potentialDepth!,
      p.barrierHeight!,
      p.wellWidth!,
      m,
    ],
    requiredParams: ["potentialDepth", "barrierHeight", "wellWidth"],
  },

  [PotentialType.ECKART]: {
    solutionClass: EckartPotentialSolution,
    potentialClass: EckartPotential,
    extractArgs: (p, m) => [
      p.potentialDepth!,
      p.barrierHeight!,
      p.wellWidth!,
      m,
    ],
    requiredParams: ["potentialDepth", "barrierHeight", "wellWidth"],
  },

  [PotentialType.ASYMMETRIC_TRIANGLE]: {
    solutionClass: AsymmetricTrianglePotentialSolution,
    potentialClass: AsymmetricTrianglePotential,
    extractArgs: (p, m) => [p.slope!, p.wellWidth!, m],
    requiredParams: ["slope", "wellWidth"],
  },

  [PotentialType.COULOMB_1D]: {
    solutionClass: Coulomb1DPotentialSolution,
    potentialClass: Coulomb1DPotential,
    extractArgs: (p, m) => [p.coulombStrength!, m],
    requiredParams: ["coulombStrength"],
  },

  [PotentialType.COULOMB_3D]: {
    solutionClass: Coulomb3DPotentialSolution,
    potentialClass: Coulomb3DPotential,
    extractArgs: (p, m) => [p.coulombStrength!, m],
    requiredParams: ["coulombStrength"],
  },

  [PotentialType.TRIANGULAR]: {
    solutionClass: TriangularPotentialSolution,
    potentialClass: TriangularPotential,
    extractArgs: (p, m) => [p.wellDepth!, p.wellWidth!, p.energyOffset!, m],
    requiredParams: ["wellDepth", "wellWidth", "energyOffset"],
  },
};

/**
 * Factory for creating analytical solutions and potential instances.
 * Eliminates code duplication in Schrodinger1DSolver.
 */
export class PotentialFactory {
  /**
   * Check if all required parameters are present.
   */
  private static hasRequiredParams(
    params: WellParameters,
    required: (keyof WellParameters)[],
  ): boolean {
    return required.every((key) => params[key] !== undefined);
  }

  /**
   * Create an analytical solution instance.
   * Returns null if the potential type doesn't support analytical solutions
   * or required parameters are missing.
   */
  static createAnalyticalSolution(
    wellParams: WellParameters,
    mass: number,
  ): AnalyticalSolution | null {
    const config = POTENTIAL_CONFIGS[wellParams.type];

    if (!config) {
      return null; // No analytical solution for this type
    }

    if (!this.hasRequiredParams(wellParams, config.requiredParams)) {
      return null; // Missing required parameters
    }

    const args = config.extractArgs(wellParams, mass);
    return new config.solutionClass(...args);
  }

  /**
   * Create a potential class instance.
   * Returns null if the potential type isn't supported
   * or required parameters are missing.
   */
  static createPotential(
    wellParams: WellParameters,
    mass: number,
  ): BasePotential | null {
    const config = POTENTIAL_CONFIGS[wellParams.type];

    if (!config) {
      return null; // No potential class for this type
    }

    if (!this.hasRequiredParams(wellParams, config.requiredParams)) {
      return null; // Missing required parameters
    }

    const args = config.extractArgs(wellParams, mass);
    return new config.potentialClass(...args);
  }
}
```

**Then simplify `Schrodinger1DSolver.ts`**:

```typescript
import { PotentialFactory } from "./PotentialFactory.js";

// BEFORE: 160 lines of switch statements
// AFTER: 3 lines
private createAnalyticalSolution(
  wellParams: WellParameters,
  mass: number,
): AnalyticalSolution | null {
  return PotentialFactory.createAnalyticalSolution(wellParams, mass);
}

public createPotential(
  wellParams: WellParameters,
  mass: number,
): BasePotential | null {
  return PotentialFactory.createPotential(wellParams, mass);
}
```

**Benefits**:
- Eliminates 300+ lines of duplication
- Single source of truth for potential configurations
- Easier to add new potential types
- More testable

---

### 6. Add Basic Accessibility (2-3 days)

**Step 1**: Add keyboard navigation to interactive elements.

**File**: `src/common/view/WaveFunctionChartNode.ts` (around line 88)

**Current**:
```typescript
this.leftMarkerHandle = new Circle(8, {
  fill: QPPWColors.markerFillProperty,
  stroke: QPPWColors.markerStrokeProperty,
  cursor: "ew-resize",
});
```

**Improved**:
```typescript
this.leftMarkerHandle = new Circle(8, {
  fill: QPPWColors.markerFillProperty,
  stroke: QPPWColors.markerStrokeProperty,
  cursor: "ew-resize",
  // Accessibility improvements
  tagName: "button",
  ariaLabel: "Left boundary marker for area measurement",
  innerContent: "Drag to adjust left boundary, or use arrow keys",
  focusable: true,
});

// Add keyboard support
const keyStep = 0.05; // nm per keypress
this.leftMarkerHandle.addInputListener({
  keydown: (event: KeyboardEvent) => {
    if (event.key === "ArrowLeft") {
      this.leftMarkerXProperty.value = Math.max(
        this.xMinProperty.value,
        this.leftMarkerXProperty.value - keyStep,
      );
      event.preventDefault();
    } else if (event.key === "ArrowRight") {
      this.leftMarkerXProperty.value = Math.min(
        this.rightMarkerXProperty.value - keyStep,
        this.leftMarkerXProperty.value + keyStep,
      );
      event.preventDefault();
    }
  },
});
```

**Step 2**: Add live region announcer.

**New File**: `src/common/utils/LiveRegionAnnouncer.ts`

```typescript
/**
 * Announcer for screen reader users using ARIA live regions.
 * Provides non-visual feedback for state changes.
 */
export class LiveRegionAnnouncer {
  private liveRegion: HTMLElement;
  private timeoutId: number | null = null;

  constructor() {
    // Create off-screen live region
    this.liveRegion = document.createElement("div");
    this.liveRegion.setAttribute("aria-live", "polite");
    this.liveRegion.setAttribute("aria-atomic", "true");
    this.liveRegion.className = "screen-reader-only";

    // Position off-screen but accessible to screen readers
    Object.assign(this.liveRegion.style, {
      position: "absolute",
      left: "-10000px",
      width: "1px",
      height: "1px",
      overflow: "hidden",
    });

    document.body.appendChild(this.liveRegion);
  }

  /**
   * Announce a message to screen readers.
   * Messages are debounced to prevent flooding.
   */
  announce(message: string, debounceMs = 100): void {
    if (this.timeoutId !== null) {
      clearTimeout(this.timeoutId);
    }

    this.timeoutId = window.setTimeout(() => {
      this.liveRegion.textContent = message;
      this.timeoutId = null;
    }, debounceMs);
  }

  /**
   * Announce urgent message (assertive).
   * Use sparingly, only for critical information.
   */
  announceUrgent(message: string): void {
    this.liveRegion.setAttribute("aria-live", "assertive");
    this.liveRegion.textContent = message;

    // Reset to polite after announcement
    setTimeout(() => {
      this.liveRegion.setAttribute("aria-live", "polite");
    }, 100);
  }

  /**
   * Clear the live region.
   */
  clear(): void {
    if (this.timeoutId !== null) {
      clearTimeout(this.timeoutId);
      this.timeoutId = null;
    }
    this.liveRegion.textContent = "";
  }
}

// Singleton instance
export const announcer = new LiveRegionAnnouncer();
```

**Usage** - Update energy level selection:

```typescript
import { announcer } from "../utils/LiveRegionAnnouncer.js";

// When energy level changes:
this.model.selectedEnergyLevelIndexProperty.link((index) => {
  const energy = this.model.getEnergyLevel(index + 1);
  announcer.announce(
    `Energy level ${index + 1} selected. Energy: ${energy.toFixed(3)} electron volts`,
  );
});
```

---

## Testing Improvements

### 7. Add Unit Tests for View Components

**New File**: `tests/view/test-wavefunction-chart-node.ts`

```typescript
import { describe, it, expect, beforeEach } from "node:test";
import { WaveFunctionChartNode } from "../../src/common/view/WaveFunctionChartNode.js";
import { OneWellModel } from "../../src/one-well/model/OneWellModel.js";
import { OneWellViewState } from "../../src/one-well/view/OneWellViewState.js";

describe("WaveFunctionChartNode", () => {
  let model: OneWellModel;
  let viewState: OneWellViewState;
  let chartNode: WaveFunctionChartNode;

  beforeEach(() => {
    model = new OneWellModel();
    viewState = new OneWellViewState();
    chartNode = new WaveFunctionChartNode(model, viewState);
  });

  it("should create chart with default dimensions", () => {
    expect(chartNode).toBeDefined();
    expect(chartNode.bounds.width).toBeGreaterThan(0);
    expect(chartNode.bounds.height).toBeGreaterThan(0);
  });

  it("should update when model energy level changes", () => {
    const initialLevel = model.selectedEnergyLevelIndexProperty.value;
    model.selectedEnergyLevelIndexProperty.value = initialLevel + 1;

    // Chart should trigger update
    // (Implementation depends on how updates are tracked)
  });

  it("should show area tool when property is true", () => {
    chartNode.showAreaToolProperty.value = false;
    expect(chartNode.showAreaToolProperty.value).toBe(false);

    chartNode.showAreaToolProperty.value = true;
    expect(chartNode.showAreaToolProperty.value).toBe(true);
  });

  it("should constrain marker positions within valid range", () => {
    chartNode.showAreaToolProperty.value = true;

    // Try to set left marker beyond right marker
    const rightPos = chartNode["rightMarkerXProperty"].value;
    chartNode["leftMarkerXProperty"].value = rightPos + 1;

    // Should be constrained
    expect(chartNode["leftMarkerXProperty"].value).toBeLessThan(rightPos);
  });
});
```

---

## Performance Monitoring

### 8. Add Performance Measurement

**New File**: `src/common/utils/PerformanceMonitor.ts`

```typescript
import Logger from "./Logger.js";

interface PerformanceMetrics {
  count: number;
  totalTime: number;
  minTime: number;
  maxTime: number;
  avgTime: number;
}

/**
 * Performance monitoring utility for tracking operation durations.
 * Useful for identifying performance bottlenecks.
 */
export class PerformanceMonitor {
  private static metrics = new Map<string, number[]>();
  private static enabled = process.env.NODE_ENV === "development";

  /**
   * Enable or disable performance monitoring.
   */
  static setEnabled(enabled: boolean): void {
    this.enabled = enabled;
  }

  /**
   * Measure the execution time of a synchronous function.
   */
  static measure<T>(name: string, fn: () => T): T {
    if (!this.enabled) {
      return fn();
    }

    const start = performance.now();
    try {
      return fn();
    } finally {
      const duration = performance.now() - start;
      this.recordMeasurement(name, duration);
    }
  }

  /**
   * Measure the execution time of an async function.
   */
  static async measureAsync<T>(
    name: string,
    fn: () => Promise<T>,
  ): Promise<T> {
    if (!this.enabled) {
      return fn();
    }

    const start = performance.now();
    try {
      return await fn();
    } finally {
      const duration = performance.now() - start;
      this.recordMeasurement(name, duration);
    }
  }

  /**
   * Record a measurement manually.
   */
  static recordMeasurement(name: string, durationMs: number): void {
    if (!this.enabled) {
      return;
    }

    const measurements = this.metrics.get(name) || [];
    measurements.push(durationMs);
    this.metrics.set(name, measurements);

    // Warn if slow (>1 frame at 60fps)
    if (durationMs > 16.67) {
      Logger.warn(`Slow operation: ${name}`, {
        duration: `${durationMs.toFixed(2)}ms`,
      });
    }
  }

  /**
   * Get performance statistics for an operation.
   */
  static getStats(name: string): PerformanceMetrics | null {
    const measurements = this.metrics.get(name);
    if (!measurements || measurements.length === 0) {
      return null;
    }

    const totalTime = measurements.reduce((sum, t) => sum + t, 0);
    return {
      count: measurements.length,
      totalTime,
      minTime: Math.min(...measurements),
      maxTime: Math.max(...measurements),
      avgTime: totalTime / measurements.length,
    };
  }

  /**
   * Get all recorded metrics.
   */
  static getAllStats(): Record<string, PerformanceMetrics> {
    const stats: Record<string, PerformanceMetrics> = {};
    for (const [name, _measurements] of this.metrics) {
      const metric = this.getStats(name);
      if (metric) {
        stats[name] = metric;
      }
    }
    return stats;
  }

  /**
   * Clear all recorded metrics.
   */
  static clear(): void {
    this.metrics.clear();
  }

  /**
   * Log all performance statistics.
   */
  static logStats(): void {
    const stats = this.getAllStats();
    const entries = Object.entries(stats).sort(
      (a, b) => b[1].totalTime - a[1].totalTime,
    );

    if (entries.length === 0) {
      Logger.info("No performance metrics recorded");
      return;
    }

    Logger.info("=== Performance Statistics ===");
    for (const [name, metric] of entries) {
      Logger.info(`${name}:`, {
        calls: metric.count,
        total: `${metric.totalTime.toFixed(2)}ms`,
        avg: `${metric.avgTime.toFixed(2)}ms`,
        min: `${metric.minTime.toFixed(2)}ms`,
        max: `${metric.maxTime.toFixed(2)}ms`,
      });
    }
  }
}

export default PerformanceMonitor;
```

**Usage**:

```typescript
import PerformanceMonitor from "../utils/PerformanceMonitor.js";

// Measure solver performance
const result = PerformanceMonitor.measure("DVR-solve", () => {
  return solveDVR(potential, mass, numStates, gridConfig, energiesOnly);
});

// Log stats on demand
if (process.env.NODE_ENV === "development") {
  window.addEventListener("keydown", (e) => {
    if (e.key === "p" && e.ctrlKey) {
      PerformanceMonitor.logStats();
    }
  });
}
```

---

## Summary

These actionable suggestions are prioritized for:
1. **Immediate value** (Quick wins in 15min - 1 hour)
2. **Code quality** (Reduce duplication, improve maintainability)
3. **User experience** (Accessibility, error handling)
4. **Developer experience** (Testing, monitoring, logging)

Start with the quick wins to build momentum, then tackle the larger refactoring efforts.
