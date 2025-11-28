/**
 * Performance monitoring utility for tracking execution times of operations.
 * Collects metrics and warns about slow operations that exceed frame budget (16.67ms at 60fps).
 */

import Logger from "./Logger.js";

type PerformanceMetrics = {
  count: number;
  mean: number;
  max: number;
  min: number;
  total: number;
};

export class PerformanceMonitor {
  private static metrics: Map<string, number[]> = new Map();
  private static enabled: boolean = true;

  /**
   * Frame time threshold in milliseconds (60fps = 16.67ms per frame).
   * Operations exceeding this threshold will trigger a warning.
   */
  private static readonly FRAME_TIME_THRESHOLD = 16.67;

  /**
   * Enable or disable performance monitoring.
   */
  static setEnabled(enabled: boolean): void {
    PerformanceMonitor.enabled = enabled;
  }

  /**
   * Check if performance monitoring is enabled.
   */
  static isEnabled(): boolean {
    return PerformanceMonitor.enabled;
  }

  /**
   * Measure the execution time of a synchronous function.
   * Returns the function's result and logs a warning if execution exceeds frame time.
   *
   * @param name - Identifier for this operation (used in logs and metrics)
   * @param fn - Function to measure
   * @returns The result of the function
   *
   * @example
   * const result = PerformanceMonitor.measure('solve-schrodinger', () => {
   *   return solver.solveNumerical(potential, mass, numStates, grid);
   * });
   */
  static measure<T>(name: string, fn: () => T): T {
    if (!PerformanceMonitor.enabled) {
      return fn();
    }

    const start = performance.now();
    try {
      return fn();
    } finally {
      const duration = performance.now() - start;

      // Store the metric
      const existing = PerformanceMonitor.metrics.get(name) || [];
      existing.push(duration);
      PerformanceMonitor.metrics.set(name, existing);

      // Warn if operation is slow (exceeds one frame at 60fps)
      if (duration > PerformanceMonitor.FRAME_TIME_THRESHOLD) {
        Logger.warn(
          `Slow operation: ${name} took ${duration.toFixed(2)}ms (>${PerformanceMonitor.FRAME_TIME_THRESHOLD.toFixed(2)}ms threshold)`,
        );
      }
    }
  }

  /**
   * Measure the execution time of an async function.
   * Returns a promise that resolves to the function's result.
   *
   * @param name - Identifier for this operation
   * @param fn - Async function to measure
   * @returns Promise resolving to the function's result
   *
   * @example
   * const result = await PerformanceMonitor.measureAsync('load-data', async () => {
   *   return await fetch('/api/data').then(r => r.json());
   * });
   */
  static async measureAsync<T>(
    name: string,
    fn: () => Promise<T>,
  ): Promise<T> {
    if (!PerformanceMonitor.enabled) {
      return fn();
    }

    const start = performance.now();
    try {
      return await fn();
    } finally {
      const duration = performance.now() - start;

      // Store the metric
      const existing = PerformanceMonitor.metrics.get(name) || [];
      existing.push(duration);
      PerformanceMonitor.metrics.set(name, existing);

      // Warn if operation is slow
      if (duration > PerformanceMonitor.FRAME_TIME_THRESHOLD) {
        Logger.warn(
          `Slow async operation: ${name} took ${duration.toFixed(2)}ms (>${PerformanceMonitor.FRAME_TIME_THRESHOLD.toFixed(2)}ms threshold)`,
        );
      }
    }
  }

  /**
   * Get statistics for a specific operation.
   * Returns count, mean, max, min, and total time.
   *
   * @param name - Identifier for the operation
   * @returns Performance statistics or null if no data exists
   */
  static getStats(name: string): PerformanceMetrics | null {
    const values = PerformanceMonitor.metrics.get(name);
    if (!values || values.length === 0) {
      return null;
    }

    const total = values.reduce((a, b) => a + b, 0);
    return {
      count: values.length,
      mean: total / values.length,
      max: Math.max(...values),
      min: Math.min(...values),
      total,
    };
  }

  /**
   * Get statistics for all tracked operations.
   *
   * @returns Map of operation names to their statistics
   */
  static getAllStats(): Map<string, PerformanceMetrics> {
    const allStats = new Map<string, PerformanceMetrics>();

    for (const [name, values] of PerformanceMonitor.metrics.entries()) {
      if (values.length > 0) {
        const total = values.reduce((a, b) => a + b, 0);
        allStats.set(name, {
          count: values.length,
          mean: total / values.length,
          max: Math.max(...values),
          min: Math.min(...values),
          total,
        });
      }
    }

    return allStats;
  }

  /**
   * Clear metrics for a specific operation.
   *
   * @param name - Identifier for the operation
   */
  static clearStats(name: string): void {
    PerformanceMonitor.metrics.delete(name);
  }

  /**
   * Clear all collected metrics.
   */
  static clearAllStats(): void {
    PerformanceMonitor.metrics.clear();
  }

  /**
   * Log a summary of all performance metrics to the console.
   */
  static logSummary(): void {
    const allStats = PerformanceMonitor.getAllStats();

    if (allStats.size === 0) {
      Logger.info("No performance metrics collected yet.");
      return;
    }

    Logger.info("=== Performance Metrics Summary ===");
    for (const [name, stats] of allStats.entries()) {
      Logger.info(
        `${name}: count=${stats.count}, mean=${stats.mean.toFixed(2)}ms, ` +
          `min=${stats.min.toFixed(2)}ms, max=${stats.max.toFixed(2)}ms, ` +
          `total=${stats.total.toFixed(2)}ms`,
      );
    }
    Logger.info("===================================");
  }
}

export default PerformanceMonitor;
