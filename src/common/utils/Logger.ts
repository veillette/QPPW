/**
 * Simple logging utility for the quantum mechanics simulation.
 * Provides methods for different log levels: debug, info, warn, error.
 */

export enum LogLevel {
  DEBUG = 0,
  INFO = 1,
  WARN = 2,
  ERROR = 3,
  NONE = 4,
}

export class Logger {
  private static currentLevel: LogLevel = LogLevel.INFO;

  /**
   * Set the minimum log level. Messages below this level won't be displayed.
   */
  static setLevel(level: LogLevel): void {
    Logger.currentLevel = level;
  }

  /**
   * Get the current log level.
   */
  static getLevel(): LogLevel {
    return Logger.currentLevel;
  }

  /**
   * Log a debug message.
   */
  static debug(message: string, ...args: unknown[]): void {
    if (Logger.currentLevel <= LogLevel.DEBUG) {
      console.log(`[DEBUG] ${message}`, ...args);
    }
  }

  /**
   * Log an info message.
   */
  static info(message: string, ...args: unknown[]): void {
    if (Logger.currentLevel <= LogLevel.INFO) {
      console.log(`[INFO] ${message}`, ...args);
    }
  }

  /**
   * Log a warning message.
   */
  static warn(message: string, ...args: unknown[]): void {
    if (Logger.currentLevel <= LogLevel.WARN) {
      console.warn(`[WARN] ${message}`, ...args);
    }
  }

  /**
   * Log an error message.
   */
  static error(message: string, ...args: unknown[]): void {
    if (Logger.currentLevel <= LogLevel.ERROR) {
      console.error(`[ERROR] ${message}`, ...args);
    }
  }
}

export default Logger;
