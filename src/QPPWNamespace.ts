/**
 * Main namespace file for the Quantum Bound States (QPPW) Simulation
 *
 * This provides a simple namespace implementation compatible with
 * scenerystack's Namespace API for registration purposes.
 */

/**
 * Simple namespace class for organizing and registering objects.
 * API compatible with scenerystack/phet-core Namespace.
 */
class Namespace {
  public name: string;
  protected registry: Map<string, unknown>;

  constructor(name: string) {
    this.name = name;
    this.registry = new Map();
  }

  /**
   * Register an object with this namespace.
   *
   * @param key - The name to register under
   * @param value - The object to register
   * @returns The registered value
   */
  register<T>(key: string, value: T): T {
    this.registry.set(key, value);
    return value;
  }

  /**
   * Get a registered object by name.
   *
   * @param key - The name to look up
   * @returns The registered value, or undefined
   */
  get(key: string): unknown {
    return this.registry.get(key);
  }

  /**
   * Get the namespace name.
   */
  getName(): string {
    return this.name;
  }
}

// Create and export the namespace object
const qppw = new Namespace("qppw");

export default qppw;
