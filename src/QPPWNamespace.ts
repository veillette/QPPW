/**
 * Main namespace file for the Quantum Bound States (QPPW) Simulation
 */

import type { Namespace } from "scenerystack";

// Check if we're in a Node.js environment
// @ts-expect-error - process is only available in Node.js
const isNode = typeof process !== 'undefined' && process.versions != null && process.versions.node != null;

// Conditional namespace creation based on environment
let qppw: Namespace;

if (!isNode) {
  // Browser/Vite environment - import from scenerystack
  const scenerystack = await import("scenerystack");
  const NamespaceConstructor = scenerystack.Namespace;
  qppw = new NamespaceConstructor("qppw");
} else {
  // Node.js environment (for running tests in terminal)
  // Create a mock namespace object that matches the Namespace interface
  qppw = {
    name: "qppw",
    register: <T>(_name: string, obj: T): T => {
      // No-op in Node.js environment - just return the object for registration purposes
      return obj;
    }
  } as Namespace;
}

export default qppw;
