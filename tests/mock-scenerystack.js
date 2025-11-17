/**
 * Mock implementation of scenerystack Namespace for Node.js testing
 * This allows tests to run in Node.js without the full scenerystack library
 */

export class Namespace {
  constructor(name) {
    this.name = name;
  }

  register(name, obj) {
    // Simply return the object for registration purposes
    return obj;
  }
}

// Export as default to match the real scenerystack module structure
export default { Namespace };
