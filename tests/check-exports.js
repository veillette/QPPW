#!/usr/bin/env node
// Set up browser globals
globalThis.self = globalThis;
globalThis.window = globalThis;
globalThis.document = { createElement: () => ({ getContext: () => null, style: {} }), body: { appendChild: () => {} } };
globalThis.phet = globalThis.phet || {};
globalThis.phet.chipper = globalThis.phet.chipper || {};
globalThis.phet.joist = globalThis.phet.joist || {};
globalThis.phet.chipper.packageObject = { name: 'scenerystack' };
globalThis.requestAnimationFrame = (cb) => setTimeout(cb, 16);

// Try to import and check exports
try {
  const dotModule = await import('scenerystack/dot');
  console.log('Exports from scenerystack/dot:');
  console.log(Object.keys(dotModule).sort().join('\n'));
} catch (e) {
  console.error('Error importing scenerystack/dot:', e.message);
}
