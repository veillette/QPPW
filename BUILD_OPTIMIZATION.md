# Build Optimization and Tree Shaking

## Overview

This document describes the tree shaking and build optimization strategies implemented in the QPPW project to reduce bundle size and improve load times.

## What is Tree Shaking?

Tree shaking is a term commonly used in JavaScript for dead-code elimination. It relies on the static structure of ES2015 module syntax (import/export) to remove unused code from the final bundle.

## Implemented Optimizations

### 1. Terser Minification

We use Terser for advanced minification with the following optimizations:

- **Dead code elimination**: Removes unreachable code
- **Unused variable removal**: Eliminates variables that are never used
- **Multiple optimization passes**: Runs 2 passes for better compression
- **Comment removal**: Strips all comments from production builds
- **Name mangling**: Shortens variable and function names

### 2. Rollup Tree Shaking

Vite uses Rollup under the hood, and we've configured it to:

- Use the "recommended" preset for tree shaking
- Respect `sideEffects` declarations in package.json
- Automatically detect and remove unused exports

### 3. Bundle Analysis

We've added the `rollup-plugin-visualizer` to help analyze what's in the bundle:

```bash
# Build with bundle analysis
npm run build:analyze

# The analysis report will be generated at dist/stats.html
# Open it in a browser to see a visual breakdown of the bundle
```

## Bundle Size Results

### Before Optimization
- vendor-scenery: 1,912.98 KB (gzip: 636.21 KB)
- vendor-other: 551.05 KB (gzip: 199.53 KB)
- **Total vendor size**: 2,463 KB (gzip: 835 KB)

### After Optimization
- vendor-scenery: 1,907.43 KB (gzip: 627.51 KB)
- vendor-other: 537.08 KB (gzip: 190.48 KB)
- **Total vendor size**: 2,444 KB (gzip: 817 KB)

### Savings
- **Raw size reduction**: ~19 KB (~0.8%)
- **Gzipped reduction**: ~18 KB (~2.2%)

## Why the Savings Are Modest

The improvements are modest because:

1. **SceneryStack is already well-optimized**: The scenerystack package is already minified and includes proper `sideEffects` declarations in its package.json
2. **No unused large dependencies**: Paper.js is NOT included in scenerystack (only `PaperAirplaneNode` component exists)
3. **Most code is actually used**: The simulation uses most of the imported scenerystack modules

## How to Further Reduce Bundle Size

If you need to reduce the bundle size further, consider:

### 1. Dynamic Imports
Load screens on demand instead of bundling everything upfront:

```typescript
// Instead of:
import IntroScreen from './intro/IntroScreen';

// Use:
const IntroScreen = () => import('./intro/IntroScreen');
```

### 2. Selective Imports
Import only what you need from scenerystack modules:

```typescript
// More specific imports help bundlers understand what's used
import { Node, Text } from 'scenerystack/scenery';
// vs
import * from 'scenerystack/scenery';
```

### 3. Code Splitting by Route
Split your application by routes or features so users only download what they need.

### 4. Analyze the Bundle
Use the bundle analyzer to identify large dependencies:

```bash
npm run build:analyze
open dist/stats.html
```

Look for:
- Large modules that could be lazy-loaded
- Duplicate dependencies
- Unused exports from large libraries

## Configuration Files

### vite.config.js
The main build configuration including:
- Terser minification settings
- Tree shaking configuration
- Manual chunk splitting
- Bundle visualization

### package.json
Added scripts:
- `build:analyze`: Build with bundle analysis enabled

## Additional Notes

### SceneryStack Side Effects
SceneryStack has many initialization side effects (like setting up globals, splash screens, etc.). The package.json properly declares these, so tree shaking respects them automatically.

### Chunk Splitting Strategy
The build is split into logical chunks:
- `vendor-scenery`: SceneryStack library code
- `vendor-katex`: KaTeX math rendering
- `vendor-other`: Other dependencies
- `quantum-core`: Core quantum simulation code
- `analytical`: Analytical solution algorithms
- `common-view`: Shared UI components
- `screen-*`: Individual screen implementations

This strategy allows:
- Better caching (vendor code rarely changes)
- Parallel loading
- Potential for lazy loading in the future

## Maintenance

When adding new dependencies:
1. Check if they're tree-shakeable (look for `sideEffects` in their package.json)
2. Prefer ES modules over CommonJS
3. Use specific imports rather than wildcard imports
4. Run `npm run build:analyze` to verify the impact on bundle size
