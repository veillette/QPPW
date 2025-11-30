# QPPW Architecture Guide for Claude

## Overview

QPPW (Quantum Particle in a Periodic Well) is a TypeScript-based physics simulation built using the PhET framework. This document outlines the architectural principles and conventions for AI assistants working on this codebase.

### Technology Stack

- **Language**: TypeScript
- **Build Tool**: Vite
- **UI Framework**: Custom SceneryStack (scenery, sun, scenery-phet)
- **State Management**: Axon (properties/observables)
- **Math**: DOT library
- **Chart Library**: Bamboo (ChartTransform, ChartRectangle, AxisLine, TickMarkSet, TickLabelSet)

## Core Architectural Principle: Separation of Concerns

The codebase maintains a **strict separation between view and model layers**. This separation is fundamental to the architecture and must be preserved in all contributions.

### Model Layer Responsibilities

The model layer (`src/*/model/`) is responsible for:

- **All physics calculations** including:
  - Quantum mechanics computations (wavefunction evolution, probability densities, energy levels)
  - Numerical solving (Numerov, spectral methods, DVR)
  - Time evolution of quantum states
  - Interpolation and integration of physical quantities
- **State management** using Axon properties
- **Simulation logic** (stepping, resetting, configuration)
- **Data provision** for visualization

**Key principle**: The model performs all calculations and exposes results through well-defined interfaces. Views request computed data from the model rather than performing calculations themselves.

### View Layer Responsibilities

The view layer (`src/*/view/`) is responsible for:

- **Rendering** visual representations of model state
- **Layout** and positioning of UI components
- **User interaction** handling (forwarding to model)
- **Observing model properties** and updating displays accordingly

**Key principle**: Views are purely presentational. They transform data from the model into visual representations but never perform physics calculations or contain business logic.

## Directory Structure

```
src/
├── common/          # Shared code across all screens
│   ├── model/       # Shared physics and solver code
│   │   ├── BaseModel.ts                # Abstract base for all screen models
│   │   ├── Schrodinger1DSolver.ts      # Core quantum solver
│   │   ├── NumerovSolver.ts            # Standard Numerov method
│   │   ├── MatrixNumerovSolver.ts      # Matrix-based Numerov
│   │   ├── WavefunctionNumerovSolver.ts # Wavefunction-specific Numerov
│   │   ├── SpectralSolver.ts           # Spectral method solver
│   │   ├── FGHSolver.ts                # Fourier Grid Hamiltonian solver
│   │   ├── DVRSolver.ts                # Discrete Variable Representation solver
│   │   ├── PotentialFunction.ts        # Potential type definitions
│   │   ├── PotentialFactory.ts         # Factory for creating potentials
│   │   ├── WavefunctionStandardization.ts # Normalization utilities
│   │   ├── LinearAlgebraUtils.ts       # Matrix operations
│   │   ├── QuantumConstants.ts         # Physical constants
│   │   ├── analytical-solutions/       # Analytical solutions for standard potentials
│   │   │   ├── infinite-well-potential.ts
│   │   │   ├── harmonic-oscillator.ts
│   │   │   ├── finite-well-potential.ts
│   │   │   ├── morse-potential.ts
│   │   │   ├── poschl-teller-potential.ts
│   │   │   ├── rosen-morse-potential.ts
│   │   │   ├── eckart-potential.ts
│   │   │   ├── coulomb-1d-potential.ts
│   │   │   ├── coulomb-3d-potential.ts
│   │   │   ├── triangular-potential.ts
│   │   │   └── asymmetric-triangle-potential.ts
│   │   └── ...
│   ├── view/        # Shared UI components
│   │   ├── BaseScreenView.ts           # Abstract base for all screen views
│   │   ├── BaseChartNode.ts            # Base class for all charts
│   │   ├── WaveFunctionChartNode.ts    # Wavefunction/probability visualization
│   │   ├── EnergyChartNode.ts          # Energy level visualization
│   │   ├── WavenumberChartNode.ts      # Fourier transform (momentum) visualization
│   │   ├── ControlPanelNode.ts         # Main control panel
│   │   ├── SimulationControlBar.ts     # Play/pause/step controls
│   │   ├── SuperpositionDialog.ts      # Superposition state builder dialog
│   │   ├── RMSIndicatorUtils.ts        # RMS calculation utilities
│   │   ├── ChartConstants.ts           # Shared chart configuration
│   │   └── ...
│   └── utils/       # Utility functions
├── intro/           # Intro screen (simplified educational interface)
│   ├── model/
│   │   └── IntroModel.ts
│   └── view/
│       ├── IntroScreenView.ts
│       ├── IntroControlPanelNode.ts    # Simplified controls with extra tools
│       └── IntroViewState.ts
├── one-well/        # Single quantum well screen
│   ├── model/
│   └── view/
├── two-wells/       # Double quantum well screen
│   ├── model/
│   └── view/
├── many-wells/      # Multiple quantum wells screen
│   ├── model/
│   └── view/
├── i18n/            # Internationalization/strings
│   └── StringManager.ts
└── QPPWColors.ts    # Color scheme definitions (all use properties for theme support)
```

**Key Subdirectories**:

- `analytical-solutions/` - Each potential type has analytical solutions for eigenstates and classical probability distributions
- `view/` - All chart nodes extend `BaseChartNode` for consistency
- `intro/view/` - Contains specialized educational tools not present in other screens

## Historical Context: The View/Model Refactoring

In commit `3f30d06`, a significant refactoring established clear separation between view and model layers:

### What Changed

**Before**: `WaveFunctionChartNode` (view) contained ~200 lines of physics calculations:

- Time evolution of quantum superpositions
- Probability density integration
- Wavefunction interpolation and derivatives

**After**: These calculations moved to `BaseModel`, and `WaveFunctionChartNode` now:

- Calls `model.getTimeEvolvedSuperposition()` for quantum state evolution
- Calls `model.getProbabilityInRegion()` for probability calculations
- Calls `model.getWavefunctionAtPosition()` for interpolated values

### Benefits Achieved

1. **Eliminated duplication**: ~200 lines of physics code removed from view layer
2. **Better testability**: Physics logic centralized and easier to unit test
3. **Clearer responsibilities**: Views focus solely on rendering
4. **Easier maintenance**: Changes to physics calculations happen in one place
5. **Better encapsulation**: Quantum mechanics implementation details hidden from views

## Guidelines for Contributors

### DO ✅

- Keep all physics calculations in model classes
- Have views request computed data from models via method calls
- Use Axon properties for reactive state management
- Maintain the existing directory structure (model/ and view/ separation)
- Follow the established pattern: views observe model properties and render accordingly
- Centralize shared logic in `BaseModel` and `BaseScreenView`

### DON'T ❌

- Put physics calculations in view classes
- Have models directly manipulate view elements
- Duplicate calculation logic between model and view
- Mix rendering logic with physics computations
- Create circular dependencies between view and model

### Example: Correct Separation

```typescript
// ❌ WRONG: View performing physics calculation
class MyChartNode extends Node {
  private updateWavefunction() {
    // DON'T: Calculate physics in the view
    const evolved = this.coefficients.map(
      (c, i) => c * Math.exp(-I * this.energies[i] * time),
    );
  }
}

// ✅ CORRECT: View requesting data from model
class MyChartNode extends Node {
  private updateWavefunction() {
    // DO: Request computed data from model
    const evolved = this.model.getTimeEvolvedSuperposition(time);
    this.renderWavefunction(evolved);
  }
}
```

## Key Model Methods for View Consumption

The model layer exposes these methods for view consumption (see `BaseModel.ts`):

- `getTimeEvolvedSuperposition(time)`: Returns evolved quantum state
- `getProbabilityInRegion(xMin, xMax)`: Returns integrated probability
- `getWavefunctionAtPosition(x)`: Returns interpolated wavefunction values
- Various property getters for energies, coefficients, potential parameters, etc.

## CRITICAL BUGS & LEARNINGS (Read This First!)

### ⚠️ Unit Conversion Bug (Fixed in 8bf740d)

**The Problem**: Converting wavefunctions from meters to nanometers for display **MUST preserve normalization**.

**WRONG**:

```typescript
// ❌ This BREAKS normalization - creates astronomical values (~10^18)
const psiInNm = psi * Math.sqrt(M_TO_NM);
```

**CORRECT**:

```typescript
// ✅ Divide by sqrt(M_TO_NM) to preserve ∫|ψ|² dx = 1
const psiInNm = psi / Math.sqrt(M_TO_NM);
```

**Why**: When converting position units from m → nm (multiply by 10^9), wavefunction amplitude must be divided by √(10^9) to preserve the normalization condition ∫|ψ|² dx = 1.

**Location**: `BaseModel.ts:getWavefunctionInNmUnits()`, `getSuperpositionInNmUnits()`

### ⚠️ Async Initialization Pattern (Fixed in 8bf740d)

**The Problem**: Chart nodes calling `update()` synchronously in constructors can **hang the entire simulation**.

**WRONG**:

```typescript
// ❌ Synchronous property.link() causes immediate execution
this.someProperty.link((value) => {
  this.update(); // Can cause infinite loops or hangs
});
```

**CORRECT**:

```typescript
// ✅ Use lazyLink() to prevent immediate execution
this.someProperty.lazyLink((value) => {
  this.update();
});

// ✅ Defer initial update to next tick
setTimeout(() => {
  this.update();
}, 0);
```

**Why**: Property links fire immediately on attachment. If multiple charts cross-trigger updates during construction, you get infinite loops or page hangs.

**Locations**: `WaveFunctionChartNode.ts`, `EnergyChartNode.ts`, `WavenumberChartNode.ts`

**Add safety checks**:

```typescript
// Prevent runaway updates
if (this.updateCounter++ > 10) {
  console.warn("Infinite update loop detected");
  return;
}

// Validate chart ranges
if (Math.abs(yMax - yMin) > 1000) {
  console.warn("Invalid chart range detected");
  return;
}
```

### ⚠️ Classical Probability Distribution Singularities (Fixed in 4d403b6)

**The Problem**: Classical probability distributions P(x) = 1/√(E - V(x)) go to **infinity at turning points** where E = V(x).

**Solution**: Use epsilon-based regularization with **relative threshold** (1% of max kinetic energy):

```typescript
// ✅ CORRECT: Relative epsilon scales with energy
const maxKE = energy - minPotential;
const epsilon = 0.01 * maxKE; // 1% of max kinetic energy
const probability = 1 / Math.sqrt(Math.max(energy - V(x), epsilon));
```

**WRONG**:

```typescript
// ❌ Absolute epsilon doesn't scale properly
const epsilon = 1e-10; // Too small, still causes infinities
```

**Special Case - Harmonic Oscillator**: The analytical formula P(x) = 1/(π√(A²-x²)) is **unusable** due to singularities. Use the epsilon method instead.

**Affected Files**: All analytical solution files in `src/common/model/analytical-solutions/`

### ⚠️ Normalization: Numerical > Analytical (Fixed in 873feb1, 74a212b)

**Key Learning**: For complex potentials (Pöschl-Teller, 3D Coulomb), **numerical normalization is more reliable** than analytical formulas involving Gamma functions and special polynomial orthogonality.

**Pattern**:

```typescript
// ✅ CORRECT: Numerical integration is robust
const norm = Math.sqrt(
  trapezoidalIntegration(
    x,
    psi.map((p) => p * p),
  ),
);
const normalizedPsi = psi.map((p) => p / norm);
```

**Why analytical formulas fail**:

- Complex Gamma function evaluations accumulate errors
- Grid discretization breaks orthogonality assumptions
- 3D → 1D radial conversion is non-trivial

**Successfully fixed**:

- Pöschl-Teller: Analytical formula had 2-5x errors → now 1.000000
- 3D Coulomb: Radial normalization gave ~10^20 errors → now 1.000000
- Harmonic oscillator, finite well: Multiple fixes needed

## Essential Commands

### Building

```bash
npm run build
```

- Runs TypeScript type checking (`tsc --noEmit`)
- Builds production bundle with Vite
- **Always run this before committing** to ensure no type errors

### Linting & Formatting

```bash
npm run lint
```

- Runs ESLint on the entire codebase
- **Must pass before committing**
- ESLint automatically applies formatting rules (similar to Prettier)
- No separate format command needed - linting handles formatting

### Development

```bash
npm start
```

- Starts Vite dev server
- Hot module replacement enabled
- Available at http://localhost:5173 (typically)

### Type Checking

```bash
npx tsc --noEmit
```

- Type-checks without emitting files
- Useful for quick type validation

### Testing

```bash
# Run main comprehensive test suite
npm test

# Run specific test suites
npm run test:coulomb        # Coulomb potential verification
npm run test:double-well    # Double well tests (23 stringent tests)
npm run test:comprehensive  # Wavefunction tests (7 methods × 6 potentials)
```

**Test Files** (`tests/`):

- `AccuracyTests.ts` - Main comprehensive test suite
- `test-wavefunction-comprehensive.ts` - 7 solving methods × 6 potentials
- `test-double-well.ts` - 23 stringent double well tests
- `verify-coulomb.ts` - Coulomb verification (use with npm run test:coulomb)
- `test-coulomb-wrapper.ts` - Coulomb wrapper validation
- **Note**: Vestigial debugging files removed in commit 0e737aa

## Code Architecture & Patterns

### Screen Structure

Each screen follows this pattern:

1. **Model**: Extends base model, manages state and physics
2. **ScreenView**: Renders charts and UI
3. **ControlPanelNode**: Control panel with settings

### Two Control Panel Types

- **ControlPanelNode** (`src/common/view/ControlPanelNode.ts`)
  - Used by: OneWellScreenView, TwoWellsScreenView, ManyWellsScreenView
  - Features: Full controls including superposition, phase color, wave function views

- **IntroControlPanelNode** (`src/intro/view/IntroControlPanelNode.ts`)
  - Used by: IntroScreenView
  - Features: Simplified controls, no superposition or phase color options
  - **Has measure area tool** (other screens do not)

### Screen View Patterns

- **BaseScreenView**: Base class with `createStandardLayout()` method
  - Creates `WaveFunctionChartNode` and `EnergyChartNode`
  - Creates `ControlPanelNode` (without area tool)
  - Used by OneWell, TwoWells, ManyWells screens

- **IntroScreenView**: Custom implementation
  - Creates charts manually
  - Uses `IntroControlPanelNode` (with specialized features)
  - Simpler interface for educational introduction

### Chart Types (All extend BaseChartNode)

1. **WaveFunctionChartNode** - Displays wavefunction/probability density
   - Modes: probability density, wavefunction (real/imaginary), phase color
   - Features:
     - Area measurement tool (intro screen only)
     - RMS position (Δx) indicator with double-headed arrow
     - Average position vertical line indicator (<x>)
     - Derivative tool (dψ/dx) with tracking circle
     - Curvature tool (d²ψ/dx²) with tracking circle
   - File: `src/common/view/WaveFunctionChartNode.ts`

2. **EnergyChartNode** - Displays energy levels and eigenstates
   - Shows potential well shape
   - Displays quantized energy levels as horizontal lines
   - Interactive eigenstate selection
   - File: `src/common/view/EnergyChartNode.ts`

3. **WavenumberChartNode** - Displays Fourier transform |φ(k)|²
   - Shows momentum distribution of quantum state
   - Features:
     - Average wavenumber <k> indicator
     - RMS wavenumber k_rms indicator with double-headed arrow
   - Only visible when enabled via control panel
   - File: `src/common/view/WavenumberChartNode.ts`

### Chart Utilities

- **RMSIndicatorUtils** (`src/common/view/RMSIndicatorUtils.ts`)
  - `createDoubleArrowShape(x1, x2, y)` - Creates double-headed arrow for RMS indicators
  - `calculateRMSStatistics(xData, yData)` - Computes average and RMS spread via trapezoidal integration
  - Used by both WaveFunctionChartNode and WavenumberChartNode

### Control Panel Types

- **ControlPanelNode** (`src/common/view/ControlPanelNode.ts`)
  - Used by: OneWellScreenView, TwoWellsScreenView, ManyWellsScreenView
  - Features: Full controls including superposition, phase color, wave function views

- **IntroControlPanelNode** (`src/intro/view/IntroControlPanelNode.ts`)
  - Used by: IntroScreenView
  - Features:
    - Simplified controls (no superposition or phase color options)
    - Area measurement tool toggle (unique to intro)
    - Particle mass slider (unique to intro)
    - RMS indicator toggle (unique to intro)
    - Derivative/curvature tool toggles (unique to intro)

## Important Conventions

### File Naming

- TypeScript files use PascalCase: `WaveFunctionChartNode.ts`
- Test files: `*.test.ts`
- Import with `.js` extension: `import { Foo } from "./Foo.js"`
  - TypeScript requires .js even though source files are .ts

### Code Style

- **Indentation**: 2 spaces
- **Quotes**: Double quotes for strings
- **Semicolons**: Required
- **Line Length**: Managed by ESLint
- **Comments**: JSDoc style for public methods

### TypeScript

- Strict mode enabled
- Explicit types preferred
- Use type imports when importing only types:
  ```typescript
  import type { SomeType } from "./module.js";
  ```

### Property Pattern (Axon)

The codebase uses Axon properties extensively:

```typescript
// Create property
public readonly someProperty: Property<number>;

// In constructor
this.someProperty = new Property(initialValue);

// Link to property changes
this.someProperty.link((value) => {
  // React to changes
});
```

### Conditional UI Elements

When adding optional UI elements (checkboxes, controls):

```typescript
// 1. Create content conditionally
const checkboxContent = condition
  ? new Checkbox(property, label, options)
  : null;

// 2. Wrap in Node if exists
const checkbox = checkboxContent
  ? new Node({ children: [checkboxContent], x: 20 })
  : null;

// 3. Add to children array conditionally
if (checkbox) {
  children.push(checkbox);
}
```

## Common Tasks

### Adding a Feature to a Specific Screen Only

Example: Area tool is only on intro screen

1. **If feature needs chart reference:**
   - Pass chart to control panel in ScreenView constructor
   - Update control panel to accept optional parameter

2. **Add feature to control panel:**
   - Create UI elements conditionally based on presence of chart
   - Add enable/disable logic based on mode
   - Add to children array conditionally

3. **Ensure other screens don't have it:**
   - Don't pass chart reference in other screen views
   - Feature won't appear due to conditional rendering

### Modifying Display Modes

Display modes are managed in model:

- `probabilityDensity`: Shows probability density
- `waveFunction`: Shows real/imaginary parts
- `phaseColor`: Color-coded phase (not in intro screen)

Enable/disable controls based on mode:

```typescript
this.model.displayModeProperty.link((mode: string) => {
  checkbox.enabled = mode === "probabilityDensity";
});
```

### Working with Potential Wells

Potential types defined in `src/common/model/PotentialFunction.ts`:

- INFINITE_WELL
- FINITE_WELL
- HARMONIC_OSCILLATOR
- MORSE
- POSCHL_TELLER
- ROSEN_MORSE
- ECKART
- ASYMMETRIC_TRIANGLE
- TRIANGULAR
- COULOMB_1D
- COULOMB_3D

Each potential has specific parameters (width, depth, barrier height, offset).

## Git Workflow

### Branches

- **main**: Primary development branch
- Use descriptive branch names for features

### Before Committing

1. Run `npm run build` - ensure no type errors
2. Run `npm run lint` - ensure code style compliance
3. Test functionality manually if UI changes

### Commit Messages

- Use conventional commits style
- Be descriptive but concise
- Example: "Add area measurement tool to intro screen only"

## Special Features (Intro Screen)

The intro screen has several unique educational tools not present on other screens:

### 1. Area Measurement Tool

- **Location**: Intro screen only
- **File**: `WaveFunctionChartNode.ts` (conditional rendering)
- **Control**: `IntroControlPanelNode.ts` - "Show Measure Area" checkbox
- **Functionality**:
  - Two draggable vertical markers (left and right)
  - Real-time probability calculation: ∫|ψ(x)|² dx between markers
  - Uses trapezoidal integration (`model.getProbabilityInRegion()`)
  - Works with both eigenstates and superpositions
  - Only enabled in probability density mode
  - Displays probability as percentage with "Probability in shaded region" label

### 2. RMS Position Indicator (Δx)

- **Control**: "Show RMS Indicator" checkbox
- **Visual**: Double-headed horizontal arrow showing spread
- **Calculation**: Uses `calculateRMSStatistics()` from `RMSIndicatorUtils.ts`
- **Display**: Shows average position <x> (vertical line) and RMS spread Δx (arrow)
- **Only enabled in probability density mode**

### 3. Derivative Tool (dψ/dx)

- **Control**: Checkbox in intro control panel
- **Visual**: Plots first derivative with tracking circle
- **Important**: Derivative scaling uses correct unit conversion (fixed in commit 33dc948)
- **Pattern**: Falls back to finite differences if analytical derivative unavailable

### 4. Curvature Tool (d²ψ/dx²)

- **Control**: Checkbox in intro control panel
- **Visual**: Plots second derivative with tracking circle
- **Uses**: Shows curvature of wavefunction (related to kinetic energy via Schrödinger equation)

### 5. Particle Mass Slider

- **Added**: Commit 7142b58
- **Location**: Intro screen only
- **Range**: Allows users to vary particle mass and observe effects on energy levels
- **Effect**: Heavier particles → more closely spaced energy levels

### Chart Constants

- Defined in `src/common/view/ChartConstants.ts`
- Shared margins and ranges across charts
- Use these constants for consistency

## Debugging Tips

### Type Errors

- Run `npx tsc --noEmit` for detailed type information
- Check import paths use `.js` extension
- Verify property types match Axon patterns

### Build Errors

- Clear dist folder: `rm -rf dist`
- Reinstall dependencies: `npm install`
- Check Vite config if module resolution issues

### Runtime Errors

- Use browser DevTools
- Check property links for memory leaks
- Verify dispose methods called for dynamic components

## String Management

Strings are managed through `StringManager` in `src/i18n/StringManager.ts`:

```typescript
import stringManager from "../../i18n/StringManager.js";

// Use string properties
new Text(stringManager.someStringProperty, options);
```

This supports internationalization (i18n).

## Colors

All colors defined in `src/QPPWColors.ts` as properties:

```typescript
import QPPWColors from "../../QPPWColors.js";

{
  fill: QPPWColors.textFillProperty,
  stroke: QPPWColors.panelStrokeProperty
}
```

Use properties (not raw values) for theme support.

## Testing

Currently no automated test suite configured.
Manual testing workflow:

1. Run `npm start`
2. Navigate to each screen
3. Test features in different modes
4. Verify controls enable/disable correctly

## Performance Considerations

- Wave function calculations can be expensive
- Use requestAnimationFrame for animations
- Unlink property listeners in dispose methods
- Avoid creating new objects in animation loops

## Testing and Quality

- **Linting**: Run `npm run lint` before committing
- **Building**: Run `npm run build` to verify TypeScript compilation
- **Separation verification**: If adding view code, ensure no physics calculations
- **Separation verification**: If adding model code, ensure no rendering logic

## Documentation

Recent documentation additions (commits d945199, e118ce5, 5bc49e3):

- **`doc/`** directory (renamed from `docs/` to match PhET convention)
  - `teacher-guide.md` - Comprehensive guide for educators
  - `implementation-notes.md` - Technical implementation details
  - `accessibility-plan.md` - Accessibility implementation roadmap
  - `reports/` - Performance and analysis reports

## Recent Major Improvements (Last 3 Months)

### Bug Fixes

1. **Unit conversion normalization** (8bf740d) - Critical fix preventing astronomical probability values
2. **Async initialization** (8bf740d) - Prevents simulation hangs from chart update loops
3. **Classical probability singularities** (4d403b6) - Proper handling of turning point infinities
4. **Wavefunction normalization** (873feb1, 74a212b) - Fixed Pöschl-Teller and 3D Coulomb normalization
5. **Derivative scaling** (33dc948) - Corrected unit conversion in derivative calculations

### Features Added

1. **Wavenumber chart** (12a0ed0, 75a18ae) - Fourier transform visualization with RMS indicators
2. **RMS position indicator** (c8d18b5, 271828a, 04f06ac) - Shows quantum uncertainty Δx
3. **Particle mass slider** (4756245, 7142b58) - Educational tool for intro screen
4. **Derivative tools** (ff7132d, 1a1d2d0) - First and second derivative visualization
5. **Position tracking circles** (1a1d2d0) - Visual feedback for derivative tools

### Refactoring

1. **Test consolidation** (0e737aa) - Removed redundant test files, kept comprehensive suites
2. **Prettier integration** - Automated code formatting in build pipeline
3. **Documentation reorganization** (d945199) - Migrated to `doc/` directory

## Solver Methods

The simulation supports **7 different numerical solving methods** (tested in `test-wavefunction-comprehensive.ts`):

1. **Numerov** - Classic shooting method (good for most potentials)
2. **Matrix Numerov** - Matrix-based eigenvalue approach
3. **Wavefunction Numerov** - Direct wavefunction evolution
4. **Spectral** - Fourier-based spectral method
5. **FGH (Fourier Grid Hamiltonian)** - Hybrid Fourier/grid approach
6. **DVR (Discrete Variable Representation)** - Variational method
7. **Analytical** - Exact solutions for standard potentials (when available)

**Recommendation**: Numerov method is the default and most reliable for general use.

## Additional Resources

- PhET framework documentation: https://github.com/phetsims
- Axon (property system): Part of scenerystack
- Bamboo (charting library): Part of scenerystack
- See `BaseModel.ts` for core model patterns
- See `BaseScreenView.ts` for core view patterns
- See `BaseChartNode.ts` for chart implementation patterns
- See `RMSIndicatorUtils.ts` for statistical calculation patterns

## Summary

The separation of concerns between view and model is the cornerstone of QPPW's architecture. When contributing:

1. **Models compute, views render** - Never mix physics calculations into view code
2. **Keep physics in the model layer** - All quantum mechanics, integration, and calculations belong in model
3. **Keep layout and rendering in the view layer** - Views transform model data into visual representations
4. **Follow the established patterns in Base classes** - Extend BaseModel, BaseScreenView, BaseChartNode
5. **Watch for critical bugs** - Unit conversions, async initialization, normalization, singularities
6. **Use numerical normalization** - More reliable than analytical formulas for complex potentials
7. **Test comprehensively** - Run test suites before committing major physics changes

This architecture makes the codebase maintainable, testable, and easier to reason about.

## Quick Reference: Common Pitfalls

❌ **DON'T**:

- Multiply wavefunctions by `sqrt(M_TO_NM)` when converting to nm
- Use `property.link()` in chart constructors (causes immediate execution)
- Use absolute epsilon for classical probability singularities
- Trust analytical normalization formulas for complex potentials (use numerical)
- Mix physics calculations into view classes

✅ **DO**:

- Divide wavefunctions by `sqrt(M_TO_NM)` to preserve normalization
- Use `property.lazyLink()` and defer initial updates with `setTimeout(..., 0)`
- Use relative epsilon (1% of max KE) for singularity handling
- Validate normalization with `∫|ψ|² dx ≈ 1.0` after calculations
- Keep all calculations in model, expose via methods like `getTimeEvolvedSuperposition()`
