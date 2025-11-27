# QPPW Architecture Guide for Claude

## Overview

QPPW (Quantum Particle in a Periodic Well) is a TypeScript-based physics simulation built using the PhET framework. This document outlines the architectural principles and conventions for AI assistants working on this codebase.

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
├── common/
│   ├── model/           # Shared physics and solver code
│   │   ├── BaseModel.ts          # Abstract base for all screen models
│   │   ├── Schrodinger1DSolver.ts # Core quantum solver
│   │   ├── NumerovSolver.ts      # Numerical method implementations
│   │   └── ...
│   └── view/            # Shared UI components
│       ├── BaseScreenView.ts     # Abstract base for all screen views
│       ├── WaveFunctionChartNode.ts # Wavefunction visualization
│       ├── EnergyChartNode.ts    # Energy level visualization
│       └── ...
├── intro/
│   ├── model/           # Intro screen model
│   └── view/            # Intro screen views
├── one-well/
│   ├── model/           # Single well model
│   └── view/            # Single well views
├── two-wells/
│   ├── model/           # Double well model
│   └── view/            # Double well views
└── many-wells/
    ├── model/           # Periodic potential model
    └── view/            # Periodic potential views
```

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
    const evolved = this.coefficients.map((c, i) =>
      c * Math.exp(-I * this.energies[i] * time)
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

## Testing and Quality

- **Linting**: Run `npm run lint` before committing
- **Building**: Run `npm run build` to verify TypeScript compilation
- **Separation verification**: If adding view code, ensure no physics calculations
- **Separation verification**: If adding model code, ensure no rendering logic

## Additional Resources

- PhET framework documentation: https://github.com/phetsims
- Axon (property system): Part of scenerystack
- See `BaseModel.ts` for core model patterns
- See `BaseScreenView.ts` for core view patterns

## Summary

The separation of concerns between view and model is the cornerstone of QPPW's architecture. When contributing:

1. **Models compute, views render**
2. **Keep physics in the model layer**
3. **Keep layout and rendering in the view layer**
4. **Follow the established patterns in Base classes**

This architecture makes the codebase maintainable, testable, and easier to reason about.
