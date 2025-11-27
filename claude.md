# Claude Code Reference - QPPW (Quantum Particle Potential Wells)

This file contains important information about the QPPW codebase to help Claude Code work more effectively on this project.

## Project Overview

QPPW is a quantum mechanics simulation built with TypeScript and Vite. It demonstrates quantum particle behavior in various potential wells with interactive visualizations.

### Technology Stack
- **Language**: TypeScript
- **Build Tool**: Vite
- **UI Framework**: Custom SceneryStack (scenery, sun, scenery-phet)
- **State Management**: Axon (properties/observables)
- **Math**: DOT library

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

## Project Structure

```
src/
├── common/          # Shared code across all screens
│   ├── model/       # Core quantum mechanics models
│   ├── view/        # Shared UI components
│   └── utils/       # Utility functions
├── intro/           # Intro screen (simplified interface)
│   ├── model/
│   └── view/
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
└── QPPWColors.ts    # Color scheme definitions
```

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
  - Uses `IntroControlPanelNode` (with area tool)
  - Simpler interface for educational introduction

### Key Components
- **WaveFunctionChartNode**: Displays wave function/probability density
  - Has `showAreaToolProperty` for area measurement feature
  - Area tool only visible on intro screen

- **EnergyChartNode**: Displays energy levels and eigenstates

- **ControlPanelNode/IntroControlPanelNode**: Control panels
  - Different feature sets per screen type

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
  import type { SomeType } from "./module.js"
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

## Special Features

### Area Measurement Tool
- **Location**: Intro screen only
- **File**: `src/common/view/WaveFunctionChartNode.ts` (lines 78-1827)
- **Control**: `src/intro/view/IntroControlPanelNode.ts`
- **Functionality**:
  - Two draggable vertical markers
  - Calculates probability in region between markers
  - Uses trapezoidal integration
  - Works with both eigenstates and superpositions
  - Only enabled in probability density mode

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
new Text(stringManager.someStringProperty, options)
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

## Future Reference

When working on this codebase:
1. Always run lint and build before committing
2. Follow existing patterns for consistency
3. Use conditional rendering for optional features
4. Pass chart references only where needed
5. Test across all screens when making common changes
6. Document significant architectural changes in this file
