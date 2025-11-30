# QPPW Accessibility Implementation Plan

**Date:** November 30, 2025
**Version:** 1.0
**Status:** Planning Phase
**Target Framework:** PhET Scenery Stack PDOM (Parallel DOM)

---

## Table of Contents

1. [Executive Summary](#executive-summary)
2. [Current State Assessment](#current-state-assessment)
3. [Accessibility Architecture Overview](#accessibility-architecture-overview)
4. [Implementation Phases](#implementation-phases)
5. [Detailed Component Specifications](#detailed-component-specifications)
6. [Testing & Validation Strategy](#testing--validation-strategy)
7. [Resources & References](#resources--references)
8. [Appendix](#appendix)

---

## Executive Summary

### Project Overview

**QPPW (Quantum Particle in a Periodic Well)** is an interactive educational physics simulation built with SceneryStack 3.0 that visualizes quantum bound states in various potential wells. The simulation currently provides rich visual interactions but lacks accessibility features for users with disabilities.

### Current Accessibility Status

- ❌ **No PDOM (Parallel DOM)** implementation
- ❌ **No screen reader support**
- ❌ **No keyboard navigation**
- ❌ **No ARIA attributes**
- ❌ **No alternative text for visualizations**
- ❌ **Pure canvas rendering** with no semantic HTML

### Proposed Solution

Implement comprehensive accessibility using PhET's Parallel DOM (PDOM) architecture, enabling:

- ✅ Full keyboard navigation
- ✅ Screen reader compatibility (NVDA, JAWS, VoiceOver, TalkBack)
- ✅ Dynamic content descriptions
- ✅ Live announcements for state changes
- ✅ Accessible interactive controls
- ✅ Semantic HTML structure alongside canvas rendering

### Implementation Scope

- **4 Screens:** Intro, One Well, Two Wells, Many Wells
- **20+ Interactive Components:** Sliders, buttons, checkboxes, combo boxes
- **6 Interactive Tools:** Area measurement, curvature, derivative, zeros, phase, classical probability
- **3 Primary Charts:** Energy chart, wavefunction chart, wavenumber chart
- **12+ Potential Types:** Each requiring unique accessible descriptions

---

## Current State Assessment

### Project Architecture

**Technology Stack:**

- **Framework:** SceneryStack 3.0.0 (PhET graphics framework)
- **Language:** TypeScript (strict mode)
- **Build Tool:** Vite
- **UI Libraries:**
  - Scenery (graphics engine)
  - Sun (UI widgets)
  - Scenery-Phet (physics widgets)
  - Bamboo (charting)
  - Axon (reactive properties)

**Codebase Metrics:**

- 106 TypeScript source files
- Clean model-view separation
- Internationalization support (English, French)
- Comprehensive test suite for physics accuracy

### User Interface Components

#### Interactive Controls

| Component             | Count                  | Current State   | Accessibility Need                       |
| --------------------- | ---------------------- | --------------- | ---------------------------------------- |
| Sliders               | 6-8 per screen         | Mouse-only      | Keyboard control, value announcements    |
| Buttons               | 3-5 per screen         | Visual only     | ARIA roles, focus management             |
| Checkboxes            | 5-8 per screen         | Canvas-rendered | Semantic HTML, checked state             |
| Radio Buttons         | 1 group                | Canvas-rendered | Proper grouping, selection announcements |
| ComboBox              | 1 (potential selector) | Visual dropdown | Accessible dropdown, navigation          |
| Energy Level Selector | 1-12 levels            | Click-only      | Keyboard selection, state announcements  |

#### Visualization Components

| Component           | Purpose                                              | Accessibility Challenge                              |
| ------------------- | ---------------------------------------------------- | ---------------------------------------------------- |
| Energy Chart        | Shows potential energy landscape and discrete levels | Needs textual summary of energy spectrum             |
| Wavefunction Chart  | Displays ψ(x) in various modes                       | Requires statistical summaries (RMS, average, nodes) |
| Wavenumber Chart    | Shows momentum distribution                          | Needs description of momentum statistics             |
| Interactive Markers | Draggable measurement tools                          | Keyboard drag functionality needed                   |

#### Screens

| Screen     | Complexity | Key Features Needing Accessibility               |
| ---------- | ---------- | ------------------------------------------------ |
| Intro      | Simplified | Fixed visualization modes, area measurement tool |
| One Well   | Standard   | Full control panel, all visualization options    |
| Two Wells  | Advanced   | Barrier height controls, tunneling descriptions  |
| Many Wells | Complex    | Multiple well configuration, band structure      |

### Gap Analysis

#### Critical Gaps (P0 - Must Have)

1. No `supportsInteractiveDescription` flag in package.json
2. No PDOM structure (screen summary, play area, control area)
3. No keyboard navigation for any controls
4. No semantic HTML elements for interactive components
5. No screen reader announcements for state changes

#### Important Gaps (P1 - Should Have)

1. No accessible names/labels for UI elements
2. No ARIA attributes for complex widgets
3. No focus indicators for keyboard users
4. No alternative descriptions for charts
5. No live region announcements

#### Enhancement Gaps (P2 - Nice to Have)

1. No keyboard shortcuts
2. No keyboard help dialog
3. No high contrast mode considerations
4. No audio descriptions for concepts
5. No accessible documentation

---

## Accessibility Architecture Overview

### PhET PDOM Pattern

PhET simulations use a **Parallel DOM** architecture that maintains a semantic HTML structure alongside the visual canvas rendering. This dual representation ensures:

1. **Visual Users:** See rich canvas graphics with smooth animations
2. **Screen Reader Users:** Navigate semantic HTML with proper ARIA attributes
3. **Keyboard Users:** Access all functionality without mouse
4. **Both Structures Stay Synchronized:** Model changes update both representations

### Three-Section Layout

All PhET screens follow this standard structure:

```
┌─────────────────────────────────────┐
│ Screen Summary                      │
│ (Dynamic overview of current state) │
├─────────────────────────────────────┤
│ Play Area                           │
│ (Interactive simulation content)    │
├─────────────────────────────────────┤
│ Control Area                        │
│ (Settings and parameter controls)   │
└─────────────────────────────────────┘
```

#### 1. Screen Summary

- **Purpose:** Provides accessible overview of simulation state
- **Content:** Dynamic descriptions that update with model changes
- **Updates:** Automatically reflects current potential, energy level, parameters
- **For QPPW:** Summarizes potential type, selected energy level, particle properties, key statistics

#### 2. Play Area

- **Purpose:** Houses interactive visualization elements
- **Content:** Charts, graphs, and visual representations
- **For QPPW:** Energy chart, wavefunction chart, wavenumber chart, interactive tools
- **Accessibility:** Textual descriptions of visual content, interactive element access

#### 3. Control Area

- **Purpose:** Contains all user controls and settings
- **Content:** Sliders, buttons, checkboxes, dropdowns
- **For QPPW:** Potential selector, parameter sliders, display options, simulation controls
- **Accessibility:** Fully keyboard-navigable with proper ARIA attributes

### PDOM Implementation Strategy

```typescript
// Conceptual structure for each screen
class AccessibleScreenView extends Node {
  constructor(model: BaseModel) {
    super();

    // Define PDOM order explicitly
    this.pdomOrder = [
      this.createScreenSummaryNode(), // Section 1
      this.createPlayAreaNode(), // Section 2
      this.createControlAreaNode(), // Section 3
    ];
  }
}
```

### Key Technologies

#### AccessibleNode (Scenery)

Core accessibility properties available on all Nodes:

- `tagName` - HTML element type
- `accessibleName` - Primary label
- `accessibleHelpText` - Supplementary guidance
- `ariaLabel`, `ariaLabelledby` - ARIA labeling
- `ariaRole` - Semantic role override
- `pdomOrder` - Navigation sequence control
- `focusable` - Focus capability

#### UtteranceQueue (Scenery)

Manages live announcements to screen readers:

- Non-visual feedback for state changes
- Priority-based queuing
- Duplicate suppression
- Timing control

#### KeyboardDragListener (Scenery)

Enables keyboard interaction with draggable elements:

- Arrow key movement
- Shift modifier for fine control
- Page Up/Down for large steps
- Customizable step sizes

---

## Implementation Phases

### Phase 0: Foundation Setup

**Estimated Effort:** 1 hour
**Priority:** P0 (Blocking)
**Dependencies:** None

#### Tasks

1. **Enable Interactive Description Support**

   **File:** `package.json`

   Add configuration:

   ```json
   {
     "phet": {
       "simFeatures": {
         "supportsInteractiveDescription": true
       }
     }
   }
   ```

2. **Verify Package Installation**

   Ensure all required Scenery accessibility features are available in SceneryStack 3.0.0

3. **Create Accessibility Directory Structure**

   ```
   src/common/view/accessibility/
   ├── QPPWDescriber.ts
   ├── QPPWAlerter.ts
   ├── ScreenSummaryNode.ts
   └── QPPWKeyboardHelpContent.ts
   ```

4. **Development Environment Setup**
   - Enable keyboard navigation in OS (especially macOS)
   - Install screen reader for testing (NVDA recommended for development)
   - Bookmark a11y-view testing URL

**Success Criteria:**

- ✅ Package.json contains interactive description flag
- ✅ Accessibility directory created
- ✅ Dev environment configured for a11y testing

---

### Phase 1: Screen Structure Implementation

**Estimated Effort:** 8-12 hours
**Priority:** P0 (Blocking)
**Dependencies:** Phase 0

#### Tasks

1. **Create Base Screen Summary Components**

   **File:** `src/common/view/accessibility/ScreenSummaryNode.ts`

   ```typescript
   import { Node } from "scenery";
   import { DerivedProperty } from "axon";
   import type BaseModel from "../../model/BaseModel.js";

   export class ScreenSummaryNode extends Node {
     constructor(model: BaseModel, screenName: string) {
       super({
         tagName: "div",

         // pdom
         labelTagName: "h2",
         labelContent: "Screen Summary",
         descriptionTagName: "p",
         descriptionContent: `${screenName} screen for exploring quantum bound states.`,
       });

       // Current state description
       this.addChild(
         new Node({
           tagName: "p",
           innerContent: new DerivedProperty(
             [model.potentialTypeProperty, model.selectedEnergyLevelProperty],
             (potentialType, energyLevel) => {
               const energy = model.energyLevels[energyLevel];
               return (
                 `Currently exploring a ${potentialType.name} potential well. ` +
                 `Selected energy level ${energyLevel + 1} with energy ${energy.toFixed(3)} eV.`
               );
             },
           ),
         }),
       );

       // Particle properties
       this.addChild(
         new Node({
           tagName: "p",
           innerContent: new DerivedProperty(
             [model.massProperty, model.widthProperty],
             (mass, width) =>
               `Particle mass: ${mass.toFixed(2)} electron masses. Well width: ${width.toFixed(2)} nm.`,
           ),
         }),
       );

       // Statistics summary
       this.addChild(
         new Node({
           tagName: "p",
           innerContent: new DerivedProperty(
             [model.rmsPositionProperty, model.averagePositionProperty],
             (rmsPosition, avgPosition) =>
               `Position statistics: Average ${avgPosition.toFixed(2)} nm, ` +
               `uncertainty (RMS) ${rmsPosition.toFixed(2)} nm.`,
           ),
         }),
       );
     }
   }
   ```

2. **Implement Play Area Node Structure**

   **Modifications:** Each screen's ScreenView file (e.g., `OneWellScreenView.ts`)

   ```typescript
   private createPlayAreaNode(): Node {
     const playAreaNode = new Node({
       tagName: 'div',

       // pdom
       labelTagName: 'h2',
       labelContent: 'Play Area',
       descriptionContent: 'Interactive visualization showing energy levels, ' +
                          'wavefunctions, and probability distributions.'
     });

     // Add charts with PDOM support
     playAreaNode.addChild(this.energyChartNode);
     playAreaNode.addChild(this.waveFunctionChartNode);
     if (this.wavenumberChartNode) {
       playAreaNode.addChild(this.wavenumberChartNode);
     }

     return playAreaNode;
   }
   ```

3. **Implement Control Area Node Structure**

   ```typescript
   private createControlAreaNode(): Node {
     const controlAreaNode = new Node({
       tagName: 'div',

       // pdom
       labelTagName: 'h2',
       labelContent: 'Control Area',
       descriptionContent: 'Controls for adjusting potential parameters, ' +
                          'particle properties, and visualization options.'
     });

     // Add control panel
     controlAreaNode.addChild(this.controlPanelNode);

     // Add simulation controls (play/pause, reset, speed)
     controlAreaNode.addChild(this.simulationControlBar);

     return controlAreaNode;
   }
   ```

4. **Update Screen View PDOM Order**

   ```typescript
   export class OneWellScreenView extends BaseScreenView {
     constructor(model: OneWellModel) {
       super(model);

       // Create accessibility structure
       this.screenSummaryNode = new ScreenSummaryNode(model, "One Well");
       this.playAreaNode = this.createPlayAreaNode();
       this.controlAreaNode = this.createControlAreaNode();

       // Set PDOM order
       this.pdomOrder = [
         this.screenSummaryNode,
         this.playAreaNode,
         this.controlAreaNode,
       ];
     }
   }
   ```

5. **Apply to All Four Screens**
   - IntroScreenView
   - OneWellScreenView
   - TwoWellsScreenView
   - ManyWellsScreenView

**Success Criteria:**

- ✅ All screens have three-section PDOM structure
- ✅ Screen summaries update with model state
- ✅ Keyboard tab navigation follows logical order
- ✅ A11y-view shows proper HTML hierarchy

---

### Phase 2: Interactive Controls - Basic Elements

**Estimated Effort:** 12-16 hours
**Priority:** P1 (Critical)
**Dependencies:** Phase 1

#### 2.1 Play/Pause Button

**File:** `src/common/view/SimulationControlBar.ts`

```typescript
private createPlayPauseButton(): Node {
  const playPauseButton = new PlayPauseButton(this.model.isPlayingProperty, {
    // visual options...

    // pdom
    tagName: 'button',
    innerContent: new DerivedProperty(
      [this.model.isPlayingProperty],
      isPlaying => isPlaying ? 'Pause' : 'Play'
    ),
    accessibleName: 'Play/Pause Simulation',
    helpText: 'Start or stop time evolution of the wavefunction. ' +
              'Keyboard shortcut: Space bar.'
  });

  // Add announcement on toggle
  this.model.isPlayingProperty.link((isPlaying, wasPlaying) => {
    if (wasPlaying !== null) {
      utteranceQueue.addToBack(
        isPlaying ? 'Simulation playing' : 'Simulation paused'
      );
    }
  });

  return playPauseButton;
}
```

#### 2.2 Reset Button

```typescript
private createResetButton(): Node {
  const resetButton = new ResetButton({
    listener: () => {
      this.model.reset();
      utteranceQueue.addToBack('Simulation reset to initial state');
    },

    // pdom
    tagName: 'button',
    innerContent: 'Reset',
    accessibleName: 'Reset Simulation',
    helpText: 'Return all parameters to their initial values. Keyboard shortcut: R.'
  });

  return resetButton;
}
```

#### 2.3 Checkboxes

**File:** `src/common/view/ControlPanelNode.ts`

```typescript
private createShowProbabilityDensityCheckbox(): Node {
  return new Checkbox(
    this.model.showProbabilityDensityProperty,
    new Text('Probability Density'),
    {
      // visual options...

      // pdom
      tagName: 'input',
      inputType: 'checkbox',
      labelContent: 'Show Probability Density',
      helpText: 'Toggle visibility of probability density plot |ψ(x)|². ' +
                'Shows where the particle is most likely to be found.',
      accessibleChecked: this.model.showProbabilityDensityProperty
    }
  );
}

private createShowPhaseColorCheckbox(): Node {
  return new Checkbox(
    this.model.showPhaseColorProperty,
    new Text('Phase Color'),
    {
      // pdom
      tagName: 'input',
      inputType: 'checkbox',
      labelContent: 'Show Phase Color Visualization',
      helpText: 'Color-code the wavefunction by phase angle. ' +
                'Helps visualize time evolution of quantum state.',
      accessibleChecked: this.model.showPhaseColorProperty
    }
  );
}

private createShowZerosCheckbox(): Node {
  return new Checkbox(
    this.model.showZerosProperty,
    new Text('Show Zeros'),
    {
      // pdom
      tagName: 'input',
      inputType: 'checkbox',
      labelContent: 'Show Wavefunction Zeros',
      helpText: 'Mark positions where wavefunction crosses zero. ' +
                'Number of zeros equals energy level minus one.',
      accessibleChecked: this.model.showZerosProperty
    }
  );
}
```

#### 2.4 Radio Button Groups

**Speed Control:**

```typescript
private createSpeedRadioButtons(): Node {
  const speedOptions = [
    { value: 1.0, label: 'Normal Speed' },
    { value: 0.25, label: 'Slow Speed' }
  ];

  return new RadioButtonGroup(
    this.model.timeScaleProperty,
    speedOptions.map(opt => ({
      value: opt.value,
      createNode: () => new Text(opt.label),

      // pdom
      labelContent: opt.label,
      descriptionContent: `Set animation speed to ${opt.label.toLowerCase()}`
    })),
    {
      // visual options...

      // pdom
      labelContent: 'Animation Speed',
      helpText: 'Control time evolution speed. Use arrow keys to navigate options, ' +
                'Space or Enter to select.'
    }
  );
}
```

**Success Criteria:**

- ✅ All buttons keyboard accessible
- ✅ Checkboxes toggle with Space key
- ✅ Radio buttons navigate with arrow keys
- ✅ State changes announced to screen readers
- ✅ Clear focus indicators visible

---

### Phase 3: Interactive Controls - Sliders

**Estimated Effort:** 10-14 hours
**Priority:** P1 (Critical)
**Dependencies:** Phase 2

Sliders are the most common control in QPPW and require careful accessibility implementation.

#### 3.1 Mass Slider

**File:** `src/common/view/ControlPanelNode.ts`

```typescript
import { AccessibleSlider } from 'scenery';

private createMassSlider(): Node {
  const slider = new HSlider(this.model.massProperty, this.massRange, {
    // visual options
    trackSize: new Dimension2(150, 5),
    thumbSize: new Dimension2(15, 30),

    // pdom
    tagName: 'input',
    inputType: 'range',
    labelContent: 'Particle Mass',
    labelTagName: 'label',

    // Accessible value text (updates as user drags)
    accessibleValueText: new DerivedProperty(
      [this.model.massProperty],
      mass => `${mass.toFixed(2)} electron masses`
    ),

    // Keyboard interaction configuration
    keyboardStep: 0.05,           // Arrow keys
    shiftKeyboardStep: 0.01,      // Shift+Arrow (fine control)
    pageKeyboardStep: 0.2,        // Page Up/Down (coarse control)

    // ARIA attributes
    ariaValueMin: this.massRange.min,
    ariaValueMax: this.massRange.max,
    ariaValueNow: this.model.massProperty,
    ariaValueText: new DerivedProperty(
      [this.model.massProperty],
      mass => `${mass.toFixed(2)} electron masses`
    ),

    // Help text
    helpText: 'Adjust particle mass between 0.5 and 1.1 electron masses. ' +
              'Use Left/Right arrow keys for small changes, ' +
              'Shift+Arrow for fine control, ' +
              'Page Up/Down for large steps, ' +
              'Home for minimum, End for maximum.'
  });

  // Alert on significant changes (debounced)
  let alertTimeout: number | null = null;
  this.model.massProperty.link((mass, oldMass) => {
    if (oldMass !== null && Math.abs(mass - oldMass) > 0.001) {
      if (alertTimeout) clearTimeout(alertTimeout);
      alertTimeout = setTimeout(() => {
        utteranceQueue.addToBack(
          `Mass changed to ${mass.toFixed(2)} electron masses. ` +
          `This affects energy level spacing.`
        );
      }, 500); // Debounce: only alert after user stops dragging
    }
  });

  return slider;
}
```

#### 3.2 Width Slider

```typescript
private createWidthSlider(): Node {
  const slider = new HSlider(this.model.widthProperty, this.widthRange, {
    // visual options...

    // pdom
    tagName: 'input',
    inputType: 'range',
    labelContent: 'Well Width',

    accessibleValueText: new DerivedProperty(
      [this.model.widthProperty],
      width => `${width.toFixed(2)} nanometers`
    ),

    keyboardStep: 0.1,
    shiftKeyboardStep: 0.01,
    pageKeyboardStep: 0.5,

    helpText: 'Adjust potential well width. Wider wells have more closely ' +
              'spaced energy levels. Use arrow keys to adjust.'
  });

  // Alert on change
  this.model.widthProperty.lazyLink(width => {
    const numLevels = this.model.energyLevels.length;
    utteranceQueue.addToBack(
      `Well width ${width.toFixed(2)} nm. Found ${numLevels} bound states.`
    );
  });

  return slider;
}
```

#### 3.3 Depth Slider

```typescript
private createDepthSlider(): Node {
  return new HSlider(this.model.depthProperty, this.depthRange, {
    // visual options...

    // pdom
    tagName: 'input',
    inputType: 'range',
    labelContent: 'Well Depth',

    accessibleValueText: new DerivedProperty(
      [this.model.depthProperty],
      depth => `${depth.toFixed(2)} electron volts`
    ),

    keyboardStep: 0.5,
    shiftKeyboardStep: 0.1,
    pageKeyboardStep: 2.0,

    helpText: 'Adjust potential well depth. Deeper wells support more bound states. ' +
              'Use arrow keys to adjust.'
  });
}
```

#### 3.4 Barrier Height Slider (Two Wells Screen)

```typescript
private createBarrierHeightSlider(): Node {
  return new HSlider(this.model.barrierHeightProperty, this.barrierRange, {
    // pdom
    tagName: 'input',
    inputType: 'range',
    labelContent: 'Barrier Height',

    accessibleValueText: new DerivedProperty(
      [this.model.barrierHeightProperty],
      height => `${height.toFixed(2)} electron volts`
    ),

    keyboardStep: 0.1,
    shiftKeyboardStep: 0.01,
    pageKeyboardStep: 0.5,

    helpText: 'Adjust barrier height between wells. Higher barriers reduce tunneling. ' +
              'Affects energy level splitting.'
  });
}
```

#### 3.5 Well Separation Slider (Two Wells & Many Wells)

```typescript
private createSeparationSlider(): Node {
  return new HSlider(this.model.separationProperty, this.separationRange, {
    // pdom
    tagName: 'input',
    inputType: 'range',
    labelContent: 'Well Separation',

    accessibleValueText: new DerivedProperty(
      [this.model.separationProperty],
      separation => `${separation.toFixed(2)} nanometers`
    ),

    keyboardStep: 0.1,
    shiftKeyboardStep: 0.01,
    pageKeyboardStep: 0.5,

    helpText: 'Distance between potential wells. Affects wavefunction overlap and tunneling.'
  });
}
```

#### Implementation Pattern for All Sliders

**Standard accessibility configuration:**

1. **Input type:** `<input type="range">`
2. **Label:** Clear descriptive text
3. **Value text:** Human-readable current value with units
4. **Keyboard steps:**
   - Arrow keys: Medium steps
   - Shift+Arrow: Fine control
   - Page Up/Down: Large steps
   - Home/End: Min/max
5. **ARIA attributes:** Min, max, now, valuetext
6. **Help text:** Explains purpose and keyboard controls
7. **Announcements:** Debounced alerts on significant changes

**Success Criteria:**

- ✅ All sliders keyboard controllable
- ✅ Current values announced to screen readers
- ✅ Multiple step sizes available (fine/medium/coarse)
- ✅ Value changes trigger appropriate alerts
- ✅ Help text explains physics and controls

---

### Phase 4: Interactive Controls - Complex Components

**Estimated Effort:** 16-20 hours
**Priority:** P1-P2 (Critical to Important)
**Dependencies:** Phase 3

#### 4.1 Energy Level Selection

Currently implemented as clickable visual elements in `EnergyChartNode`. Needs keyboard-accessible alternative.

**File:** `src/common/view/EnergyChartNode.ts`

**Approach A: Button List (Recommended)**

```typescript
export class EnergyChartNode extends Node {
  private createEnergyLevelSelector(): Node {
    const selectorNode = new Node({
      tagName: "div",
      ariaRole: "radiogroup",

      // pdom
      labelContent: "Energy Level Selection",
      descriptionContent: new DerivedProperty(
        [this.model.energyLevelsProperty],
        (levels) =>
          `${levels.length} available bound states. Select to view wavefunction.`,
      ),
    });

    // Create button for each energy level
    this.model.energyLevels.forEach((energy, index) => {
      const button = this.createEnergyLevelButton(index, energy);
      selectorNode.addChild(button);
    });

    return selectorNode;
  }

  private createEnergyLevelButton(levelIndex: number, energy: number): Node {
    const button = new Node({
      tagName: "button",
      ariaRole: "radio",

      // pdom
      innerContent: `Level ${levelIndex + 1}`,
      accessibleName: `Energy Level ${levelIndex + 1}`,
      descriptionContent: `Energy: ${energy.toFixed(3)} eV. ${levelIndex} nodes.`,

      // Selected state
      ariaPressed: new DerivedProperty(
        [this.model.selectedEnergyLevelProperty],
        (selected) => selected === levelIndex,
      ),

      // Visual styling based on selection
      cursor: "pointer",
    });

    // Click handler
    button.addInputListener({
      click: () => {
        this.selectEnergyLevel(levelIndex, energy);
      },

      // Keyboard activation (Enter/Space)
      keydown: (event: KeyboardEvent) => {
        if (event.key === "Enter" || event.key === " ") {
          this.selectEnergyLevel(levelIndex, energy);
          event.preventDefault();
        }
      },
    });

    return button;
  }

  private selectEnergyLevel(level: number, energy: number): void {
    this.model.selectedEnergyLevelProperty.value = level;

    // Announce selection
    utteranceQueue.addToBack(
      `Selected energy level ${level + 1}. ` +
        `Energy: ${energy.toFixed(3)} eV. ` +
        `Wavefunction has ${level} node${level !== 1 ? "s" : ""}.`,
      { priority: Utterance.MEDIUM_PRIORITY },
    );
  }
}
```

**Approach B: Listbox (Alternative)**

```typescript
private createEnergyLevelListbox(): Node {
  const listbox = new Node({
    tagName: 'div',
    ariaRole: 'listbox',

    // pdom
    labelContent: 'Energy Levels',
    ariaActiveDescendant: new DerivedProperty(
      [this.model.selectedEnergyLevelProperty],
      level => `level-${level}`
    )
  });

  // Arrow key navigation
  listbox.addInputListener({
    keydown: (event: KeyboardEvent) => {
      const current = this.model.selectedEnergyLevelProperty.value;
      const max = this.model.energyLevels.length - 1;

      if (event.key === 'ArrowUp' && current > 0) {
        this.model.selectedEnergyLevelProperty.value = current - 1;
        event.preventDefault();
      } else if (event.key === 'ArrowDown' && current < max) {
        this.model.selectedEnergyLevelProperty.value = current + 1;
        event.preventDefault();
      } else if (event.key === 'Home') {
        this.model.selectedEnergyLevelProperty.value = 0;
        event.preventDefault();
      } else if (event.key === 'End') {
        this.model.selectedEnergyLevelProperty.value = max;
        event.preventDefault();
      }
    }
  });

  return listbox;
}
```

#### 4.2 Potential Type ComboBox

Most complex interactive component - dropdown with 12+ options.

**File:** `src/common/view/ControlPanelNode.ts`

```typescript
private createPotentialTypeComboBox(listParent: Node): Node {
  // Prepare accessible items
  const items = this.potentialTypes.map(potentialType => ({
    value: potentialType,
    createNode: () => new Text(potentialType.name),

    // pdom - accessible name and description for each option
    a11yLabel: potentialType.name,
    a11yDescription: this.getPotentialDescription(potentialType)
  }));

  const comboBox = new ComboBox(
    this.model.potentialTypeProperty,
    items,
    listParent,
    {
      // visual options...
      buttonTouchAreaXDilation: 10,
      buttonTouchAreaYDilation: 5,

      // pdom
      tagName: 'div', // Container
      ariaRole: 'combobox',

      accessibleName: 'Potential Type',
      helpText: 'Select quantum potential well type. ' +
                'Press Enter to open menu, use arrow keys to navigate, ' +
                'Enter to select, Escape to close.',

      // Current selection announcement
      accessibleNameBehavior: (node, options, accessibleName) => {
        const current = this.model.potentialTypeProperty.value;
        return `${accessibleName}: ${current.name}`;
      }
    }
  );

  // Announce when potential type changes
  this.model.potentialTypeProperty.link((newType, oldType) => {
    if (oldType !== null) {
      const numLevels = this.model.energyLevels.length;
      utteranceQueue.addToBack(
        `Potential changed to ${newType.name}. ` +
        `Found ${numLevels} bound state${numLevels !== 1 ? 's' : ''}. ` +
        `${this.getPotentialDescription(newType)}`,
        { priority: Utterance.HIGH_PRIORITY }
      );
    }
  });

  return comboBox;
}

private getPotentialDescription(potentialType: PotentialType): string {
  const descriptions: Record<string, string> = {
    'Infinite Square Well':
      'Particle confined in rigid box with infinite barriers. ' +
      'Exactly solvable with uniform energy spacing.',

    'Finite Square Well':
      'Square well with finite barrier height. ' +
      'Realistic model with exponentially decaying wavefunctions outside well.',

    'Harmonic Oscillator':
      'Quadratic potential resembling mass on spring. ' +
      'Energy levels uniformly spaced.',

    'Morse Potential':
      'Models molecular vibrations with anharmonic oscillations. ' +
      'Energy spacing decreases at higher levels.',

    'Pöschl-Teller Potential':
      'Exactly solvable hyperbolic secant potential. ' +
      'Important in quantum scattering theory.',

    'Rosen-Morse Potential':
      'Variant of Pöschl-Teller with different asymptotic behavior.',

    'Eckart Potential':
      'Barrier potential used in molecular physics. ' +
      'Models quantum tunneling through barriers.',

    'Asymmetric Triangle':
      'Tilted potential well with linear slope. ' +
      'Models particle in electric field.',

    'Triangular Potential':
      'V-shaped potential well. ' +
      'Related to Airy functions.',

    'Coulomb 1D':
      'One-dimensional hydrogen-like attractive potential. ' +
      'Energy levels follow 1/n² pattern.',

    'Coulomb 3D':
      'Three-dimensional hydrogen atom potential. ' +
      'Includes angular momentum quantum numbers.',

    'Double Square Well':
      'Two square wells separated by barrier. ' +
      'Demonstrates quantum tunneling and energy level splitting.'
  };

  return descriptions[potentialType.name] || 'Custom quantum potential.';
}
```

**Success Criteria:**

- ✅ Energy levels selectable via keyboard
- ✅ Level selection announced with energy and node count
- ✅ Potential type ComboBox keyboard navigable
- ✅ Each potential type has descriptive text
- ✅ Changes trigger informative announcements

---

### Phase 5: Chart Descriptions & Visualizations

**Estimated Effort:** 14-18 hours
**Priority:** P2 (Important)
**Dependencies:** Phase 4

#### 5.1 Energy Chart Accessible Description

**File:** `src/common/view/EnergyChartNode.ts`

```typescript
export class EnergyChartNode extends Node {
  constructor(model: BaseModel, ...) {
    super({
      tagName: 'div',
      ariaRole: 'region',

      // pdom
      labelContent: 'Energy Level Diagram',

      // Dynamic description
      descriptionContent: new DerivedProperty(
        [
          model.potentialTypeProperty,
          model.energyLevelsProperty,
          model.selectedEnergyLevelProperty,
          model.widthProperty,
          model.depthProperty
        ],
        (potentialType, energyLevels, selected, width, depth) => {
          return this.createEnergyChartDescription(
            potentialType, energyLevels, selected, width, depth
          );
        }
      )
    });

    // ... rest of implementation
  }

  private createEnergyChartDescription(
    potentialType: PotentialType,
    energyLevels: number[],
    selected: number,
    width: number,
    depth: number
  ): string {
    const numLevels = energyLevels.length;
    const selectedEnergy = energyLevels[selected];
    const groundEnergy = energyLevels[0];

    let description = `${potentialType.name} potential well. `;
    description += `Width: ${width.toFixed(2)} nm. `;

    if (depth !== undefined) {
      description += `Depth: ${depth.toFixed(2)} eV. `;
    }

    description += `\n\n`;
    description += `Found ${numLevels} bound state${numLevels !== 1 ? 's' : ''}. `;
    description += `Ground state energy: ${groundEnergy.toFixed(3)} eV. `;

    if (numLevels > 1) {
      const firstExcited = energyLevels[1];
      const spacing = firstExcited - groundEnergy;
      description += `First excited state: ${firstExcited.toFixed(3)} eV. `;
      description += `Energy spacing: ${spacing.toFixed(3)} eV. `;
    }

    description += `\n\n`;
    description += `Currently viewing level ${selected + 1} ` +
                   `with energy ${selectedEnergy.toFixed(3)} eV.`;

    // Classical turning points
    const turningPoints = this.calculateTurningPoints(selected);
    if (turningPoints) {
      description += ` Classical turning points at ` +
                     `${turningPoints.left.toFixed(2)} nm and ` +
                     `${turningPoints.right.toFixed(2)} nm.`;
    }

    return description;
  }
}
```

#### 5.2 Wavefunction Chart Description

**File:** `src/common/view/WaveFunctionChartNode.ts`

```typescript
export class WaveFunctionChartNode extends Node {
  constructor(model: BaseModel, ...) {
    super({
      tagName: 'div',
      ariaRole: 'region',

      // pdom
      labelContent: 'Wavefunction Visualization',
      descriptionContent: new DerivedProperty(
        [
          model.displayModeProperty,
          model.selectedEnergyLevelProperty,
          model.rmsPositionProperty,
          model.averagePositionProperty,
          model.wavefunctionProperty
        ],
        (...args) => this.createWavefunctionDescription(...args)
      )
    });
  }

  private createWavefunctionDescription(
    displayMode: string,
    level: number,
    rmsPosition: number,
    avgPosition: number,
    wavefunction: ComplexArray
  ): string {
    let desc = `Wavefunction for energy level ${level + 1}. `;

    // Display mode
    if (displayMode === 'probability-density') {
      desc += `Showing probability density |ψ(x)|². `;
      desc += `This indicates where the particle is most likely to be found. `;
    } else if (displayMode === 'magnitude') {
      desc += `Showing wavefunction magnitude |ψ(x)|. `;
    } else if (displayMode === 'real-imaginary') {
      desc += `Showing real and imaginary components of ψ(x). `;
      desc += `Real part in blue, imaginary part in red. `;
    } else if (displayMode === 'phase') {
      desc += `Showing phase angle of complex wavefunction. `;
    }

    desc += `\n\n`;

    // Statistical properties
    desc += `Position statistics: `;
    desc += `Average position: ${avgPosition.toFixed(2)} nanometers. `;
    desc += `Position uncertainty (RMS, Δx): ${rmsPosition.toFixed(2)} nm. `;

    // Heisenberg uncertainty relation
    const momentumUncertainty = this.calculateMomentumUncertainty();
    const uncertaintyProduct = rmsPosition * momentumUncertainty;
    const hbar = 1.054571817e-34; // Reduced Planck constant
    desc += `Heisenberg uncertainty product Δx·Δp: ${uncertaintyProduct.toFixed(2)}. `;
    desc += `Minimum allowed: ${(hbar/2).toExponential(2)}. `;

    desc += `\n\n`;

    // Nodes (zeros)
    const zeros = this.countWavefunctionZeros(wavefunction);
    desc += `Wavefunction has ${zeros} node${zeros !== 1 ? 's' : ''} ` +
            `(zero crossing${zeros !== 1 ? 's' : ''}). `;
    desc += `For quantum number n=${level + 1}, expected ${level} nodes. `;

    // Symmetry
    const symmetry = this.analyzeSymmetry(wavefunction);
    if (symmetry === 'even') {
      desc += `Wavefunction is even (symmetric). `;
    } else if (symmetry === 'odd') {
      desc += `Wavefunction is odd (antisymmetric). `;
    }

    return desc;
  }

  private countWavefunctionZeros(wavefunction: ComplexArray): number {
    let count = 0;
    for (let i = 1; i < wavefunction.length; i++) {
      const prev = wavefunction[i-1].real;
      const curr = wavefunction[i].real;
      if (prev * curr < 0) count++;
    }
    return count;
  }

  private analyzeSymmetry(wavefunction: ComplexArray): string {
    const mid = Math.floor(wavefunction.length / 2);
    let evenScore = 0;
    let oddScore = 0;

    for (let i = 0; i < mid; i++) {
      const left = wavefunction[i].real;
      const right = wavefunction[wavefunction.length - 1 - i].real;

      if (Math.abs(left - right) < 0.01) evenScore++;
      if (Math.abs(left + right) < 0.01) oddScore++;
    }

    if (evenScore > oddScore * 1.5) return 'even';
    if (oddScore > evenScore * 1.5) return 'odd';
    return 'neither';
  }
}
```

#### 5.3 Wavenumber Chart Description (Intro Screen)

**File:** `src/common/view/WavenumberChartNode.ts`

```typescript
export class WavenumberChartNode extends Node {
  constructor(model: BaseModel) {
    super({
      tagName: "div",
      ariaRole: "region",

      // pdom
      labelContent: "Momentum Distribution",
      descriptionContent: new DerivedProperty(
        [
          model.rmsWavenumberProperty,
          model.averageWavenumberProperty,
          model.momentumDistributionProperty,
        ],
        (rmsK, avgK, distribution) => {
          let desc = `Momentum space representation showing |φ(k)|². `;
          desc += `This is the Fourier transform of the position wavefunction. `;
          desc += `\n\n`;

          desc += `Average wavenumber: ${avgK.toFixed(3)} inverse nanometers. `;
          desc += `Wavenumber uncertainty (Δk): ${rmsK.toFixed(3)} nm⁻¹. `;

          const avgMomentum = avgK * hbar;
          desc += `\n\n`;
          desc += `Average momentum: ${avgMomentum.toExponential(2)} kg·m/s. `;

          // Uncertainty relation in momentum space
          const positionUncertainty = model.rmsPositionProperty.value;
          const uncertaintyProduct = positionUncertainty * rmsK;
          desc += `Position-momentum uncertainty: Δx·Δk = ${uncertaintyProduct.toFixed(2)}. `;
          desc += `Minimum allowed: 0.5 (dimensionless). `;

          return desc;
        },
      ),
    });
  }
}
```

**Success Criteria:**

- ✅ All charts have meaningful text descriptions
- ✅ Descriptions include key physics concepts
- ✅ Statistical properties clearly stated
- ✅ Descriptions update dynamically with model state
- ✅ Screen reader users can understand visualizations without seeing them

---

### Phase 6: Interactive Tools Accessibility

**Estimated Effort:** 12-16 hours
**Priority:** P3 (Enhancement)
**Dependencies:** Phase 5

#### 6.1 Area Measurement Tool

**File:** `src/common/view/chart-tools/AreaMeasurementTool.ts`

```typescript
export class AreaMeasurementTool extends Node {
  private leftMarker: DraggableMarkerNode;
  private rightMarker: DraggableMarkerNode;
  private integratedProbabilityProperty: Property<number>;

  constructor(model: BaseModel, chartTransform: ChartTransform) {
    super({
      tagName: "div",

      // pdom
      labelContent: "Area Measurement Tool",
      descriptionContent:
        "Drag markers to measure probability between two positions. " +
        "Use keyboard to fine-tune marker positions.",
    });

    // Create draggable markers
    this.leftMarker = this.createDraggableMarker(
      "Left",
      model.leftMarkerProperty,
    );
    this.rightMarker = this.createDraggableMarker(
      "Right",
      model.rightMarkerProperty,
    );

    // Probability display
    this.integratedProbabilityProperty = new DerivedProperty(
      [
        model.leftMarkerProperty,
        model.rightMarkerProperty,
        model.wavefunctionProperty,
      ],
      (left, right, wavefunction) =>
        this.calculateProbability(left, right, wavefunction),
    );

    // Accessible readout
    this.addChild(
      new Node({
        tagName: "div",
        ariaRole: "status",
        ariaLive: "polite",
        innerContent: new DerivedProperty(
          [
            model.leftMarkerProperty,
            model.rightMarkerProperty,
            this.integratedProbabilityProperty,
          ],
          (left, right, prob) =>
            `Measuring from ${left.toFixed(2)} to ${right.toFixed(2)} nm. ` +
            `Integrated probability: ${(prob * 100).toFixed(1)} percent.`,
        ),
      }),
    );

    this.addChild(this.leftMarker);
    this.addChild(this.rightMarker);
  }

  private createDraggableMarker(
    label: string,
    positionProperty: Property<number>,
  ): Node {
    const marker = new Node({
      tagName: "div",
      ariaRole: "slider",

      // pdom
      accessibleName: `${label} Measurement Marker`,
      labelContent: `${label} boundary for probability integration`,

      // Current position
      ariaValueText: new DerivedProperty(
        [positionProperty],
        (position) => `Position: ${position.toFixed(2)} nanometers`,
      ),

      // Range
      ariaValueMin: -5, // Adjust based on actual range
      ariaValueMax: 5,
      ariaValueNow: positionProperty,

      helpText:
        "Use Left/Right arrow keys to move marker. " +
        "Shift+Arrow for fine control (0.01 nm steps). " +
        "Page Up/Down for large steps (0.5 nm). " +
        "Home/End for range limits.",
    });

    // Keyboard drag listener
    const keyboardDragListener = new KeyboardDragListener({
      positionProperty: positionProperty,
      transform: this.chartTransform,

      dragDelta: 0.1, // nm per arrow key press
      shiftDragDelta: 0.01, // nm per shift+arrow

      // Announce on drag end
      end: () => {
        const position = positionProperty.value;
        const probability = this.integratedProbabilityProperty.value;
        utteranceQueue.addToBack(
          `${label} marker at ${position.toFixed(2)} nanometers. ` +
            `Integrated probability: ${(probability * 100).toFixed(1)} percent.`,
        );
      },
    });

    marker.addInputListener(keyboardDragListener);

    // Visual rendering
    // ... existing visual code

    return marker;
  }

  private calculateProbability(
    left: number,
    right: number,
    wavefunction: ComplexArray,
  ): number {
    // Trapezoidal integration of |ψ(x)|² between markers
    // ... existing calculation code
    return probability;
  }
}
```

#### 6.2 Curvature Tool

**File:** `src/common/view/chart-tools/CurvatureTool.ts`

```typescript
export class CurvatureTool extends Node {
  constructor(model: BaseModel, enabledProperty: Property<boolean>) {
    super({
      tagName: "div",

      // pdom
      labelContent: "Curvature Visualization",
      descriptionContent: new DerivedProperty([enabledProperty], (enabled) =>
        enabled
          ? "Showing second derivative d²ψ/dx². Curvature is proportional to (V(x) - E)ψ(x) " +
            "according to the Schrödinger equation. Positive curvature where V > E, negative where V < E."
          : "Curvature visualization disabled.",
      ),

      // Visibility tied to enabled state
      accessibleVisible: enabledProperty,
    });

    // Toggle checkbox (in control panel)
    // ... handled separately in ControlPanelNode
  }
}
```

#### 6.3 Derivative Tool

Similar pattern to Curvature Tool with first derivative information.

#### 6.4 Zeros Visualization

**File:** `src/common/view/chart-tools/ZerosVisualization.ts`

```typescript
export class ZerosVisualization extends Node {
  constructor(model: BaseModel, enabledProperty: Property<boolean>) {
    super({
      tagName: "div",
      ariaRole: "status",
      ariaLive: "polite",

      // pdom
      labelContent: "Wavefunction Zeros",
      descriptionContent: new DerivedProperty(
        [model.selectedEnergyLevelProperty, model.wavefunctionZerosProperty],
        (level, zeros) => {
          if (!enabledProperty.value) return "Zeros visualization disabled.";

          const positions = zeros.map((z) => z.toFixed(2)).join(", ");
          return (
            `Energy level ${level + 1} has ${zeros.length} node${zeros.length !== 1 ? "s" : ""} ` +
            `(zero crossing${zeros.length !== 1 ? "s" : ""}) ` +
            `at positions: ${positions} nanometers.`
          );
        },
      ),
    });
  }
}
```

**Success Criteria:**

- ✅ Area measurement markers keyboard draggable
- ✅ Marker positions announced clearly
- ✅ Integrated probability value accessible
- ✅ Visualization tools have descriptive text
- ✅ Tool states announced when enabled/disabled

---

### Phase 7: Live Alerts & Dynamic Feedback

**Estimated Effort:** 8-12 hours
**Priority:** P2 (Important)
**Dependencies:** All previous phases

#### 7.1 Create Centralized Alerter

**File:** `src/common/view/accessibility/QPPWAlerter.ts`

```typescript
import { Alerter } from "scenery-phet";
import utteranceQueue from "utteranceQueue.js";
import type BaseModel from "../../model/BaseModel.js";

export class QPPWAlerter extends Alerter {
  private model: BaseModel;

  constructor(model: BaseModel) {
    super();
    this.model = model;

    // Set up listeners for important state changes
    this.setupAlerts();
  }

  private setupAlerts(): void {
    // Energy level changes
    this.model.selectedEnergyLevelProperty.lazyLink((level, oldLevel) => {
      this.alertEnergyLevelChange(level, oldLevel);
    });

    // Potential type changes
    this.model.potentialTypeProperty.lazyLink((newType, oldType) => {
      this.alertPotentialTypeChange(newType, oldType);
    });

    // Play/pause state
    this.model.isPlayingProperty.lazyLink((isPlaying) => {
      this.alertPlaybackStateChange(isPlaying);
    });

    // Significant parameter changes
    this.model.massProperty.lazyLink((mass) => {
      this.alertMassChange(mass);
    });

    this.model.widthProperty.lazyLink((width) => {
      this.alertWidthChange(width);
    });
  }

  private alertEnergyLevelChange(level: number, oldLevel: number): void {
    const energy = this.model.energyLevels[level];
    const nodes = level; // Quantum number n-1

    const alert =
      `Selected energy level ${level + 1}. ` +
      `Energy: ${energy.toFixed(3)} electron volts. ` +
      `Wavefunction has ${nodes} node${nodes !== 1 ? "s" : ""}.`;

    utteranceQueue.addToBack(alert, {
      priority: Utterance.MEDIUM_PRIORITY,
    });
  }

  private alertPotentialTypeChange(
    newType: PotentialType,
    oldType: PotentialType,
  ): void {
    const numLevels = this.model.energyLevels.length;
    const groundEnergy = this.model.energyLevels[0];

    const alert =
      `Potential changed to ${newType.name}. ` +
      `Found ${numLevels} bound state${numLevels !== 1 ? "s" : ""}. ` +
      `Ground state energy: ${groundEnergy.toFixed(3)} eV.`;

    utteranceQueue.addToBack(alert, {
      priority: Utterance.HIGH_PRIORITY,
    });
  }

  private alertPlaybackStateChange(isPlaying: boolean): void {
    const alert = isPlaying
      ? "Simulation playing. Wavefunction evolving in time."
      : "Simulation paused.";

    utteranceQueue.addToBack(alert);
  }

  private alertMassChange(mass: number): void {
    // Debounced - only alert after user stops adjusting
    this.debouncedAlert(
      `Particle mass changed to ${mass.toFixed(2)} electron masses. ` +
        `Energy levels recalculated.`,
      500,
    );
  }

  private alertWidthChange(width: number): void {
    this.debouncedAlert(
      `Well width changed to ${width.toFixed(2)} nanometers. ` +
        `Found ${this.model.energyLevels.length} bound states.`,
      500,
    );
  }

  // Helper for debounced alerts
  private debouncedAlertTimer: number | null = null;
  private debouncedAlert(message: string, delay: number): void {
    if (this.debouncedAlertTimer) {
      clearTimeout(this.debouncedAlertTimer);
    }

    this.debouncedAlertTimer = setTimeout(() => {
      utteranceQueue.addToBack(message);
      this.debouncedAlertTimer = null;
    }, delay);
  }

  // Context-aware alerts
  public alertSuperstateChange(superpositionType: string): void {
    const alerts: Record<string, string> = {
      single:
        "Single eigenstate selected. Stationary state, no time evolution.",
      "two-state":
        "Two-state superposition. Wavefunction oscillates between wells.",
      wavepacket: "Gaussian wavepacket created. Evolves as localized particle.",
      coherent: "Coherent state superposition. Minimal uncertainty wavepacket.",
      custom: "Custom superposition configured.",
    };

    utteranceQueue.addToBack(
      alerts[superpositionType] || "Superposition state changed.",
      { priority: Utterance.MEDIUM_PRIORITY },
    );
  }

  public alertResetAllSimulation(): void {
    utteranceQueue.addToBack(
      "Simulation reset. All parameters returned to initial values.",
      { priority: Utterance.HIGH_PRIORITY },
    );
  }
}
```

#### 7.2 Alert Priorities & Timing

**Guidelines for alert usage:**

| Event Type          | Priority | Timing             | Example                   |
| ------------------- | -------- | ------------------ | ------------------------- |
| Critical errors     | HIGH     | Immediate          | "Calculation failed"      |
| Major state changes | HIGH     | Immediate          | "Potential type changed"  |
| User actions        | MEDIUM   | Immediate          | "Energy level selected"   |
| Parameter changes   | MEDIUM   | Debounced (500ms)  | "Mass changed to..."      |
| Automatic updates   | LOW      | Debounced (1000ms) | "Wavefunction normalized" |
| Informational       | LOW      | On focus/request   | Help text, descriptions   |

**Anti-patterns to avoid:**

- ❌ Alerting every frame during animation
- ❌ Announcing minor numerical changes (e.g., every 0.01 eV change)
- ❌ Duplicate alerts for same information
- ❌ Interrupting user mid-action with low-priority alerts
- ❌ Over-verbose announcements (keep under 2-3 sentences)

**Success Criteria:**

- ✅ Important state changes announced clearly
- ✅ Alerts don't interrupt user workflow
- ✅ No redundant or excessive announcements
- ✅ Priorities set appropriately
- ✅ Timing optimized for comprehension

---

### Phase 8: Keyboard Shortcuts & Navigation

**Estimated Effort:** 6-10 hours
**Priority:** P3 (Enhancement)
**Dependencies:** Phases 1-6

#### 8.1 Global Keyboard Shortcuts

**File:** `src/common/view/BaseScreenView.ts`

```typescript
import { GlobalKeyboardListener } from "scenery";

export class BaseScreenView extends ScreenView {
  private keyboardListener: GlobalKeyboardListener;

  constructor(model: BaseModel) {
    super();

    // ... other initialization

    this.setupKeyboardShortcuts();
  }

  private setupKeyboardShortcuts(): void {
    this.keyboardListener = new GlobalKeyboardListener({
      // Play/Pause - Space bar
      " ": () => {
        this.model.isPlayingProperty.toggle();
        return true; // Handled
      },

      // Reset - R key
      r: () => {
        this.reset();
        utteranceQueue.addToBack("Simulation reset to initial state");
        return true;
      },

      // Energy level navigation (when chart focused)
      ArrowUp: () => {
        if (this.energyChartHasFocus()) {
          this.incrementEnergyLevel();
          return true;
        }
        return false; // Not handled - allow default tab navigation
      },

      ArrowDown: () => {
        if (this.energyChartHasFocus()) {
          this.decrementEnergyLevel();
          return true;
        }
        return false;
      },

      Home: () => {
        if (this.energyChartHasFocus()) {
          this.model.selectedEnergyLevelProperty.value = 0;
          utteranceQueue.addToBack("Ground state selected");
          return true;
        }
        return false;
      },

      End: () => {
        if (this.energyChartHasFocus()) {
          const maxLevel = this.model.energyLevels.length - 1;
          this.model.selectedEnergyLevelProperty.value = maxLevel;
          utteranceQueue.addToBack(
            `Highest bound state selected, level ${maxLevel + 1}`,
          );
          return true;
        }
        return false;
      },

      // Speed toggle - S key
      s: () => {
        const current = this.model.timeScaleProperty.value;
        this.model.timeScaleProperty.value = current === 1.0 ? 0.25 : 1.0;
        utteranceQueue.addToBack(
          `Speed set to ${current === 1.0 ? "slow" : "normal"}`,
        );
        return true;
      },

      // Help - ? or H key
      "?": () => {
        this.showKeyboardHelp();
        return true;
      },
      h: () => {
        this.showKeyboardHelp();
        return true;
      },
    });

    this.addInputListener(this.keyboardListener);
  }

  private energyChartHasFocus(): boolean {
    // Check if energy chart or its children have focus
    return (
      this.energyChartNode.focused ||
      this.energyChartNode.getDescendantsWithFocus().length > 0
    );
  }

  private incrementEnergyLevel(): void {
    const current = this.model.selectedEnergyLevelProperty.value;
    const max = this.model.energyLevels.length - 1;
    if (current < max) {
      this.model.selectedEnergyLevelProperty.value = current + 1;
    }
  }

  private decrementEnergyLevel(): void {
    const current = this.model.selectedEnergyLevelProperty.value;
    if (current > 0) {
      this.model.selectedEnergyLevelProperty.value = current - 1;
    }
  }

  private showKeyboardHelp(): void {
    // Open keyboard help dialog
    const helpDialog = new KeyboardHelpDialog(this.keyboardHelpContent);
    helpDialog.show();
  }
}
```

#### 8.2 Keyboard Help Content

**File:** `src/common/view/accessibility/QPPWKeyboardHelpContent.ts`

```typescript
import { KeyboardHelpSection, KeyboardHelpIconFactory } from "scenery-phet";
import { VBox } from "scenery";

export class QPPWKeyboardHelpContent extends VBox {
  constructor() {
    // General controls section
    const generalControlsSection = new KeyboardHelpSection("General Controls", [
      KeyboardHelpIconFactory.iconRow("Space", "Play or pause simulation"),
      KeyboardHelpIconFactory.iconRow("R", "Reset simulation to initial state"),
      KeyboardHelpIconFactory.iconRow(
        "S",
        "Toggle between normal and slow speed",
      ),
      KeyboardHelpIconFactory.iconRow(
        ["?", "or", "H"],
        "Show this keyboard help",
      ),
      KeyboardHelpIconFactory.tabRow("Move to next control"),
      KeyboardHelpIconFactory.shiftPlusTabRow("Move to previous control"),
    ]);

    // Energy level navigation section
    const energyLevelSection = new KeyboardHelpSection(
      "Energy Level Selection",
      [
        KeyboardHelpIconFactory.iconRow(
          "Up Arrow",
          "Select higher energy level",
        ),
        KeyboardHelpIconFactory.iconRow(
          "Down Arrow",
          "Select lower energy level",
        ),
        KeyboardHelpIconFactory.iconRow(
          "Home",
          "Select ground state (level 1)",
        ),
        KeyboardHelpIconFactory.iconRow("End", "Select highest bound state"),
        KeyboardHelpIconFactory.iconRow(
          ["Enter", "or", "Space"],
          "Activate selected level",
        ),
      ],
    );

    // Slider controls section
    const sliderControlsSection = new KeyboardHelpSection("Slider Controls", [
      KeyboardHelpIconFactory.iconRow(
        ["Left Arrow", "or", "Right Arrow"],
        "Adjust value in small steps",
      ),
      KeyboardHelpIconFactory.iconRow(
        ["Shift", "+", "Arrow"],
        "Adjust value in fine steps",
      ),
      KeyboardHelpIconFactory.iconRow(
        ["Page Up", "or", "Page Down"],
        "Adjust value in large steps",
      ),
      KeyboardHelpIconFactory.iconRow("Home", "Jump to minimum value"),
      KeyboardHelpIconFactory.iconRow("End", "Jump to maximum value"),
    ]);

    // Combo box and dropdowns section
    const comboBoxSection = new KeyboardHelpSection("Dropdown Menus", [
      KeyboardHelpIconFactory.iconRow("Enter", "Open menu"),
      KeyboardHelpIconFactory.iconRow(
        ["Up Arrow", "or", "Down Arrow"],
        "Navigate options",
      ),
      KeyboardHelpIconFactory.iconRow("Enter", "Select option"),
      KeyboardHelpIconFactory.iconRow("Escape", "Close menu without selecting"),
    ]);

    // Measurement tool section
    const measurementToolSection = new KeyboardHelpSection(
      "Area Measurement Tool",
      [
        KeyboardHelpIconFactory.iconRow(
          ["Left Arrow", "or", "Right Arrow"],
          "Move marker (medium steps)",
        ),
        KeyboardHelpIconFactory.iconRow(
          ["Shift", "+", "Arrow"],
          "Move marker (fine steps, 0.01 nm)",
        ),
        KeyboardHelpIconFactory.iconRow(
          ["Page Up", "or", "Page Down"],
          "Move marker (large steps, 0.5 nm)",
        ),
      ],
    );

    super({
      children: [
        generalControlsSection,
        energyLevelSection,
        sliderControlsSection,
        comboBoxSection,
        measurementToolSection,
      ],
      spacing: 15,
      align: "left",

      // pdom
      tagName: "div",
      labelContent: "Keyboard Shortcuts Reference",
      descriptionContent:
        "Complete list of keyboard commands for navigating " +
        "and controlling the quantum particle simulation.",
    });
  }
}
```

**Success Criteria:**

- ✅ All major actions have keyboard shortcuts
- ✅ Shortcuts follow standard conventions
- ✅ Help dialog clearly documents all shortcuts
- ✅ Shortcuts don't conflict with browser/screen reader keys
- ✅ Visual feedback for keyboard actions

---

### Phase 9: Testing & Validation

**Estimated Effort:** 16-20 hours
**Priority:** P0 (Critical)
**Dependencies:** All implementation phases

#### 9.1 Manual Testing Checklist

**Screen Reader Testing:**

```markdown
## NVDA (Windows) Testing Checklist

### Navigation

- [ ] Tab order follows logical sequence (summary → play area → controls)
- [ ] All interactive elements reachable via tab
- [ ] Shift+Tab navigates backwards correctly
- [ ] Focus indicators clearly visible

### Content

- [ ] Screen summary reads correctly on page load
- [ ] Chart descriptions are announced when focused
- [ ] Dynamic content updates announced appropriately
- [ ] No unnecessary or redundant announcements

### Controls

- [ ] Buttons: Name and role announced correctly
- [ ] Checkboxes: State (checked/unchecked) announced
- [ ] Sliders: Current value and range announced
- [ ] ComboBox: Closed/expanded state clear, options navigable
- [ ] Energy levels: Selection announced with energy value

### Interactions

- [ ] Slider value changes announced (debounced)
- [ ] Energy level changes announced clearly
- [ ] Potential type changes announced with context
- [ ] Play/pause state changes announced
- [ ] Reset action announced

### Live Regions

- [ ] Aria-live alerts trigger at appropriate times
- [ ] Alert content is clear and concise
- [ ] No alert spam during continuous changes
- [ ] High-priority alerts interrupt appropriately

## JAWS (Windows) Testing Checklist

[Same structure as NVDA]

## VoiceOver (macOS/iOS) Testing Checklist

[Same structure adapted for VoiceOver commands]

## TalkBack (Android) Testing Checklist

[Same structure adapted for TalkBack gestures]
```

#### 9.2 Keyboard-Only Testing

```markdown
## Keyboard Navigation Testing

### Basic Navigation

- [ ] Can reach all controls using only keyboard
- [ ] Tab order matches visual/logical flow
- [ ] Focus visible at all times
- [ ] No keyboard traps (can escape all widgets)
- [ ] Skip links or landmarks available for efficient navigation

### Control Interaction

- [ ] Space: Play/pause works
- [ ] R: Reset works
- [ ] All sliders: Arrow keys adjust values
- [ ] All sliders: Shift+arrow for fine control
- [ ] All sliders: Page up/down for coarse control
- [ ] All sliders: Home/end for min/max
- [ ] Checkboxes: Space toggles state
- [ ] Radio buttons: Arrow keys navigate, space selects
- [ ] Buttons: Enter or space activates
- [ ] ComboBox: Enter opens, arrows navigate, enter selects

### Energy Levels

- [ ] Can select any energy level via keyboard
- [ ] Up/down arrows work when focused
- [ ] Home selects ground state
- [ ] End selects highest state

### Measurement Tool

- [ ] Markers draggable with arrow keys
- [ ] Fine/medium/coarse steps work correctly
- [ ] Position announced after movement

### Shortcuts

- [ ] All documented shortcuts work as expected
- [ ] Shortcuts don't interfere with screen reader commands
- [ ] Help dialog accessible via keyboard
```

#### 9.3 Automated Testing

**File:** `tests/accessibility-tests.ts`

```typescript
import { describe, it, assert } from "@testing";
import { OneWellModel } from "../src/one-well/model/OneWellModel.js";
import { OneWellScreenView } from "../src/one-well/view/OneWellScreenView.js";

describe("PDOM Structure Tests", () => {
  let model: OneWellModel;
  let screenView: OneWellScreenView;

  beforeEach(() => {
    model = new OneWellModel();
    screenView = new OneWellScreenView(model);
  });

  it("should have three main PDOM sections", () => {
    const pdomOrder = screenView.pdomOrder;
    assert.equal(pdomOrder.length, 3, "Should have 3 main sections");

    const labels = pdomOrder.map((node) => node.labelContent);
    assert.include(labels, "Screen Summary");
    assert.include(labels, "Play Area");
    assert.include(labels, "Control Area");
  });

  it("should have accessible names for all interactive elements", () => {
    const interactiveNodes = screenView.getInteractiveNodes();

    interactiveNodes.forEach((node) => {
      assert.isTrue(
        node.accessibleName !== null ||
          node.innerContent !== null ||
          node.labelContent !== null,
        `Node ${node.constructor.name} missing accessible name`,
      );
    });
  });

  it("should have proper ARIA roles for semantic elements", () => {
    // Check buttons
    const playPauseButton = screenView.findNodeByName("Play/Pause");
    assert.equal(playPauseButton.tagName, "button");

    // Check sliders
    const massSlider = screenView.findNodeByLabel("Particle Mass");
    assert.equal(massSlider.tagName, "input");
    assert.equal(massSlider.inputType, "range");
  });

  it("should announce energy level changes", () => {
    // Clear utterance queue
    utteranceQueue.clear();

    // Change energy level
    model.selectedEnergyLevelProperty.value = 2;

    // Check that announcement was queued
    const alerts = utteranceQueue.queue;
    assert.isTrue(alerts.length > 0, "Should have queued alert");
    assert.include(alerts[0].toString(), "Energy level 3");
    assert.include(alerts[0].toString(), "eV");
  });

  it("should update dynamic descriptions when model changes", () => {
    const screenSummary = screenView.screenSummaryNode;
    const initialDescription = screenSummary.descriptionContent;

    // Change potential type
    model.potentialTypeProperty.value = potentialTypes[1];

    const updatedDescription = screenSummary.descriptionContent;
    assert.notEqual(
      initialDescription,
      updatedDescription,
      "Description should update when potential changes",
    );
  });

  it("should have keyboard step sizes configured for sliders", () => {
    const massSlider = screenView.findNodeByLabel("Particle Mass");

    assert.isNumber(massSlider.keyboardStep);
    assert.isNumber(massSlider.shiftKeyboardStep);
    assert.isNumber(massSlider.pageKeyboardStep);
    assert.isTrue(
      massSlider.shiftKeyboardStep < massSlider.keyboardStep,
      "Shift step should be finer than regular step",
    );
  });
});

describe("Keyboard Navigation Tests", () => {
  it("should focus next element on Tab", () => {
    const view = new OneWellScreenView(new OneWellModel());
    const focusableElements = view.getFocusableElements();

    // Simulate tab navigation
    focusableElements[0].focus();
    simulateKeyPress("Tab");

    assert.equal(
      document.activeElement,
      focusableElements[1].pdomElement,
      "Tab should move focus to next element",
    );
  });

  it("should activate button on Enter key", () => {
    const model = new OneWellModel();
    const view = new OneWellScreenView(model);
    const playButton = view.findNodeByName("Play/Pause");

    const initialState = model.isPlayingProperty.value;

    playButton.focus();
    simulateKeyPress("Enter");

    assert.notEqual(
      model.isPlayingProperty.value,
      initialState,
      "Enter should toggle play state",
    );
  });
});

describe("Screen Reader Announcement Tests", () => {
  it("should announce reset action", () => {
    utteranceQueue.clear();

    const model = new OneWellModel();
    const view = new OneWellScreenView(model);

    view.reset();

    const alerts = utteranceQueue.queue;
    assert.isTrue(alerts.length > 0);
    assert.include(alerts[0].toString().toLowerCase(), "reset");
  });

  it("should debounce slider value announcements", (done) => {
    utteranceQueue.clear();

    const model = new OneWellModel();

    // Rapidly change mass value
    for (let i = 0; i < 10; i++) {
      model.massProperty.value += 0.01;
    }

    // Should not have announced 10 times
    assert.isTrue(
      utteranceQueue.queue.length < 5,
      "Should debounce rapid changes",
    );

    // Wait for debounce, then check final announcement
    setTimeout(() => {
      assert.isTrue(
        utteranceQueue.queue.length > 0,
        "Should announce after debounce",
      );
      done();
    }, 600);
  });
});
```

#### 9.4 A11y View Testing

Create accessible view HTML for visual PDOM inspection:

**File:** `tests/qppw-a11y-view.html`

```html
<!DOCTYPE html>
<html lang="en">
  <head>
    <meta charset="UTF-8" />
    <title>QPPW A11y View</title>
    <style>
      body {
        display: flex;
        margin: 0;
        font-family: Arial, sans-serif;
      }

      #sim-frame {
        width: 50%;
        height: 100vh;
        border: none;
      }

      #pdom-view {
        width: 50%;
        height: 100vh;
        overflow-y: auto;
        padding: 20px;
        background: #f5f5f5;
        font-family: monospace;
      }

      .pdom-node {
        margin: 10px 0;
        padding: 10px;
        background: white;
        border-left: 3px solid #007bff;
      }

      .pdom-node[aria-hidden="true"] {
        opacity: 0.5;
      }

      .tag-name {
        color: #d73a49;
        font-weight: bold;
      }

      .attribute {
        color: #005cc5;
      }

      .content {
        color: #24292e;
        margin-top: 5px;
      }
    </style>
  </head>
  <body>
    <iframe
      id="sim-frame"
      src="../index.html?supportsInteractiveDescription"
    ></iframe>
    <div id="pdom-view">
      <h2>Parallel DOM Structure</h2>
      <div id="pdom-tree"></div>
    </div>

    <script type="module">
      // Monitor PDOM changes and display
      function updatePDOMView() {
        const iframe = document.getElementById("sim-frame");
        const pdomRoot =
          iframe.contentDocument.querySelector("[data-pdom-root]");

        if (!pdomRoot) {
          setTimeout(updatePDOMView, 100);
          return;
        }

        const treeDiv = document.getElementById("pdom-tree");
        treeDiv.innerHTML = renderPDOMTree(pdomRoot);
      }

      function renderPDOMTree(element, level = 0) {
        let html = "";
        const indent = "  ".repeat(level);

        html += `<div class="pdom-node" style="margin-left: ${level * 20}px">`;
        html += `<span class="tag-name">&lt;${element.tagName.toLowerCase()}&gt;</span>`;

        // Show important attributes
        [
          "id",
          "role",
          "aria-label",
          "aria-labelledby",
          "aria-describedby",
          "aria-live",
        ].forEach((attr) => {
          if (element.hasAttribute(attr)) {
            html += ` <span class="attribute">${attr}="${element.getAttribute(attr)}"</span>`;
          }
        });

        if (element.textContent.trim()) {
          html += `<div class="content">${element.textContent.trim()}</div>`;
        }

        html += "</div>";

        // Recurse for children
        Array.from(element.children).forEach((child) => {
          html += renderPDOMTree(child, level + 1);
        });

        return html;
      }

      // Update every 500ms
      setInterval(updatePDOMView, 500);
    </script>
  </body>
</html>
```

**Success Criteria:**

- ✅ All automated tests pass
- ✅ Manual screen reader testing complete for NVDA, JAWS, VoiceOver
- ✅ Keyboard-only navigation verified
- ✅ A11y-view shows correct PDOM structure
- ✅ No accessibility violations in automated scanning tools
- ✅ User testing with actual AT users (if possible)

---

## Testing & Validation Strategy

### Testing Levels

#### Level 1: Developer Testing (Continuous)

- A11y-view visual inspection during development
- Keyboard navigation testing for each component
- Automated unit tests for PDOM structure
- Browser console accessibility warnings addressed

#### Level 2: QA Testing (Before Release)

- Comprehensive screen reader testing (NVDA, JAWS, VoiceOver, TalkBack)
- Full keyboard-only interaction testing
- Automated accessibility scanning (axe, WAVE)
- Cross-browser compatibility testing

#### Level 3: User Acceptance Testing (Pre-Launch)

- Testing with actual users who rely on assistive technology
- Feedback collection and iteration
- Performance testing with AT enabled
- Documentation review by accessibility experts

### Testing Tools

**Screen Readers:**

- **NVDA** (Windows, free) - Primary development tool
- **JAWS** (Windows, commercial) - Industry standard
- **VoiceOver** (macOS/iOS, built-in) - Apple ecosystem
- **TalkBack** (Android, built-in) - Mobile testing

**Automated Testing:**

- **axe-core** - Accessibility rules engine
- **WAVE** - Web accessibility evaluation tool
- **Lighthouse** - Chrome DevTools accessibility audit
- **Pa11y** - Automated testing dashboard

**Manual Testing:**

- **A11y-view** - Custom PDOM visualization
- **Browser DevTools** - Accessibility tree inspection
- **Keyboard navigation** - Tab, arrow keys, shortcuts

### Success Metrics

| Metric                       | Target                           | Measurement      |
| ---------------------------- | -------------------------------- | ---------------- |
| PDOM coverage                | 100% interactive elements        | Automated scan   |
| Keyboard accessibility       | 100% functionality               | Manual testing   |
| Screen reader compatibility  | Works with NVDA, JAWS, VoiceOver | User testing     |
| WCAG 2.1 Level AA compliance | 0 violations                     | axe-core scan    |
| User task completion         | >90% success rate                | User testing     |
| AT user satisfaction         | >4/5 rating                      | Post-test survey |

---

## Resources & References

### PhET Documentation

- [PhET Interactive Description Technical Guide](https://github.com/phetsims/phet-info/blob/main/doc/interactive-description-technical-guide.md) - Comprehensive implementation guide
- [Scenery Accessibility Documentation](https://phas.ubc.ca/~sqilabs/phetsims/scenery/doc/accessibility.html) - PDOM API reference
- [PhET Accessibility Implementation](https://phet.colorado.edu/vi/accessibility/implementation) - Overview and philosophy
- [PhET Accessibility Research](https://phet.colorado.edu/en/accessibility/research) - Research findings and best practices

### Academic Papers

- [Parallel DOM Architecture for Accessible Interactive Simulations](https://dl.acm.org/doi/10.1145/3192714.3192817) - W4A 2018 conference paper

### W3C Standards

- [WCAG 2.1 Guidelines](https://www.w3.org/WAI/WCAG21/quickref/) - Web Content Accessibility Guidelines
- [ARIA Authoring Practices Guide](https://www.w3.org/WAI/ARIA/apg/) - ARIA design patterns
- [WAI-ARIA 1.2 Specification](https://www.w3.org/TR/wai-aria-1.2/) - Technical specification

### Assistive Technology Resources

- [NVDA User Guide](https://www.nvaccess.org/files/nvda/documentation/userGuide.html)
- [JAWS Keyboard Shortcuts](https://www.freedomscientific.com/training/jaws/hotkeys/)
- [VoiceOver User Guide](https://support.apple.com/guide/voiceover/welcome/mac)
- [TalkBack Help](https://support.google.com/accessibility/android/topic/3529932)

### Testing Tools

- [axe DevTools](https://www.deque.com/axe/devtools/)
- [WAVE Browser Extension](https://wave.webaim.org/extension/)
- [Lighthouse](https://developers.google.com/web/tools/lighthouse)

---

## Appendix

### A. ARIA Attributes Quick Reference

| Attribute          | Purpose                       | Example                                |
| ------------------ | ----------------------------- | -------------------------------------- |
| `aria-label`       | Accessible name               | `aria-label="Particle Mass Slider"`    |
| `aria-labelledby`  | Reference to labeling element | `aria-labelledby="mass-label"`         |
| `aria-describedby` | Reference to description      | `aria-describedby="mass-help"`         |
| `aria-live`        | Live region politeness        | `aria-live="polite"`                   |
| `aria-valuemin`    | Minimum slider value          | `aria-valuemin="0.5"`                  |
| `aria-valuemax`    | Maximum slider value          | `aria-valuemax="1.1"`                  |
| `aria-valuenow`    | Current slider value          | `aria-valuenow="0.8"`                  |
| `aria-valuetext`   | Human-readable value          | `aria-valuetext="0.8 electron masses"` |
| `aria-pressed`     | Toggle button state           | `aria-pressed="true"`                  |
| `aria-checked`     | Checkbox state                | `aria-checked="true"`                  |
| `aria-role`        | Semantic role override        | `aria-role="slider"`                   |
| `aria-hidden`      | Hide from AT                  | `aria-hidden="true"`                   |

### B. Common PDOM Patterns

**Button:**

```typescript
{
  tagName: 'button',
  innerContent: 'Click me',
  accessibleName: 'Descriptive name'
}
```

**Checkbox:**

```typescript
{
  tagName: 'input',
  inputType: 'checkbox',
  labelContent: 'Option label',
  accessibleChecked: booleanProperty
}
```

**Slider:**

```typescript
{
  tagName: 'input',
  inputType: 'range',
  labelContent: 'Parameter name',
  accessibleValueText: derivedProperty,
  ariaValueMin: 0,
  ariaValueMax: 100,
  keyboardStep: 1,
  shiftKeyboardStep: 0.1,
  pageKeyboardStep: 10
}
```

**Heading:**

```typescript
{
  tagName: 'h2',
  innerContent: 'Section Title'
}
```

**Description:**

```typescript
{
  tagName: 'p',
  innerContent: 'Descriptive text'
}
```

**Live Region:**

```typescript
{
  tagName: 'div',
  ariaRole: 'status',
  ariaLive: 'polite',
  innerContent: dynamicProperty
}
```

### C. Potential Type Descriptions Reference

Complete descriptions for all 12 quantum potentials:

1. **Infinite Square Well**: Particle in rigid box with infinite barriers. Exactly solvable. Energy levels: E_n = n²π²ℏ²/(2mL²)

2. **Finite Square Well**: Square well with finite barrier. Wavefunctions decay exponentially outside well. Transcendental equation for energies.

3. **Harmonic Oscillator**: Quadratic potential V(x) = ½mω²x². Equally spaced energy levels: E_n = ℏω(n + ½)

4. **Morse Potential**: V(x) = D(1 - e^(-ax))². Models molecular vibrations. Anharmonic oscillator.

5. **Pöschl-Teller**: V(x) = -V₀sech²(x/a). Exactly solvable hyperbolic potential. Important in scattering theory.

6. **Rosen-Morse**: V(x) = -V₀tanh²(x/a) + V₁tanh(x/a). Variant of Pöschl-Teller.

7. **Eckart Potential**: Barrier potential for tunneling studies. Used in molecular physics.

8. **Asymmetric Triangle**: V(x) = F·x inside well (electric field). Airy function solutions.

9. **Triangular Potential**: V(x) = |x|. V-shaped well. Related to Airy functions.

10. **Coulomb 1D**: V(x) = -k/|x|. One-dimensional hydrogen-like atom. E_n ∝ 1/n²

11. **Coulomb 3D**: V(r) = -k/r. Full hydrogen atom with angular momentum. E_n = -13.6eV/n²

12. **Double Square Well**: Two square wells with barrier. Demonstrates tunneling and level splitting.

### D. Physics Concepts for Descriptions

**Key quantum concepts to explain accessibly:**

- **Wavefunction (ψ)**: Mathematical description of quantum state. Contains all information about particle.

- **Probability Density (|ψ|²)**: Shows where particle is likely to be found upon measurement.

- **Energy Levels**: Discrete allowed energies for bound states. Quantum number n = 1, 2, 3...

- **Nodes**: Points where wavefunction crosses zero. Number of nodes = n - 1

- **RMS (Root Mean Square)**: Measure of spread or uncertainty. Position uncertainty Δx.

- **Heisenberg Uncertainty**: Δx · Δp ≥ ℏ/2. Cannot know position and momentum precisely.

- **Tunneling**: Quantum phenomenon where particle penetrates classically forbidden regions.

- **Superposition**: Linear combination of eigenstates. Enables time evolution and interference.

- **Time Evolution**: ψ(t) = ψ(0)e^(-iEt/ℏ). Phase rotation at frequency E/ℏ.

### E. Implementation Timeline Estimate

| Phase                        | Duration      | Dependencies | Team Size               |
| ---------------------------- | ------------- | ------------ | ----------------------- |
| Phase 0: Foundation          | 1 day         | None         | 1 developer             |
| Phase 1: Screen Structure    | 2-3 days      | Phase 0      | 1 developer             |
| Phase 2: Basic Controls      | 3-4 days      | Phase 1      | 1 developer             |
| Phase 3: Sliders             | 2-3 days      | Phase 2      | 1 developer             |
| Phase 4: Complex Controls    | 4-5 days      | Phase 3      | 1-2 developers          |
| Phase 5: Chart Descriptions  | 3-4 days      | Phase 4      | 1 developer             |
| Phase 6: Interactive Tools   | 3-4 days      | Phase 5      | 1 developer             |
| Phase 7: Live Alerts         | 2-3 days      | Phases 1-6   | 1 developer             |
| Phase 8: Keyboard Shortcuts  | 2-3 days      | Phases 1-6   | 1 developer             |
| Phase 9: Testing             | 4-5 days      | All phases   | 2-3 testers             |
| **Total Estimated Duration** | **6-8 weeks** | -            | **1-2 developers + QA** |

**Note:** Timeline assumes experienced developer familiar with PhET/Scenery stack. First-time implementation may take 50-100% longer.

---

## Conclusion

Implementing comprehensive accessibility for QPPW is a substantial but achievable undertaking. The PhET Scenery PDOM architecture provides a robust foundation, and the simulation's clean architecture supports accessibility enhancement without major refactoring.

**Key Success Factors:**

1. Follow PhET's established patterns and best practices
2. Prioritize keyboard navigation and screen reader support
3. Provide meaningful, physics-rich descriptions
4. Test extensively with actual assistive technology
5. Iterate based on user feedback

**Expected Impact:**

- Enable access for blind and low-vision students
- Support keyboard-only users
- Improve usability for all users
- Demonstrate commitment to inclusive education
- Set foundation for future accessibility work

**Next Steps:**

1. Review and approve this implementation plan
2. Begin Phase 0 (foundation setup)
3. Implement incrementally, testing each phase
4. Gather user feedback throughout
5. Document lessons learned for future projects

---

**Document Version:** 1.0
**Last Updated:** November 30, 2025
**Author:** Claude (AI Assistant)
**Status:** Awaiting Review & Approval
