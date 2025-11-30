# Accessibility

## Overview

The Quantum Bound States Simulation (QPPW) is committed to providing an accessible educational experience for all students, including those who use assistive technologies such as screen readers or rely on keyboard navigation.

## Current Status

🚀 **Active Development - Phases 1 & 2 Complete**

We have completed the following:

### Phase 0: Foundation ✅ Complete
- ✅ **Interactive Description Support Enabled** - Added `supportsInteractiveDescription: true` to package.json
- ✅ **Directory Structure Created** - Established `src/common/view/accessibility/` for accessibility components
- ✅ **Implementation Plan Developed** - Comprehensive technical roadmap documented

### Phase 1: Screen Structure ✅ Complete
- ✅ **Three-Section PDOM Layout** - Implemented for all four screens (Intro, One Well, Two Wells, Many Wells)
- ✅ **Screen Summary Components** - Dynamic descriptions that update with model state
- ✅ **Navigation Order** - Logical tab order: Screen Summary → Play Area → Control Area

### Phase 2: Basic Interactive Controls ✅ Complete
- ✅ **Play/Pause Button** - Accessible with keyboard, announces state changes
- ✅ **Reset Button** - PDOM support with screen reader announcements
- ✅ **Checkboxes** - Wave function views, classical probability (keyboard accessible)
- ✅ **Radio Button Groups** - Display mode, animation speed (arrow key navigation)
- ✅ **Live Announcements** - Screen reader feedback for energy level, potential type, and parameter changes

## Planned Features

Our accessibility implementation will include:

### Screen Reader Support
- **Compatible with**: NVDA, JAWS, VoiceOver, TalkBack
- **Parallel DOM (PDOM)**: Semantic HTML structure alongside visual canvas
- **Dynamic Descriptions**: Real-time updates for quantum state changes
- **Live Announcements**: Screen reader feedback for user actions

### Keyboard Navigation
- **Full Keyboard Access**: All interactive elements controllable via keyboard
- **Keyboard Shortcuts**: Quick access to common actions (Play/Pause, Reset, etc.)
- **Focus Management**: Clear visual focus indicators and logical tab order
- **Help Dialog**: Comprehensive keyboard shortcuts reference

### Accessible Content
- **Screen Summaries**: Overview of current simulation state
- **Chart Descriptions**: Textual descriptions of energy levels, wavefunctions, and momentum distributions
- **Control Labels**: Clear, descriptive labels for all UI controls
- **Physics Context**: Accessible explanations of quantum concepts

## Implementation Phases

Our accessibility implementation follows a phased approach:

1. **Phase 0: Foundation** (✅ Complete)
   - Enable interactive description support
   - Create directory structure
   - Set up development environment

2. **Phase 1: Screen Structure** (✅ Complete)
   - Implement three-section PDOM layout
   - Add screen summaries for all four screens
   - Establish navigation order

3. **Phase 2: Basic Interactive Controls** (✅ Complete)
   - Make buttons, checkboxes, and radio buttons accessible
   - Add PDOM attributes and help text
   - Implement screen reader announcements

4. **Phase 3: Sliders** (🚧 Next)
   - Add keyboard navigation to sliders
   - Implement accessible value announcements
   - Configure keyboard step sizes

4. **Phase 4: Complex Components** (Planned)
   - Energy level selection via keyboard
   - Accessible potential type dropdown
   - Interactive tool accessibility

5. **Phase 5: Visualizations** (Planned)
   - Accessible descriptions for energy charts
   - Wavefunction statistical summaries
   - Momentum distribution descriptions

6. **Phase 6: Testing & Validation** (Planned)
   - Screen reader testing
   - Keyboard-only navigation verification
   - User acceptance testing

## Technical Details

Our accessibility implementation uses the **PhET Parallel DOM (PDOM)** architecture, which:

- Maintains semantic HTML alongside visual canvas rendering
- Provides screen reader users with equivalent information
- Enables keyboard users to access all functionality
- Keeps visual and accessible representations synchronized

### Key Technologies
- **Scenery AccessibleNode**: Core accessibility properties
- **UtteranceQueue**: Manages screen reader announcements
- **KeyboardDragListener**: Keyboard interaction with draggable elements
- **ARIA Attributes**: Semantic roles and states for complex widgets

## Documentation

For detailed technical information, see:

- **[Accessibility Implementation Plan](docs/ACCESSIBILITY_IMPLEMENTATION_PLAN.md)** - Complete technical specification with code examples and implementation timeline
- **[Accessibility Components README](src/common/view/accessibility/README.md)** - Directory structure and component descriptions

## Standards & Guidelines

This implementation follows:

- **WCAG 2.1 Level AA** - Web Content Accessibility Guidelines
- **WAI-ARIA 1.2** - Accessible Rich Internet Applications
- **PhET Accessibility Best Practices** - Proven patterns from PhET Interactive Simulations

## Contributing

We welcome contributions to improve accessibility! If you:

- Use assistive technology and want to provide feedback
- Have expertise in accessibility implementation
- Encounter accessibility barriers
- Want to help with testing

Please see [CONTRIBUTE.md](CONTRIBUTE.md) for contribution guidelines or [open an issue](https://github.com/veillette/QPPW/issues) with your feedback.

## Timeline

- **Phase 0 (Foundation)**: ✅ Complete - November 2025
- **Phase 1 (Screen Structure)**: ✅ Complete - November 2025
- **Phase 2 (Basic Controls)**: ✅ Complete - November 2025
- **Phase 3 (Sliders)**: 🚧 In Progress
- **Phases 4-5 (Advanced Features)**: Planned Q1 2026
- **Phase 6 (Testing & Launch)**: Planned Q1-Q2 2026

## Contact

For accessibility-specific questions or feedback, please [open an issue](https://github.com/veillette/QPPW/issues) with the "accessibility" label.

---

*Last Updated: November 30, 2025*
