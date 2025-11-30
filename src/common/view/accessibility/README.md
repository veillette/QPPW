# Accessibility Components

This directory contains components and utilities for making the QPPW simulation accessible to users with disabilities, following the PhET Parallel DOM (PDOM) architecture.

## Overview

The accessibility implementation enables:
- **Screen reader support** (NVDA, JAWS, VoiceOver, TalkBack)
- **Keyboard navigation** for all interactive elements
- **Live announcements** for state changes
- **Semantic HTML** structure alongside canvas rendering

## Components

### ScreenSummaryNode.ts

Creates the accessible screen summary following PhET PDOM patterns. This is the first section of the three-section PDOM layout.

**Features:**
- Dynamic overview of simulation state
- Updates automatically when model properties change
- Includes potential type, energy level, parameters, and statistics

**Usage:**
```typescript
const screenSummaryNode = new ScreenSummaryNode(model, {
  screenName: "One Well",
  screenDescription: "Explore quantum bound states in a single potential well."
});
```

### QPPWDescriber.ts

Provides accessible descriptions for quantum physics concepts.

**Features:**
- Potential type descriptions with physics context
- Superposition type explanations
- Display mode descriptions
- Announcement text generation
- Help text for controls

**Usage:**
```typescript
const description = QPPWDescriber.getPotentialTypeDescription(
  PotentialType.HARMONIC_OSCILLATOR
);
// Returns: "Quadratic potential resembling mass on spring. Energy levels uniformly spaced."

const announcement = QPPWDescriber.createEnergyLevelAnnouncement(2, 1.234, 5);
// Returns: "Selected energy level 3 of 5. Energy: 1.234 electron volts. Wavefunction has 2 nodes."
```

### QPPWAlerter.ts

Manages live announcements to screen readers using the UtteranceQueue.

**Features:**
- Automatic alerts for energy level changes
- Potential type change announcements
- Playback state announcements
- Debounced parameter change alerts
- Configurable alert priorities

**Usage:**
```typescript
const alerter = new QPPWAlerter(model);
// Automatically sets up listeners for model changes

// Manual announcements
alerter.alertResetAll();
```

## Implementation Status

### Phase 0: Foundation ✅ Complete
- ✅ Interactive description support enabled
- ✅ Directory structure created
- ✅ Development environment configured

### Phase 1: Screen Structure ✅ Complete
- ✅ Three-section PDOM layout implemented
- ✅ Screen summaries for all four screens
- ✅ PDOM navigation order established

### Phase 2: Interactive Controls - Basic Elements ✅ Complete
- ✅ Play/Pause button with PDOM support
- ✅ Reset button with announcements
- ✅ Checkboxes (wave function views, classical probability)
- ✅ Radio button groups (display mode, animation speed)

### Phase 3: Sliders 🚧 Planned
- Particle mass slider
- Well width slider
- Well depth slider
- Other parameter sliders

### Future Phases
- Phase 4: Complex components (combo boxes, energy level selection)
- Phase 5: Chart descriptions
- Phase 6: Interactive tools
- Phase 7: Live alerts (partially complete)
- Phase 8: Keyboard shortcuts
- Phase 9: Testing & validation

## Architecture

### Three-Section PDOM Layout

All screens follow this standard structure:

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

### Key Accessibility Features

**Semantic HTML:**
- Proper heading hierarchy (h2 for sections)
- Form controls with labels
- ARIA attributes where needed

**Keyboard Navigation:**
- Tab order follows logical flow
- All interactive elements keyboard accessible
- Consistent keyboard shortcuts

**Screen Reader Support:**
- Dynamic content descriptions
- Live region announcements
- Context-sensitive help text

## Testing

### Manual Testing
- Use NVDA (Windows) or VoiceOver (macOS) for screen reader testing
- Navigate with keyboard only (no mouse)
- Verify all controls are accessible and announced correctly

### A11y View
View the PDOM structure by running the simulation with `?supportsInteractiveDescription` query parameter.

## References

- [PhET Accessibility Guide](https://github.com/phetsims/phet-info/blob/main/doc/interactive-description-technical-guide.md)
- [WCAG 2.1 Guidelines](https://www.w3.org/WAI/WCAG21/quickref/)
- [ARIA Authoring Practices](https://www.w3.org/WAI/ARIA/apg/)
- [Accessibility Implementation Plan](../../../../docs/ACCESSIBILITY_IMPLEMENTATION_PLAN.md)

## Contributing

When adding new interactive elements:

1. **Always add PDOM properties:**
   - `labelContent` or `accessibleName`
   - `helpText` for complex controls
   - `descriptionContent` for additional context

2. **Test with screen readers:**
   - Verify announcements are clear and concise
   - Check navigation order makes sense
   - Ensure state changes are announced

3. **Follow PhET patterns:**
   - Use established PDOM structures
   - Maintain consistent language
   - Keep descriptions physics-focused

## License

Same as parent project.
