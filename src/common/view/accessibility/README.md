# Accessibility Components

This directory contains accessibility-related components for the QPPW simulation, implementing PhET's Parallel DOM (PDOM) architecture.

## Planned Components

Based on the [Accessibility Implementation Plan](../../../../docs/ACCESSIBILITY_IMPLEMENTATION_PLAN.md), this directory will contain:

- **QPPWDescriber.ts** - Generates accessible descriptions for quantum states and visualizations
- **QPPWAlerter.ts** - Manages live announcements to screen readers for state changes
- **ScreenSummaryNode.ts** - Provides dynamic screen summaries for each simulation screen
- **QPPWKeyboardHelpContent.ts** - Keyboard shortcuts reference and help dialog

## Purpose

These components enable:
- Screen reader compatibility (NVDA, JAWS, VoiceOver, TalkBack)
- Full keyboard navigation
- Dynamic content descriptions
- Live announcements for state changes
- Semantic HTML structure alongside canvas rendering

## Implementation Status

🚧 **Planning Phase** - Directory structure created. Components to be implemented according to the phased approach outlined in the implementation plan.

See [ACCESSIBILITY_IMPLEMENTATION_PLAN.md](../../../../docs/ACCESSIBILITY_IMPLEMENTATION_PLAN.md) for full details.
