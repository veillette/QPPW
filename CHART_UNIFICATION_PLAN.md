# Chart Unification Plan - Detailed Analysis & Recommendations

## Executive Summary

Based on analysis of the architecture branch, **significant progress has been made** (Phases 1-2 complete). This document outlines **remaining opportunities** for further simplification and unification of the three chart components.

**Current State:**
- **Total lines**: 3,347 lines across 3 chart files (down from 3,690 - **9% reduction achieved**)
- **BaseChartNode**: 245 lines ✅ (unified base class)
- **EnergyChartNode**: 1,032 lines ✅ (reduced from 1,464 via PotentialRenderer extraction)
- **WaveFunctionChartNode**: 1,532 lines ⚠️ (largest, needs attention)
- **WavenumberChartNode**: 538 lines ✅ (smallest, well-structured)

---

## What's Been Accomplished ✅

### Phase 1: Base Class Unification (COMPLETE)
- ✅ All three charts now extend `BaseChartNode`
- ✅ Common coordinate transformation methods extracted
- ✅ Shared chart setup (margins, transform, background, clipping) unified
- ✅ **Estimated reduction: 300-400 lines**

### Phase 2: Renderer Strategy Pattern (COMPLETE)
- ✅ Extracted 13 potential renderers into separate files
- ✅ `EnergyChartNode` now uses registry pattern for potential rendering
- ✅ Old 520-line `updatePotentialCurve()` method replaced with ~40 lines
- ✅ **Estimated reduction: 450 lines**

### Tool System Infrastructure (PARTIAL)
- ✅ `BaseChartTool` abstract class created (250 lines)
- ✅ `ChartToolRegistry` for managing tools (154 lines)
- ✅ `BaseVisualization` for passive visualizations (64 lines)
- ⚠️ Only used in `WaveFunctionChartNode` - not leveraged in other charts

---

## Current Problems & Opportunities

### 1. ⚠️ **WaveFunctionChartNode is Monolithic** (1,532 lines - CRITICAL)

**Problem**: Single class handles **three different display modes** with significant behavioral differences:

```typescript
// Three distinct display modes in one class:
- probabilityDensity: Shows |ψ|², RMS indicators, classical probability
- waveFunction: Shows real/imaginary/magnitude parts, time evolution
- phaseColor: Shows phase-colored magnitude, different coloring logic
```

**Current responsibilities** (too many for one class):
1. Probability density plotting with filled area
2. Wavefunction component plotting (real, imag, magnitude)
3. Phase color visualization with complex coloring
4. Superposition state handling with time evolution
5. RMS position calculations and indicator display
6. Tool management (area, curvature, derivative)
7. Zeros visualization
8. Classical probability overlay
9. Dynamic Y-axis range calculations (different for each mode)
10. State label formatting (different for each mode)
11. Accessibility descriptions (different for each mode)

**Key observation**: The class has multiple `if (displayMode === "...")` conditionals throughout:
- `updateViewRange()`: Lines 919-933 (different ranges for each mode)
- `updateWaveFunction()`: Lines 1215-1328 (completely different rendering)
- `updateSuperpositionWavefunction()`: Lines 1064-1173 (completely different rendering)
- `plotWaveFunctionComponents()`: Lines 1334-1397 (only for waveFunction mode)
- `plotProbabilityDensityFromArray()`: Lines 1402-1439 (only for probabilityDensity mode)
- Plus tool visibility logic, label updates, etc.

**Should we split it?**

**Option A: Split into 3 separate chart classes** ❌ NOT RECOMMENDED
```
- ProbabilityDensityChartNode (focused on |ψ|²)
- WaveFunctionChartNode (focused on ψ components)
- PhaseColorChartNode (focused on phase visualization)
```
**Pros:**
- Each class would be ~500 lines (much more manageable)
- Clear separation of concerns
- Easier to test and maintain individual modes

**Cons:**
- Need separate constructors and initialization for each
- More files to manage (though organized)
- Some duplication of common wavefunction data fetching
- Need to handle switching between chart types in parent view

**Option B: Extract display strategies (composition over splitting)** ✅ RECOMMENDED
```typescript
// Keep one chart, but extract display strategies
interface DisplayStrategy {
  updateViewRange(data: WavefunctionData): Range;
  render(data: WavefunctionData, context: RenderContext): void;
  updateLabels(data: WavefunctionData): void;
}

class ProbabilityDensityStrategy implements DisplayStrategy { ... }
class WaveFunctionStrategy implements DisplayStrategy { ... }
class PhaseColorStrategy implements DisplayStrategy { ... }
```
**Pros:**
- Reduces complexity without file proliferation
- Uses proven strategy pattern (like PotentialRenderer)
- Easier mode switching (just swap strategy)
- Common code stays in main chart class

**Cons:**
- Still one large file (but more organized)
- Need to design strategy interface carefully

**Option C: Do nothing** ❌ NOT RECOMMENDED
- WaveFunctionChartNode at 1,532 lines is hard to maintain
- Adding new features or modes becomes increasingly difficult

**RECOMMENDATION**: **Option B (Display Strategy Pattern)** - estimated **400-500 line reduction**

---

### 2. 🔄 **Duplicated Manual Grid Line Code** (MEDIUM PRIORITY)

**Problem**: Both `EnergyChartNode` and `WaveFunctionChartNode` have nearly identical manual grid line creation:

**EnergyChartNode.ts (lines 382-420):**
```typescript
// Manual Y-axis grid lines (GridLineSet causes hang)
for (let energy = yMin; energy <= yMax; energy += 5) {
  const y = this.chartMargins.top + this.chartTransform.modelToViewY(energy);
  const gridLine = new Line(...);
  axesNode.addChild(gridLine);
}

// Manual X-axis grid lines (GridLineSet causes hang)
for (let pos = -X_AXIS_RANGE_NM; pos <= X_AXIS_RANGE_NM; pos += 2) {
  const x = this.chartMargins.left + this.chartTransform.modelToViewX(pos);
  const gridLine = new Line(...);
  axesNode.addChild(gridLine);
}
```

**WaveFunctionChartNode.ts (lines 409-427):**
```typescript
// Manual X-axis grid lines (GridLineSet causes hang)
for (let pos = -X_AXIS_RANGE_NM; pos <= X_AXIS_RANGE_NM; pos += 2) {
  const x = this.chartMargins.left + this.chartTransform.modelToViewX(pos);
  const gridLine = new Line(...);
  axesNode.addChild(gridLine);
}
```

**Duplication**: ~80 lines of nearly identical code

**Solution**: Extract to `BaseChartNode` utility methods:
```typescript
// In BaseChartNode:
protected createVerticalGridLines(spacing: number, range?: Range): Line[] {
  const yMin = range?.min ?? this.yMinProperty.value;
  const yMax = range?.max ?? this.yMaxProperty.value;
  const lines: Line[] = [];

  for (let y = yMin; y <= yMax; y += spacing) {
    const viewY = this.chartMargins.top + this.chartTransform.modelToViewY(y);
    lines.push(new Line(
      this.chartMargins.left,
      viewY,
      this.chartMargins.left + this.plotWidth,
      viewY,
      { stroke: QPPWColors.gridLineProperty, lineWidth: 1 }
    ));
  }
  return lines;
}

protected createHorizontalGridLines(spacing: number, range?: Range): Line[] {
  // Similar implementation for X-axis
}
```

**Expected reduction**: ~60-80 lines (removed from subclasses, ~30 added to base)

---

### 3. 🔄 **Common Axis Creation Patterns** (LOW-MEDIUM PRIORITY)

**Problem**: All three charts create similar axis structures but with slight variations:

**Pattern observed across all charts:**
```typescript
protected createAxes(): Node {
  const axesNode = new Node();

  // Y-axis at left edge
  const yAxisLeftNode = new AxisLine(this.chartTransform, Orientation.VERTICAL, {...});

  // Y-axis at origin (center)
  const yAxisNode = new AxisLine(this.chartTransform, Orientation.VERTICAL, { value: 0, opacity: 0.3 });

  // X-axis at bottom
  const xAxisBottomNode = new AxisLine(...);

  // Tick marks
  const yTickMarksNode = new TickMarkSet(...);
  const xTickMarksNode = new TickMarkSet(...);

  // Tick labels
  const yTickLabelsNode = new TickLabelSet(...);
  const xTickLabelsNode = new TickLabelSet(...);

  // Axis labels
  const yLabelText = new Text(...);
  const xLabelText = new Text(...);

  return axesNode;
}
```

**Commonality**: ~70% of axis creation code is identical

**Solution**: Create `ChartAxisBuilder` utility class:
```typescript
class ChartAxisBuilder {
  constructor(private chartTransform: ChartTransform, private margins: ChartMargins) {}

  createStandardYAxis(position: 'left' | 'center', opacity?: number): AxisLine { ... }
  createStandardXAxis(position: 'bottom' | 'center'): AxisLine { ... }
  createTickMarks(orientation: Orientation, spacing: number, edge: string): TickMarkSet { ... }
  createTickLabels(
    orientation: Orientation,
    spacing: number,
    formatter?: (value: number) => string
  ): TickLabelSet { ... }
  createAxisLabel(text: string, orientation: 'horizontal' | 'vertical', position: {x: number, y: number}): Text { ... }
}
```

Usage in subclasses becomes much simpler:
```typescript
protected createAxes(): Node {
  const builder = new ChartAxisBuilder(this.chartTransform, this.chartMargins);
  const axesNode = new Node();

  // Add grid lines (using new base class methods)
  axesNode.children = [
    ...this.createHorizontalGridLines(2),
    ...this.createVerticalGridLines(5)
  ];

  // Add axes
  axesNode.addChild(builder.createStandardYAxis('left'));
  axesNode.addChild(builder.createStandardYAxis('center', 0.3));
  axesNode.addChild(builder.createStandardXAxis('bottom'));

  // Add tick marks and labels
  axesNode.addChild(builder.createTickMarks(Orientation.VERTICAL, 5, 'min'));
  axesNode.addChild(builder.createTickLabels(Orientation.VERTICAL, 5, v => v.toFixed(0)));

  // Add axis labels
  axesNode.addChild(builder.createAxisLabel(this.yAxisLabelText, 'vertical', {x: 15, y: chartHeight/2}));
  axesNode.addChild(builder.createAxisLabel("Position (nm)", 'horizontal', {x: chartWidth/2, y: chartHeight-15}));

  return axesNode;
}
```

**Expected reduction**: ~100-150 lines per chart × 3 charts = ~300-450 lines total
(accounting for new builder class ~150 lines)

---

### 4. 🔄 **RMS Indicator Display Logic** (LOW PRIORITY)

**Problem**: Both `WaveFunctionChartNode` and `WavenumberChartNode` have nearly identical RMS indicator display logic:

**Common pattern:**
```typescript
// Calculate RMS statistics
const { avg, rms } = calculateRMSStatistics(xGrid, probabilityDensity);

// Show labels if enabled
if (this.shouldShowRMSIndicators()) {
  this.avgPositionLabel.string = `⟨x⟩ = ${avg.toFixed(2)} nm`;
  this.rmsPositionLabel.string = `Δx = ${rms.toFixed(2)} nm`;

  // Draw double arrow indicator
  const x1 = this.dataToViewX(avg - rms);
  const x2 = this.dataToViewX(avg + rms);
  const y = this.dataToViewY(this.yMaxProperty.value * 0.8);
  this.rmsIndicator.shape = createDoubleArrowShape(x1, x2, y);
} else {
  this.avgLabel.string = "";
  this.rmsLabel.string = "";
  this.rmsIndicator.shape = null;
}
```

**Duplication**: ~80 lines across two charts

**Solution**: Extract to shared utility or base class method:
```typescript
// In BaseChartNode or new RMSIndicatorHelper utility:
protected updateRMSIndicators(
  avg: number,
  rms: number,
  yPosition: number,
  labels: { avgLabel: Text, rmsLabel: Text },
  indicator: Path,
  showProperty: BooleanProperty,
  labelPrefix: string  // "⟨x⟩" or "⟨k⟩"
): void {
  if (showProperty.value) {
    labels.avgLabel.string = `${labelPrefix} = ${avg.toFixed(2)}`;
    labels.rmsLabel.string = `Δ${labelPrefix.slice(1, -1)} = ${rms.toFixed(2)}`;

    const x1 = this.dataToViewX(avg - rms);
    const x2 = this.dataToViewX(avg + rms);
    const y = this.dataToViewY(yPosition);
    indicator.shape = createDoubleArrowShape(x1, x2, y);
  } else {
    labels.avgLabel.string = "";
    labels.rmsLabel.string = "";
    indicator.shape = null;
  }
}
```

**Expected reduction**: ~60 lines (small but eliminates duplication)

---

### 5. 🔄 **Common Plotting Utilities** (LOW PRIORITY)

**Problem**: Both `WaveFunctionChartNode` and `WavenumberChartNode` plot probability-like distributions with filled areas:

**Common pattern:**
```typescript
private plotProbabilityDensityFromArray(xGrid: number[], data: number[]): void {
  const shape = new Shape();
  const points: {x: number, y: number}[] = [];

  // Build points
  for (let i = 0; i < xGrid.length; i++) {
    points.push({
      x: this.dataToViewX(xGrid[i] * M_TO_NM),
      y: this.dataToViewY(data[i])
    });
  }

  // Create filled area
  const y0 = this.dataToViewY(0);
  shape.moveTo(points[0].x, y0);
  shape.lineTo(points[0].x, points[0].y);
  for (let i = 0; i < points.length - 1; i++) {
    shape.lineTo(points[i].x, points[i].y);
  }
  shape.lineTo(points[points.length-1].x, points[points.length-1].y);
  shape.lineTo(points[points.length-1].x, y0);
  shape.close();

  this.path.shape = shape;
}
```

**Duplication**: ~50 lines across two charts

**Solution**: Extract to `BaseChartNode`:
```typescript
protected createFilledCurveShape(
  xGrid: number[],
  yData: number[],
  xScale: number = 1  // e.g., M_TO_NM conversion
): Shape {
  // Implementation of common filled curve logic
}
```

**Expected reduction**: ~30-40 lines

---

## Should WaveFunctionChartNode Be Split? 🤔

### Analysis of WaveFunctionChartNode (1,532 lines)

**The question**: Does this chart do "multiple charts" that should be separated?

**Answer**: **Conceptually yes, practically no (with caveat)**

**Why conceptually yes:**
- It has 3 **very different display modes** that are almost like separate visualizations:
  1. **Probability Density** (`probabilityDensity`): Shows |ψ|² with RMS indicators
  2. **Wavefunction Components** (`waveFunction`): Shows Re(ψ), Im(ψ), |ψ| separately
  3. **Phase Color** (`phaseColor`): Shows color-coded phase visualization

- Each mode has:
  - Different Y-axis ranges and labels
  - Different visual elements (different paths shown/hidden)
  - Different update logic
  - Different tool behaviors

**Why practically no (mostly):**
- All modes share the **same X-axis** (position space)
- All modes share the **same model data** (wavefunction from boundStates)
- All modes share **common tools** (area, curvature, derivative)
- All modes share **common infrastructure** (axes, grid, coordinate transforms)
- Switching between modes is a core feature (user expectation)

**THE BEST APPROACH: Strategy Pattern (Option B)**

Instead of splitting into separate chart classes, extract the **mode-specific logic** into strategy objects:

```typescript
// WaveFunctionChartNode.ts (main class - reduced from 1,532 to ~800 lines)
export class WaveFunctionChartNode extends BaseChartNode {
  private displayStrategy: WaveFunctionDisplayStrategy;

  constructor(...) {
    // Common initialization
    this.displayStrategy = this.createDisplayStrategy(displayMode);
  }

  private update(): void {
    const data = this.getWavefunctionData();

    // Delegate mode-specific logic to strategy
    this.displayStrategy.updateViewRange(data);
    this.displayStrategy.render(data);
    this.displayStrategy.updateIndicators(data);

    // Common logic stays here
    this.updateZeroLine();
    this.updateTools();
  }
}

// New files (3 strategy files, ~200 lines each):
// - ProbabilityDensityDisplayStrategy.ts
// - WaveFunctionComponentsDisplayStrategy.ts
// - PhaseColorDisplayStrategy.ts
```

**Benefits:**
1. **Main chart file reduced by 400-600 lines** (to ~800-1100 lines)
2. **Each display mode isolated** in its own file (~200-250 lines each)
3. **Easy to add new modes** (just add new strategy)
4. **Common code stays shared** (tools, axes, base setup)
5. **Mode switching is trivial** (just swap strategy instance)
6. **Testing is easier** (can test each strategy independently)

**This is exactly like the PotentialRenderer pattern that worked so well for EnergyChartNode!**

---

## Recommended Implementation Plan

### Phase 3: Extract Grid Line and Axis Utilities (QUICK WIN)
**Effort**: 4-6 hours | **Impact**: Medium | **Risk**: Low

#### Actions:
1. **Add grid line methods to BaseChartNode** (~1 hour)
   ```typescript
   protected createHorizontalGridLines(spacing: number, range?: Range): Line[]
   protected createVerticalGridLines(spacing: number, range?: Range): Line[]
   ```

2. **Create ChartAxisBuilder utility class** (~2 hours)
   - Location: `src/common/view/ChartAxisBuilder.ts`
   - Methods for standard axes, tick marks, tick labels, axis labels

3. **Refactor chart createAxes() methods to use utilities** (~1-2 hours)
   - Update EnergyChartNode, WaveFunctionChartNode, WavenumberChartNode

4. **Test all three charts** (~1 hour)
   - Verify visual appearance unchanged
   - Check that dynamic range updates still work

**Expected outcome**:
- **~300-400 line reduction** across all charts
- **Eliminates ~80 lines of duplicated grid line code**
- **Simplifies axis creation** to ~20-30 lines per chart (from ~180)

---

### Phase 4: Display Strategy Pattern for WaveFunctionChartNode (MAJOR WIN)
**Effort**: 12-16 hours | **Impact**: High | **Risk**: Medium

#### Actions:

1. **Design strategy interface** (~2 hours)
   ```typescript
   // src/common/view/wavefunction-strategies/WaveFunctionDisplayStrategy.ts
   export interface WaveFunctionDisplayStrategy {
     updateViewRange(data: WavefunctionData, chartTransform: ChartTransform): void;
     render(data: WavefunctionData, context: RenderContext): void;
     updateIndicators(data: WavefunctionData, showRMS: boolean): void;
     updateLabels(data: WavefunctionData): void;
     getYAxisLabel(): string;
     shouldShowTools(): boolean;
   }
   ```

2. **Extract ProbabilityDensityDisplayStrategy** (~3 hours)
   - File: `src/common/view/wavefunction-strategies/ProbabilityDensityDisplayStrategy.ts`
   - Extract lines 1065-1113 and 1215-1264 from WaveFunctionChartNode
   - Handle probability density plotting
   - Handle RMS indicator display
   - Handle classical probability overlay

3. **Extract WaveFunctionComponentsDisplayStrategy** (~3 hours)
   - File: `src/common/view/wavefunction-strategies/WaveFunctionComponentsDisplayStrategy.ts`
   - Extract lines 1152-1173 and 1308-1397 from WaveFunctionChartNode
   - Handle real/imaginary/magnitude plotting
   - Handle time evolution phase
   - Handle component visibility toggles

4. **Extract PhaseColorDisplayStrategy** (~3 hours)
   - File: `src/common/view/wavefunction-strategies/PhaseColorDisplayStrategy.ts`
   - Extract lines 1129-1151 and 1279-1306 from WaveFunctionChartNode
   - Handle phase color visualization
   - Handle magnitude rendering

5. **Refactor WaveFunctionChartNode to use strategies** (~2-3 hours)
   - Create strategy based on display mode
   - Delegate mode-specific logic to strategy
   - Keep common logic (tools, axes, model linking) in main class

6. **Test all three display modes thoroughly** (~1-2 hours)
   - Test mode switching
   - Test superposition vs single state
   - Test time evolution
   - Test tool interactions

**Expected outcome**:
- **WaveFunctionChartNode reduced from 1,532 to ~800-900 lines**
- **Three strategy files of ~200-250 lines each** (more manageable)
- **Improved maintainability** (each mode isolated)
- **Easier to add new modes** in future

---

### Phase 5: RMS Indicator Unification (POLISH)
**Effort**: 2-3 hours | **Impact**: Low | **Risk**: Low

#### Actions:
1. **Create RMSIndicatorHelper utility** (~1 hour)
   - Location: `src/common/view/RMSIndicatorHelper.ts` OR add to `BaseChartNode`

2. **Refactor WaveFunctionChartNode and WavenumberChartNode** (~1 hour)
   - Replace duplicated RMS display logic with helper calls

3. **Test RMS indicators** (~30 minutes)

**Expected outcome**:
- **~60 line reduction** (eliminates duplication)
- **Consistent RMS indicator behavior** across charts

---

## Summary of Recommended Phases

| Phase | Effort | Impact | Risk | Line Reduction | Priority |
|-------|--------|--------|------|----------------|----------|
| **Phase 1** | Done ✅ | High | Low | -300 to -400 | - |
| **Phase 2** | Done ✅ | High | Low | -450 | - |
| **Phase 3** | 4-6h | Medium | Low | -300 to -400 | **HIGH** |
| **Phase 4** | 12-16h | High | Medium | -400 to -600 | **HIGH** |
| **Phase 5** | 2-3h | Low | Low | -60 | LOW |
| **TOTAL** | 18-25h | **Very High** | Low-Med | **-1,200 to -1,500** | - |

---

## Final Expected Impact

| Metric | Before (Original) | After Phase 2 (Current) | After All Phases (Projected) | Total Improvement |
|--------|-------------------|------------------------|------------------------------|-------------------|
| **Total Lines** | 3,690 | 3,347 | **~2,200-2,500** | **-32% to -40%** |
| **Largest File** | 1,610 | 1,532 | **~900** | **-44%** |
| **Code Duplication** | High | Medium | **Minimal** | 🎯 |
| **Testability** | Poor | Good | **Excellent** | 🎯 |
| **Maintainability** | Poor | Good | **Excellent** | 🎯 |

---

## Decision: Should We Split WaveFunctionChartNode?

### Final Recommendation: ✅ **Use Strategy Pattern (Don't Split Into Separate Classes)**

**Rationale:**
1. **Strategy pattern provides most benefits with least disruption**
   - Same benefits as splitting (isolated mode logic)
   - Avoids complexity of managing multiple chart classes
   - Uses proven pattern (like PotentialRenderer)

2. **Splitting into separate classes would create new problems:**
   - Need to manage 3 chart instances in parent view
   - More complex mode switching logic
   - Potential for drift between implementations
   - More files without clear organizational benefit

3. **Strategy pattern is the established pattern in this codebase**
   - `PotentialRenderer` strategy worked extremely well for `EnergyChartNode`
   - Reduced `updatePotentialCurve()` from 520 lines to ~40 lines
   - Same approach will work for display modes

**The sweet spot**: Keep one unified chart with swappable strategies for mode-specific behavior.

---

## Open Questions for Discussion

1. **Phase 3 vs Phase 4 priority?**
   - Phase 3 (grid lines + axis builder) is quick win, low risk
   - Phase 4 (display strategies) is bigger impact, medium risk
   - Recommendation: Do Phase 3 first, then Phase 4

2. **Should display strategies be classes or objects?**
   ```typescript
   // Option A: Classes (like PotentialRenderer)
   class ProbabilityDensityDisplayStrategy implements WaveFunctionDisplayStrategy { ... }

   // Option B: Strategy objects with functions
   const probabilityDensityStrategy = { updateViewRange: (...) => {...}, render: (...) => {...} }
   ```
   Recommendation: **Classes** (matches PotentialRenderer pattern, more extensible)

3. **Where to put strategy files?**
   - Option A: `src/common/view/wavefunction-strategies/`
   - Option B: `src/common/view/display-strategies/`
   - Recommendation: **Option A** (more specific, follows `potential-renderers/` naming)

4. **Should we extract a ChartAxisBuilder now or wait?**
   - Could be done independently as Phase 3
   - Low risk, medium payoff
   - Recommendation: **Do it** (Phase 3 - it's a quick win)

---

## Next Steps

1. **Review this plan** with team
2. **Get approval** for Phase 3 and/or Phase 4
3. **Create feature branch** from `architecture` branch
4. **Implement Phase 3** (grid lines + axis builder)
5. **Test and review**
6. **Implement Phase 4** (display strategies)
7. **Final testing and merge**

---

## Conclusion

The architecture branch has made **excellent progress** (9% reduction, strategy pattern established). The **biggest remaining opportunity** is applying the same strategy pattern to `WaveFunctionChartNode`'s display modes.

**Key insight**: `WaveFunctionChartNode` doesn't need to be split into separate chart classes. It needs its **display mode logic extracted into strategies**, just like how `EnergyChartNode` extracted potential rendering into strategies.

**Projected total improvement**:
- **~1,200-1,500 line reduction** (32-40% of original size)
- **Main files under 1,000 lines** (all easily maintainable)
- **Excellent testability** (strategies can be unit tested)
- **Clear architecture** (strategy pattern used consistently)

This approach follows the proven **incremental refactoring** methodology established in Phases 1-2. 🎯
