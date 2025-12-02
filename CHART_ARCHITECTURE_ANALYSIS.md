# Chart Architecture Analysis & Recommendations

## Executive Summary

The three chart components (**EnergyChartNode**, **WaveFunctionChartNode**, **WavenumberChartNode**) contain **3,690 lines of code** with significant duplication and inconsistent architecture. This document outlines specific refactoring strategies to reduce complexity and improve maintainability.

---

## Current State Analysis

### File Sizes
- **EnergyChartNode.ts**: 1,464 lines
- **WaveFunctionChartNode.ts**: 1,610 lines (largest)
- **WavenumberChartNode.ts**: 616 lines
- **BaseChartNode.ts**: 246 lines (base class)

### Current Problems

#### 1. **Inconsistent Inheritance Hierarchy** ⚠️ CRITICAL
```
EnergyChartNode extends BaseChartNode ✓
WaveFunctionChartNode extends Node ✗ (should extend BaseChartNode)
WavenumberChartNode extends Node ✗ (should extend BaseChartNode)
```

**Impact**: WaveFunctionChartNode and WavenumberChartNode duplicate all the functionality that BaseChartNode provides:
- Chart dimensions and margins setup
- ChartTransform creation
- Background rectangle
- Clipped plot content node
- Zero line creation
- Coordinate transformation methods (dataToViewX, dataToViewY, viewToDataX, viewToDataY)
- Property initialization

**Estimated Duplication**: ~150-200 lines per chart file

#### 2. **Massive Monolithic Classes**
- EnergyChartNode has **14 major responsibilities**:
  - Axis creation with manual grid lines
  - Potential curve rendering (13 different potential types)
  - Energy level visualization
  - Classical turning points
  - Keyboard navigation
  - PDOM accessibility
  - Legend creation
  - Hit detection for energy levels
  - Hover states
  - Update orchestration

- WaveFunctionChartNode has **10 major responsibilities**:
  - Wavefunction rendering (3 display modes)
  - Superposition handling
  - Time evolution
  - Tool management (6 different tools)
  - RMS statistics calculation and display
  - Zeros visualization
  - Phase color visualization
  - Classical probability overlay

#### 3. **Code Duplication Across Charts**

**Axis Creation** (appears in all 3 charts):
```typescript
// Similar patterns in all three files
const yAxisLeftNode = new AxisLine(this.chartTransform, Orientation.VERTICAL, {...});
const xAxisNode = new AxisLine(this.chartTransform, Orientation.HORIZONTAL, {...});
const yTickMarksNode = new TickMarkSet(...);
const xTickLabelsNode = new TickLabelSet(...);
```

**Grid Lines** (manual implementation in EnergyChartNode and WaveFunctionChartNode):
```typescript
// Duplicated in both files
for (let pos = -X_AXIS_RANGE_NM; pos <= X_AXIS_RANGE_NM; pos += 2) {
  const x = this.chartMargins.left + this.chartTransform.modelToViewX(pos);
  const gridLine = new Line(x, this.chartMargins.top, x, ...);
}
```

**Property Linking Patterns**:
```typescript
// Similar patterns repeated across all charts
this.model.potentialTypeProperty.lazyLink(() => this.update());
this.model.wellWidthProperty.lazyLink(() => this.update());
this.model.wellDepthProperty.lazyLink(() => this.update());
```

#### 4. **Complex Rendering Logic Mixed with View Logic**

**EnergyChartNode.updatePotentialCurve()**: 520 lines (lines 755-1270)
- Contains 13 separate if/else blocks for different potential types
- Each block has custom rendering logic
- Should be extracted into separate strategy classes

**WaveFunctionChartNode plotting methods**: Multiple methods scattered throughout
- plotWaveFunctionComponents()
- plotProbabilityDensityFromArray()
- plotSuperpositionComponents()
- All mixed with chart update logic

---

## Recommended Refactoring Strategy

### Phase 1: Unify Base Class Hierarchy (HIGH PRIORITY)

**Effort**: Low | **Impact**: High | **Risk**: Low

#### Actions:
1. **Extend BaseChartNode consistently**
   - Make WaveFunctionChartNode extend BaseChartNode
   - Make WavenumberChartNode extend BaseChartNode

2. **Remove duplicated code**:
   - Delete duplicate ChartTransform setup
   - Delete duplicate coordinate transformation methods
   - Delete duplicate backgroundRect and plotContentNode setup
   - Delete duplicate property initialization

3. **Enhance BaseChartNode** with common utilities:
   ```typescript
   // Add to BaseChartNode
   protected createStandardAxis(orientation: Orientation): Node {
     // Standard axis creation logic
   }

   protected createGridLines(orientation: Orientation, spacing: number): Node[] {
     // Manual grid line creation (bamboo GridLineSet causes issues)
   }

   protected linkCommonProperties(): void {
     // Common property linking
     this.model.potentialTypeProperty.lazyLink(() => this.update());
     this.model.wellWidthProperty.lazyLink(() => this.update());
     // ...
   }
   ```

**Expected Reduction**: ~300-400 lines total across the three files

---

### Phase 2: Extract Rendering Strategies (MEDIUM PRIORITY)

**Effort**: Medium | **Impact**: High | **Risk**: Medium

#### Problem:
EnergyChartNode.updatePotentialCurve() has 13 different rendering blocks for potential types, making it extremely difficult to:
- Add new potential types
- Modify existing rendering
- Test individual potentials
- Understand the code

#### Solution: Strategy Pattern

Create a **PotentialRenderer** hierarchy:

```typescript
// New file: PotentialRenderer.ts
abstract class PotentialRenderer {
  abstract render(
    xGrid: number[],
    wellWidth: number,
    wellDepth: number,
    viewTransform: {
      dataToViewX: (x: number) => number;
      dataToViewY: (y: number) => number;
    }
  ): Shape;
}

class InfiniteWellRenderer extends PotentialRenderer {
  render(...): Shape {
    // Extract lines 771-794 from EnergyChartNode
  }
}

class FiniteWellRenderer extends PotentialRenderer {
  render(...): Shape {
    // Extract lines 795-808 from EnergyChartNode
  }
}

class HarmonicOscillatorRenderer extends PotentialRenderer {
  render(...): Shape {
    // Extract lines 809-833 from EnergyChartNode
  }
}

// ... 10 more renderers
```

Usage in EnergyChartNode:
```typescript
// Before: 520 lines of if/else
private updatePotentialCurve(boundStates: BoundStateResult): void {
  if (potentialType === PotentialType.INFINITE_WELL) {
    // 100 lines
  } else if (potentialType === PotentialType.FINITE_WELL) {
    // 100 lines
  } // ... 11 more blocks
}

// After: ~20 lines
private readonly renderers = new Map<PotentialType, PotentialRenderer>([
  [PotentialType.INFINITE_WELL, new InfiniteWellRenderer()],
  [PotentialType.FINITE_WELL, new FiniteWellRenderer()],
  // ...
]);

private updatePotentialCurve(boundStates: BoundStateResult): void {
  const renderer = this.renderers.get(this.model.potentialTypeProperty.value);
  if (!renderer) return;

  const shape = renderer.render(
    boundStates.xGrid,
    this.model.wellWidthProperty.value,
    this.model.wellDepthProperty.value,
    { dataToViewX: this.dataToViewX.bind(this), dataToViewY: this.dataToViewY.bind(this) }
  );

  this.potentialPath.shape = shape;
}
```

**File Structure**:
```
src/common/view/
├── chart-renderers/
│   ├── PotentialRenderer.ts (abstract base)
│   ├── InfiniteWellRenderer.ts
│   ├── FiniteWellRenderer.ts
│   ├── HarmonicOscillatorRenderer.ts
│   ├── CoulombRenderer.ts
│   ├── MorseRenderer.ts
│   └── ... (8 more files)
```

**Benefits**:
- Each renderer is ~30-50 lines (testable in isolation)
- Easy to add new potential types (just create new renderer)
- Clear separation of concerns
- EnergyChartNode.ts reduces from 1,464 → ~1,000 lines

**Expected Reduction**: ~450 lines moved to separate renderer files

---

### Phase 3: Extract Axis/Grid Builder (LOW PRIORITY)

**Effort**: Low | **Impact**: Medium | **Risk**: Low

Create a shared **ChartAxisBuilder** utility:

```typescript
// New file: ChartAxisBuilder.ts
export class ChartAxisBuilder {
  constructor(
    private chartTransform: ChartTransform,
    private margins: ChartMargins,
    private colors: typeof QPPWColors
  ) {}

  createVerticalAxis(options: {
    position: 'left' | 'center';
    value?: number;
    opacity?: number;
  }): AxisLine {
    // Shared implementation
  }

  createHorizontalAxis(options: {
    position: 'bottom' | 'center';
    value?: number;
  }): AxisLine {
    // Shared implementation
  }

  createVerticalGridLines(spacing: number, range: Range): Line[] {
    // Manual grid line creation
  }

  createTickMarks(orientation: Orientation, spacing: number): TickMarkSet {
    // Shared tick mark creation
  }

  createTickLabels(
    orientation: Orientation,
    spacing: number,
    formatter: (value: number) => string
  ): TickLabelSet {
    // Shared tick label creation
  }
}
```

Usage:
```typescript
// In EnergyChartNode
protected createAxes(): Node {
  const builder = new ChartAxisBuilder(this.chartTransform, this.chartMargins, QPPWColors);

  const axesNode = new Node();

  // Grid lines
  axesNode.children = [
    ...builder.createVerticalGridLines(2, new Range(-4, 4)),
    ...builder.createHorizontalGridLines(5, new Range(this.yMin, this.yMax))
  ];

  // Axes
  axesNode.addChild(builder.createVerticalAxis({ position: 'left' }));
  axesNode.addChild(builder.createHorizontalAxis({ position: 'bottom' }));

  // Tick marks and labels
  axesNode.addChild(builder.createTickMarks(Orientation.VERTICAL, 5));
  axesNode.addChild(builder.createTickLabels(Orientation.VERTICAL, 5, v => v.toFixed(0)));

  return axesNode;
}
```

**Expected Reduction**: ~100-150 lines per chart

---

### Phase 4: Create Composition-Based Tool System (OPTIONAL)

**Effort**: Medium-High | **Impact**: Medium | **Risk**: Low

WaveFunctionChartNode already has good separation with tool classes, but they're tightly coupled. Consider:

```typescript
// New pattern: Tool registry
class ChartToolRegistry {
  private tools: Map<string, ChartTool> = new Map();

  registerTool(name: string, tool: ChartTool): void {
    this.tools.set(name, tool);
  }

  updateAll(context: ToolUpdateContext): void {
    this.tools.forEach(tool => {
      if (tool.isEnabled()) {
        tool.update(context);
      }
    });
  }
}

// In WaveFunctionChartNode
private toolRegistry = new ChartToolRegistry();

constructor(...) {
  // Register tools
  this.toolRegistry.registerTool('area', this.areaMeasurementTool);
  this.toolRegistry.registerTool('curvature', this.curvatureTool);
  // ...
}

private update(): void {
  // ...
  this.toolRegistry.updateAll({
    displayMode: this.getEffectiveDisplayMode(),
    boundStates,
    selectedIndex
  });
}
```

**Expected Reduction**: ~50-100 lines (mostly organization, not line count)

---

## Implementation Priority

### Immediate (Week 1-2):
1. ✅ **Make WaveFunctionChartNode extend BaseChartNode**
2. ✅ **Make WavenumberChartNode extend BaseChartNode**
3. ✅ **Remove duplicated coordinate transformation code**

**Quick Win**: Reduces total codebase by ~300-400 lines with minimal risk

### Short-term (Week 3-4):
4. ✅ **Extract PotentialRenderer strategy classes**
5. ✅ **Create ChartAxisBuilder utility**

**Major Win**: Reduces complexity significantly, makes code much more maintainable

### Long-term (Future):
6. ⚪ **Refine tool system architecture** (if needed)
7. ⚪ **Extract wavefunction rendering strategies** (if WaveFunctionChartNode is still too large)

---

## Expected Impact

| Metric | Before | After Phase 1 | After Phase 2 | Improvement |
|--------|--------|---------------|---------------|-------------|
| **Total Lines** | 3,690 | 3,400 | 3,000 | -19% |
| **Largest File** | 1,610 | 1,400 | 1,000 | -38% |
| **Code Duplication** | High | Low | Minimal | 🎯 |
| **Testability** | Poor | Good | Excellent | 🎯 |
| **Maintainability** | Poor | Good | Excellent | 🎯 |

---

## Alternative Approaches Considered

### ❌ Template Method Pattern
- **Considered**: Making BaseChartNode more prescriptive with template methods
- **Rejected**: Charts are too different; would create artificial constraints

### ❌ Single Mega Chart Component
- **Considered**: Merging all charts into one with mode switching
- **Rejected**: Would create a 4,000+ line monster; makes things worse

### ❌ Complete Rewrite Using Chart Library
- **Considered**: Using a third-party charting library
- **Rejected**: Quantum mechanics visualizations are too specialized; libraries don't support:
  - Complex wavefunction phase visualization
  - Probability density with quantum turning points
  - Interactive energy level selection
  - Custom potential rendering

---

## Conclusion

The three chart components suffer from **inconsistent inheritance** (critical issue) and **complex monolithic rendering logic** (major issue).

**Recommended approach**:
1. **Phase 1** (Unify base class) is essential and low-risk
2. **Phase 2** (Extract renderers) provides the most value for EnergyChartNode
3. **Phase 3** (Axis builder) is optional polish

This strategy follows **incremental refactoring** principles:
- ✅ Small, safe steps
- ✅ Each phase delivers value
- ✅ No big-bang rewrites
- ✅ Maintains functionality throughout

**Estimated total effort**: 3-4 weeks for Phases 1-2
**Risk level**: Low (changes are isolated and testable)
**ROI**: High (significantly improves maintainability)
