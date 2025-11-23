# Updating Tests to New Energy Convention

When the analytical solver is updated to use the **V=0 convention** (wells at V=0, barrier at V=wellDepth), apply these changes to `test-double-well-comprehensive.ts`:

## Energy Convention Changes

### OLD Convention (Current):

- Wells at V = -V₀ (negative potential)
- Barrier at V = 0 (reference)
- Bound states: **-V₀ < E < 0** (negative energies)
- Ground state has **most negative** energy

### NEW Convention:

- Wells at V = 0 (reference potential)
- Barrier at V = wellDepth (positive value)
- Bound states: **0 < E < wellDepth** (positive energies)
- Ground state has **lowest positive** energy

## Required Code Changes

### 1. Update File Header (lines 5-13)

```typescript
/**
 * Energy Convention:
 * - Wells are at V = 0 (reference potential)
 * - Barrier is at V = wellDepth (positive value)
 * - Bound states have energies: 0 < E < wellDepth
 * - Ground state has lowest positive energy (closest to 0)
 * - Higher energy states are further from 0 (approaching wellDepth)
 */
```

### 2. Update WKB Transmission Function (lines 91-115)

```typescript
/**
 * Calculate WKB transmission coefficient
 * New convention: wells at V=0, barrier at V=wellDepth
 * Energy is measured from well bottom (0 < E < wellDepth for bound states)
 */
function calculateWKBTransmission(
  energy: number, // Energy in eV (measured from well bottom)
  wellDepth: number, // Well depth in eV (barrier height relative to well)
  barrierWidth: number, // Barrier width in nm
  particleMass: number = 1.0,
): number {
  const E = energy * EV_TO_JOULES;
  const V = wellDepth * EV_TO_JOULES;
  const a = barrierWidth * NM_TO_M;
  const m = particleMass * ELECTRON_MASS;

  // For bound states, E < V (energy below barrier top)
  if (E >= V) return 1.0; // Above barrier

  // In barrier region: κ² = 2m(V - E)/ℏ²
  const kappa = Math.sqrt(2 * m * (V - E)) / HBAR;
  const transmission = Math.exp(-2 * kappa * a);

  return transmission;
}
```

### 3. Update Energy Assertions

#### Test 1 (line 214):

```typescript
// OLD:
assert(
  groundState.energy < 0 && groundState.energy > -5.0,
  "Ground state energy should be bound (-wellDepth < E < 0)",
);

// NEW:
assert(
  groundState.energy > 0 && groundState.energy < 5.0,
  "Ground state energy should be bound (0 < E < wellDepth)",
);
```

#### Extreme Tests (line 355):

```typescript
// OLD:
assert(
  groundState.energy < 0 && groundState.energy > -wellDepth,
  `Ground state should be bound (-${wellDepth} < E < 0 eV)`,
);

// NEW:
assert(
  groundState.energy > 0 && groundState.energy < wellDepth,
  `Ground state should be bound (0 < E < ${wellDepth} eV)`,
);
```

#### Parameter Sweep (line 570-572):

```typescript
// OLD:
const allBound = result.states.every(
  (s) => s.energy < 0 && s.energy > -wellDepth,
);

// NEW:
const allBound = result.states.every(
  (s) => s.energy > 0 && s.energy < wellDepth,
);
```

### 4. Update Console Messages (line 209):

```typescript
// OLD:
console.log(
  `  Ground state energy: ${groundState.energy.toFixed(6)} eV (-${5.0} < E < 0)`,
);

// NEW:
console.log(
  `  Ground state energy: ${groundState.energy.toFixed(6)} eV (0 < E < ${5.0} eV)`,
);
```

## Quick Find & Replace Guide

Run these replacements across the file:

1. `energy < 0 && energy > -wellDepth` → `energy > 0 && energy < wellDepth`
2. `energy < 0 && energy > -5.0` → `energy > 0 && energy < 5.0`
3. `(-wellDepth < E < 0` → `(0 < E < wellDepth`
4. `(-${wellDepth} < E < 0` → `(0 < E < ${wellDepth}`
5. `(-${5.0} < E < 0)` → `(0 < E < ${5.0} eV)`

## Verification

After making changes, run:

```bash
npx tsx --import ./tests/browser-globals.js tests/test-double-well-comprehensive.ts
```

All 72 assertions should still pass with the new convention.

## Expected Energy Values

With the new convention, typical ground state energies will be:

- **Old**: E₀ = -4.685 eV (for 5 eV well depth)
- **New**: E₀ = 0.315 eV (for 5 eV well depth)

The difference is exactly the well depth (shift of +5 eV).
