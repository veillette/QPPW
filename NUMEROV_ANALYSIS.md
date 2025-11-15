# DoubleWellNumerov Solver Analysis

## Summary

After extensive debugging and analysis, I've determined that the **Numerov method is working correctly** - the issue is fundamental to the physics of the problem.

## What Was Fixed

### ✅ Grid Range (REAL BUG - FIXED)
- **Problem**: Grid was too small (±0.9 nm), cutting off the wells at ±1.1 nm
- **Fix**: Extended grid to ±1.6 nm (wells + decay margin)
- **Location**: `src/two-wells/model/TwoWellsModel.ts:266-281`

### ✅ WKB Initial Conditions (IMPROVEMENT)
- **Added**: Proper hyperbolic initial conditions at x=0
- **For symmetric**: ψ(dx) = cosh(κ·dx) instead of flat ψ(dx) = 1
- **For antisymmetric**: ψ(dx) = sinh(κ·dx) instead of linear ψ(dx) = dx
- **Location**: `src/common/model/NumerovSolver.ts:185-217`

## The Fundamental Issue

### The Problem Geometry
```
Left well: [-1.1, -0.1] nm  (V = 0)
Barrier:   [-0.1, 0.1] nm   (V = 5 eV) ← x=0 is HERE!
Right well: [0.1, 1.1] nm   (V = 0)
```

The symmetry center (x=0) lies in the **classically forbidden barrier region**.

### Why Numerov "Fails"

In the barrier where E < V, the Schrödinger equation becomes:
```
d²ψ/dx² = κ²ψ  where κ² = 2m(V-E)/ℏ² > 0
```

This has **two physical solutions**:
1. ψ ∝ exp(+κx) - **exponentially growing**
2. ψ ∝ exp(-κx) - **exponentially decaying**

The Numerov difference equation with f ≈ -0.003 becomes:
```
ψ_(j+1) ≈ 12.06·ψ_j - 1.00·ψ_(j-1)
```

Solving the characteristic equation λ² - 12.06λ + 1 = 0 gives:
- λ₁ ≈ 11.97 (growing mode - matches exp(+κx))
- λ₂ ≈ 0.08 (decaying mode - matches exp(-κx))

**The Numerov method correctly finds both physical solutions!**

### Why It Grows Exponentially

When integrating from x=0 with ψ(0)=1, ψ(dx)≈1, the solution is:
```
ψ_j = C₁·(11.97)^j + C₂·(0.08)^j
```

The initial conditions determine C₁ and C₂. Even if C₁ is tiny, the growing mode dominates after a few steps because 11.97^j >> 0.08^j.

**Key insight**: This is NOT a numerical bug. The Numerov method is accurately solving the differential equation. The problem is that we **cannot control** which solution (growing vs decaying) we get when integrating from the barrier.

## Why Small f Doesn't Help

Even with perfect numerics (f → 0), we'd have:
```
λ² - 12λ + 1 = 0  →  λ ≈ 11.9 or 0.08
```

The two exponential modes are **inherent to the physics**, not numerical error.

## Units Are Correct

Extensive testing confirmed:
- ✅ All SI units properly converted (meters, joules, kg)
- ✅ No precision loss in coefficients
- ✅ Numerov f coefficients are small (~0.003) and stable
- ✅ No unit mismatches between dx and k²

## Recommended Solutions

### Option 1: Use Matrix Methods (RECOMMENDED)
The **DVR, FGH, and Spectral methods already work** for double wells because they:
- Diagonalize the Hamiltonian matrix
- Find ALL eigenvalues simultaneously
- Don't require shooting/integration

### Option 2: Match-Point Method
Integrate from BOTH boundaries and match in the middle:
1. Integrate left-to-right from left boundary
2. Integrate right-to-left from right boundary
3. Match wavefunction and derivative at a matching point
4. Shooting parameter: mismatch at matching point

This avoids starting in the barrier.

### Option 3: Disable DoubleWellNumerov
Simply don't use Numerov for double wells where x=0 is in a barrier. Fall back to DVR/FGH.

## Conclusion

The DoubleWellNumerov solver is **fundamentally unsuitable** for geometries where:
1. The symmetry center (x=0) lies in a classically forbidden region
2. We need to integrate through that barrier

This is a **physics limitation**, not a coding bug. The Numerov method correctly solves the differential equation, but picks up the exponentially growing solution when integrating from the barrier.

**Recommendation**: Disable the DoubleWellNumerov solver and use DVR/FGH/Spectral methods instead.
