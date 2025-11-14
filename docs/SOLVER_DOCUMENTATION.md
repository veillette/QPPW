# 1D Bound State Solver Documentation

## Overview

This document describes the implementation of the 1D time-independent Schrödinger equation (TISE) solver for the QPPW quantum physics simulation.

## Features

### Analytical Solutions
For well-known potentials, the solver provides exact analytical solutions:
- **Infinite Square Well**: $E_n = \frac{n^2 \pi^2 \hbar^2}{2mL^2}$
- **Harmonic Oscillator**: $E_n = \hbar\omega(n + \frac{1}{2})$

### Numerical Solutions
For arbitrary potentials, two numerical methods are available:

#### 1. Numerov Method
- Higher-order finite-difference method with $O(h^6)$ error
- Uses shooting method to find bound states
- Iterative formula:
  $$\psi_{j+1} = \frac{(12 - 10f_j)\psi_j - (1+f_{j-1})\psi_{j-1}}{1+f_{j+1}}$$
  where $f_j = \frac{h^2}{12}k^2(x_j)$ and $k^2(x) = \frac{2m(E-V(x))}{\hbar^2}$

#### 2. Discrete Variable Representation (DVR) Method
- Matrix diagonalization approach
- Constructs Hamiltonian $H = T + V$
- Potential energy: Diagonal matrix with $V_{ii} = V(x_i)$
- Kinetic energy: Colbert-Miller formula
  $$T_{ij} = \frac{\hbar^2}{2m\Delta x^2} \begin{cases} \frac{\pi^2}{3} & \text{for } i=j \\ \frac{2(-1)^{i-j}}{(i-j)^2} & \text{for } i \neq j \end{cases}$$

## File Structure

```
src/common/model/
├── QuantumConstants.ts          # Physical constants (ℏ, m_e, etc.)
├── PotentialFunction.ts         # Type definitions and interfaces
├── AnalyticalSolutions.ts       # Analytical solvers
├── NumerovSolver.ts            # Numerov method implementation
├── DVRSolver.ts                # DVR method implementation
├── Schrodinger1DSolver.ts      # Main solver class
└── SolverExamples.ts           # Usage examples and tests
```

## Usage

### Basic Usage

```typescript
import Schrodinger1DSolver, { NumericalMethod } from './common/model/Schrodinger1DSolver.js';
import { PotentialType } from './common/model/PotentialFunction.js';
import QuantumConstants from './common/model/QuantumConstants.js';

// Create solver instance
const solver = new Schrodinger1DSolver(NumericalMethod.DVR);

// Define grid
const gridConfig = {
  xMin: 0,
  xMax: 1e-9,  // 1 nm
  numPoints: 200
};

// Solve for infinite well (analytical)
const result = solver.solveAnalyticalIfPossible(
  {
    type: PotentialType.INFINITE_WELL,
    wellWidth: 1e-9
  },
  QuantumConstants.ELECTRON_MASS,
  5,  // Number of states
  gridConfig
);

// Access results
result.energies.forEach((E, n) => {
  const E_eV = Schrodinger1DSolver.joulesToEV(E);
  console.log(`E_${n+1} = ${E_eV.toFixed(4)} eV`);
});
```

### Custom Potentials

```typescript
// Define custom potential function
const potential = (x: number) => {
  // Finite square well
  if (x >= 0 && x <= 1e-9) {
    return -Schrodinger1DSolver.eVToJoules(5);  // -5 eV
  }
  return 0;
};

// Solve numerically
const result = solver.solveNumerical(
  potential,
  QuantumConstants.ELECTRON_MASS,
  5,
  gridConfig
);
```

### Helper Functions

The `Schrodinger1DSolver` class provides static helper methods:

```typescript
// Create common potential functions
const infiniteWell = Schrodinger1DSolver.createInfiniteWellPotential(1e-9);
const harmonicOsc = Schrodinger1DSolver.createHarmonicOscillatorPotential(1000);
const finiteWell = Schrodinger1DSolver.createFiniteWellPotential(1e-9, 5e-19);

// Unit conversions
const energyJoules = Schrodinger1DSolver.eVToJoules(1.0);
const energyEV = Schrodinger1DSolver.joulesToEV(1.602e-19);
```

## User Preferences

Users can select the numerical method through the preferences system:

```typescript
import QPPWPreferences from './QPPWPreferences.js';
import { NumericalMethod } from './common/model/Schrodinger1DSolver.js';

// Get current method
const method = QPPWPreferences.numericalMethodProperty.value;

// Set method
QPPWPreferences.numericalMethodProperty.value = NumericalMethod.NUMEROV;
```

## Integration with Models

The solver is integrated into the `OneWellModel` class:

```typescript
// Get energy for quantum number n
const energy = oneWellModel.getEnergyLevel(3);  // E_3 in eV

// Get wavefunction
const wavefunction = oneWellModel.getWavefunction(3);  // ψ_3(x)

// Get spatial grid
const xGrid = oneWellModel.getXGrid();  // x values in nm

// Get all bound states
const boundStates = oneWellModel.getBoundStates();
```

## Running Examples

To test the solver implementation:

```typescript
import { runAllExamples } from './common/model/SolverExamples.js';

// Run all examples in console
runAllExamples();
```

This will execute examples demonstrating:
1. Infinite well with analytical solution
2. Infinite well with DVR method
3. Harmonic oscillator with analytical solution
4. Finite well with DVR method
5. Comparison of Numerov vs DVR methods

## Performance Considerations

### Grid Resolution
- More grid points → higher accuracy, slower computation
- Typical range: 100-500 points
- DVR method scales as O(N³) due to matrix diagonalization
- Numerov method scales as O(N) per energy search

### Method Selection
- **DVR**: Better for general potentials, finds all states simultaneously
- **Numerov**: Better for deep wells, requires energy range specification

### Caching
The `OneWellModel` caches results and recalculates only when parameters change.

## Mathematical Background

### Time-Independent Schrödinger Equation
$$-\frac{\hbar^2}{2m}\frac{d^2\psi}{dx^2} + V(x)\psi = E\psi$$

### Bound State Conditions
1. $\psi(x) \to 0$ as $x \to \pm\infty$
2. $\psi$ and $d\psi/dx$ continuous
3. $\int|\psi|^2 dx = 1$ (normalization)

### DVR Basis
The DVR method uses a basis of delta functions on grid points, making the potential matrix diagonal while using the Colbert-Miller formula for exact kinetic energy in the DVR basis.

## References

1. **Numerov Method**:
   - Numerov, B. (1924). "Note on the numerical integration of d²x/dt² = f(x,t)"

2. **DVR Method**:
   - Colbert, D. T., & Miller, W. H. (1992). "A novel discrete variable representation for quantum mechanical reactive scattering via the S-matrix Kohn method." *Journal of Chemical Physics*, 96(3), 1982-1991.

3. **Quantum Mechanics**:
   - Griffiths, D. J., & Schroeter, D. F. (2018). *Introduction to Quantum Mechanics* (3rd ed.). Cambridge University Press.

## Visualization Features

The solver is integrated with comprehensive visualization tools:

### Energy Chart
- Displays potential energy curves for all potential types
- Shows discrete energy levels as horizontal lines
- Interactive selection of energy levels
- Hover to display energy values
- Color-coded visualization with phase-dependent coloring

### Wavefunction Chart
The simulation provides three visualization modes:
1. **Probability Density** (`|ψ|²`): Shows the probability distribution
2. **Wavefunction Components**: Real part, imaginary part, and magnitude
3. **Phase Color**: Rainbow-colored visualization where hue represents quantum phase

### Chart Specifications
- **Energy Chart**: 600×300 pixels with margins (left: 60, right: 20, top: 40, bottom: 50)
- **Wavefunction Chart**: 600×140 pixels with margins (left: 60, right: 20, top: 10, bottom: 40)
- Both charts share synchronized x-axis (Position in nm: -4 to +4 nm)
- Y-axis ranges automatically adjust based on potential type and selected state

## Future Enhancements

Potential improvements to consider:
- Support for 2D/3D potentials
- Time-dependent Schrödinger equation solver (partial implementation exists for phase evolution)
- More analytical solutions (currently supports 10 potential types)
- Adaptive grid refinement
- GPU acceleration for large grids
