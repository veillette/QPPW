# Quantum Bound States Simulation

An interactive quantum mechanics simulation exploring the behavior of quantum particles in potential wells. This educational tool visualizes wave functions, probability densities, and energy levels for bound states in various quantum potentials.

## Overview

This simulation allows users to explore the properties of quantum "particles" bound in potential wells and observe how the wave functions and probability densities that describe them evolve (or don't evolve) over time. Users can interactively manipulate potential wells and see real-time solutions to the 1D Schrödinger equation.

### Key Features

- **Real-time Quantum Simulation**: Solves the 1D Schrödinger equation numerically in real time with the mass set to the electron mass
- **Interactive Potential Wells**: Click and drag directly on the potential energy diagram to modify well parameters
- **Multiple Well Types**: Includes 9 analytically solvable potentials:
  - Infinite and finite square wells
  - Harmonic oscillator
  - Morse potential (with width parameterization)
  - Pöschl-Teller potential (with width parameterization)
  - Asymmetric triangle potential
  - Coulomb 1D and 3D potentials
  - Double square well (for Two Wells screen)
- **Comprehensive Visualizations**:
  - Wave functions (real and imaginary parts, magnitude)
  - Probability density with filled area visualization
  - Phase-colored visualization showing quantum phase evolution
  - Interactive energy level selection with hover labels
  - High-resolution wavefunction display (1000 points for analytical solutions)
- **Energy Level Display**: Observe discrete energy levels with color-coded selection
- **Time Evolution**: Watch quantum phase evolution in real-time animation with proper eigenstate evolution
- **Superposition States**: Create and visualize quantum superpositions:
  - Single eigenstate (ψₖ)
  - Two-state superposition (ψ₀ + ψ₁)
  - Localized wavepackets (narrow and wide)
  - Coherent states (true coherent states for harmonic oscillator, coherent-like wavepackets for other potentials)
  - Custom superpositions
- **Parameter Controls**: Adjust well width, depth, and offset interactively

### Learning Objectives

This simulation is designed to help students understand:

- **Wave-Particle Duality**: Observe how quantum particles behave as waves
- **Quantum Confinement**: See how boundary conditions lead to discrete energy levels
- **Probability Interpretation**: Understand the relationship between wave function and probability density
- **Energy Splitting**: Use the double well configuration to illustrate symmetric and antisymmetric wave functions
- **Quantum Tunneling**: Visualize tunneling effects in barrier potentials
- **Uncertainty Principle**: Explore the relationship between position and momentum uncertainty
- **Simple Harmonic Oscillator**: Study one of the fundamental quantum systems
- **Coulomb Potential**: Investigate the quantum mechanics of the hydrogen atom

## Topics Covered

- Quantum Mechanics
- Wave Functions
- Energy Levels
- Potential Energy
- Potential Wells
- Probability Density
- Simple Harmonic Oscillator
- Hydrogen Atom
- Superposition
- Uncertainty Principle

## Technical Details

### Requirements

- **Node.js**: Version 18 or higher
- **Package Manager**: npm

### Technology Stack

- **SceneryStack**: Framework for interactive physics simulations
  - Uses Bamboo for chart rendering (ChartTransform, TickMarkSet, TickLabelSet)
  - Scenery for interactive graphics and user input
  - Axon for reactive properties and data binding
- **TypeScript**: Type-safe development
- **Vite**: Fast build tool and development server
- **Quantum Solvers**: Multiple numerical methods with analytical solutions
  - **Analytical Solutions**: Exact solutions for 9 potentials evaluated at 1000 grid points
    - Infinite/finite square wells
    - Harmonic oscillator
    - Morse potential (width-parameterized)
    - Pöschl-Teller potential (width-parameterized)
    - Coulomb 1D and 3D
    - Asymmetric triangle
    - Double square well
  - **Numerical Methods**: For validation and extension
    - DVR (Discrete Variable Representation)
    - Spectral method
    - Matrix Numerov
    - FGH (Fourier Grid Hamiltonian)
- **Time Evolution**: Proper quantum dynamics with individual eigenstate phase evolution for superpositions

### Installation

```bash
# Clone the repository
git clone https://github.com/veillette/QPPW.git

# Navigate to the project directory
cd QPPW

# Install dependencies
npm install
```

### Development

```bash
# Start development server
npm run dev

# Type checking
npm run check

# Lint code
npm run lint

# Format code
npm run format
```

### Building

```bash
# Build for production
npm run build

# Build generates a single HTML file for easy deployment
```

### Testing

The project includes comprehensive test suites to verify the accuracy of quantum mechanical calculations across all numerical methods and potential types.

#### Method 1: Terminal-Based Tests (Fastest)

Run the complete test suite directly in your terminal:

```bash
npm test
```

This command:
- Runs all 45+ accuracy tests across 4 numerical methods (DVR, Spectral, Matrix Numerov, FGH)
- Tests harmonic oscillator, finite square wells, Coulomb potential, and double wells
- Displays detailed pass/fail results with error percentages
- Shows performance metrics (execution times) for each method
- Provides a summary of passed/failed tests

**Output includes:**
- ✓/✗ Pass/fail indicators for each test
- Numerical vs analytical energy comparisons
- Error percentages for each energy level
- Performance summary with timing statistics

#### Method 2: Browser-Based Tests with Visual Interface

The easiest way to run the comprehensive test suite:

```bash
# 1. Build the project (required)
npm run build

# 2. Open the test runner in your browser
# On macOS:
open tests/accuracy-tests.html
# On Linux:
xdg-open tests/accuracy-tests.html
# On Windows:
start tests/accuracy-tests.html
```

The browser test interface provides:
- **Run Full Tests** button - Comprehensive test suite (60+ tests)
- **Run Quick Check** button - Rapid validation subset
- Visual pass/fail indicators with color coding
- Detailed error percentages and performance metrics
- Side-by-side comparison of numerical vs analytical solutions

Results are displayed in the browser console with formatted output.

#### Method 3: Development Server with Interactive Console

Run tests interactively through the development server:

```bash
# 1. Start the development server
npm start

# 2. Open your browser to: http://localhost:5173/tests/accuracy-tests.html
#    (If port 5173 is in use, Vite will use the next available port)

# 3. Click "Run Full Tests" or "Run Quick Check" buttons
```

The test page displays results in a formatted console with color-coded pass/fail indicators.

**Alternative**: Run tests via browser developer console:

```javascript
# From http://localhost:5173 open the console (F12) and run:
import('/src/common/model/AccuracyTests.ts').then(m => m.runAccuracyTests());

# Or for a quick validation:
import('/src/common/model/AccuracyTests.ts').then(m => m.runQuickAccuracyCheck());
```

#### Test Coverage

The comprehensive test suite (`tests/accuracy-tests.html`) validates:

1. **Harmonic Oscillator** - Tests first 10 energy levels with 0.1% error tolerance
2. **Finite Square Wells** (4 configurations) - Tests first 10 energy levels with 0.5% error tolerance
3. **3D Coulomb/Hydrogen Atom** - Tests first 10 energy levels with 1.0% error tolerance
4. **Double Square Wells** (3 configurations) - Tests first 10 energy levels with 1.0% error tolerance

All four numerical methods are tested (DVR, Spectral, Matrix Numerov, FGH) across multiple grid sizes to ensure accuracy and robustness. The tests verify that numerical methods match analytical solutions within specified tolerances.

For more details, see [`tests/README.md`](tests/README.md).

## Deployment

The simulation is designed to be deployed as a single HTML file, making it easy to:

- Host on GitHub Pages
- Embed in educational platforms
- Distribute to students
- Use offline in classrooms

The included GitHub Actions workflow automates:
- Type checking and linting
- Building the simulation
- Deployment to GitHub Pages
- Release creation for tagged versions

## Educational Use

This simulation is ideal for:

- **Introductory Quantum Mechanics Courses**: Visualizing abstract quantum concepts
- **Advanced Physics Classes**: Exploring specific quantum systems in detail
- **Self-Directed Learning**: Interactive exploration of quantum phenomena
- **Classroom Demonstrations**: Real-time visualization of quantum mechanics

### Teaching Tips

1. **Start Simple**: Begin with a single square well to introduce basic concepts
2. **Build Intuition**: Have students manipulate the potential and observe how wave functions change
3. **Compare Systems**: Use different potential configurations to highlight similarities and differences
4. **Quantitative Analysis**: Encourage students to relate observed patterns to theoretical predictions
5. **Double Wells**: Use the two-well configuration to discuss molecular bonding and energy splitting

## Inspiration

This project was inspired by [PhET Interactive Simulations](https://phet.colorado.edu/), specifically the Quantum Bound States simulation, which has been a valuable educational resource for teaching quantum mechanics to students worldwide.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- Built with [SceneryStack](https://github.com/phetsims/scenery) framework
- Inspired by PhET Interactive Simulations at the University of Colorado Boulder
- Developed to make quantum mechanics more accessible and intuitive for students
