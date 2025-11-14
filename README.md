# Quantum Bound States Simulation

An interactive quantum mechanics simulation exploring the behavior of quantum particles in potential wells. This educational tool visualizes wave functions, probability densities, and energy levels for bound states in various quantum potentials.

## Overview

This simulation allows users to explore the properties of quantum "particles" bound in potential wells and observe how the wave functions and probability densities that describe them evolve (or don't evolve) over time. Users can interactively manipulate potential wells and see real-time solutions to the 1D Schrödinger equation.

### Key Features

- **Real-time Quantum Simulation**: Solves the 1D Schrödinger equation numerically in real time with the mass set to the electron mass
- **Interactive Potential Wells**: Click and drag directly on the potential energy diagram to modify well parameters
- **Multiple Well Types**: Includes 10 exactly solvable potentials:
  - Infinite and finite square wells
  - Harmonic oscillator
  - Morse potential
  - Pöschl-Teller potential
  - Rosen-Morse potential
  - Eckart potential
  - Asymmetric triangle potential
  - Coulomb 1D and 3D potentials
- **Comprehensive Visualizations**:
  - Wave functions (real and imaginary parts, magnitude)
  - Probability density with filled area visualization
  - Phase-colored visualization showing quantum phase evolution
  - Interactive energy level selection with hover labels
- **Energy Level Display**: Observe discrete energy levels with color-coded selection
- **Time Evolution**: Watch quantum phase evolution in real-time animation
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
  - Axon for reactive properties
- **TypeScript**: Type-safe development
- **Vite**: Fast build tool and development server
- **Real-time Numerical Solver**: Efficient computation of quantum mechanical solutions
  - DVR (Discrete Variable Representation) method
  - Numerov shooting method
  - Analytical solutions for exactly solvable potentials

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
