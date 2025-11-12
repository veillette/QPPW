# Contributing to Quantum Bound States Simulation

Thank you for your interest in contributing to the Quantum Bound States Simulation project! This document provides guidelines for contributing to the project.

## Getting Started

1. **Fork the repository** on GitHub
2. **Clone your fork** locally:
   ```bash
   git clone https://github.com/YOUR-USERNAME/QPPW.git
   cd QPPW
   ```
3. **Install dependencies**:
   ```bash
   npm install
   ```
4. **Create a branch** for your changes:
   ```bash
   git checkout -b feature/your-feature-name
   ```

## Development Workflow

### Running the Development Server

```bash
npm run dev
```

This will start a local development server with hot module reloading.

### Code Quality

Before submitting your changes, ensure your code passes all quality checks:

```bash
# Type checking
npm run check

# Linting
npm run lint

# Code formatting
npm run format
```

### Building

Test that your changes build successfully:

```bash
npm run build
```

## Contribution Guidelines

### Code Style

- Follow the existing code style and conventions
- Use TypeScript for type safety
- Write clear, descriptive variable and function names
- Add comments for complex logic or physics concepts
- Use Prettier for code formatting (configured in `.prettierrc`)
- Follow ESLint rules (configured in `eslint.config.js`)

### Physics Accuracy

- Ensure quantum mechanics implementations are physically accurate
- Reference reliable sources for physics equations and concepts
- Test numerical solvers for stability and accuracy
- Document any approximations or simplifications made

### Commit Messages

Write clear, descriptive commit messages:

- Use the present tense ("Add feature" not "Added feature")
- Use the imperative mood ("Move cursor to..." not "Moves cursor to...")
- Limit the first line to 72 characters or less
- Reference issues and pull requests liberally

Examples:
```
Add harmonic oscillator potential visualization
Fix wave function normalization bug
Update README with installation instructions
```

### Pull Request Process

1. **Update documentation** if you're changing functionality
2. **Add or update tests** if applicable
3. **Ensure all checks pass** (type checking, linting, building)
4. **Write a clear PR description**:
   - Describe what changes you made and why
   - Reference any related issues
   - Include screenshots or GIFs for visual changes
   - List any breaking changes

5. **Request a review** from a maintainer
6. **Address feedback** promptly and courteously

## Types of Contributions

### Bug Reports

When reporting bugs, please include:

- A clear, descriptive title
- Steps to reproduce the issue
- Expected behavior vs. actual behavior
- Screenshots or error messages if applicable
- Your browser and operating system
- Any relevant console errors

### Feature Requests

When suggesting features, please:

- Explain the educational value of the feature
- Describe how it would work from a user's perspective
- Consider if it aligns with the project's goals
- Provide examples or mockups if possible

### Code Contributions

Areas where contributions are especially welcome:

- **New Potential Types**: Implementing additional quantum well configurations
- **Visualization Enhancements**: Improving graphics and animations
- **Educational Features**: Adding teaching tools or explanations
- **Performance Improvements**: Optimizing numerical solvers
- **Accessibility**: Making the simulation more accessible
- **Documentation**: Improving README, comments, or adding tutorials
- **Testing**: Adding unit tests or integration tests
- **Bug Fixes**: Addressing issues from the issue tracker

### Documentation

Help improve documentation by:

- Fixing typos or clarifying explanations
- Adding examples or tutorials
- Improving code comments
- Creating educational materials for teachers

## Physics and Mathematics

When working on quantum mechanics implementations:

- **Use SI units** consistently (or clearly document unit systems)
- **Normalize wave functions** properly
- **Handle boundary conditions** correctly
- **Validate against known solutions** (e.g., particle in a box, harmonic oscillator)
- **Reference equations** from standard quantum mechanics textbooks
- **Consider numerical stability** for edge cases

## Code of Conduct

### Our Pledge

We are committed to providing a welcoming and inclusive environment for all contributors, regardless of experience level, background, or identity.

### Expected Behavior

- Be respectful and considerate
- Welcome newcomers and help them learn
- Accept constructive criticism gracefully
- Focus on what's best for the project and the community
- Show empathy towards others

### Unacceptable Behavior

- Harassment, discrimination, or offensive comments
- Personal attacks or insults
- Trolling or inflammatory comments
- Spam or off-topic content

## Questions?

If you have questions about contributing:

- Check existing issues and pull requests
- Open a new issue with the "question" label
- Reach out to the maintainers

## License

By contributing to this project, you agree that your contributions will be licensed under the MIT License.

## Recognition

Contributors will be recognized in:

- The project's README (for significant contributions)
- Release notes (for features and bug fixes)
- The GitHub contributors page

Thank you for helping make quantum mechanics education more accessible and engaging!
