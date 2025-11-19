/**
 * Advanced Quantum Bound State Solver for the Time-Independent Schrodinger Equation
 *
 * Features:
 * - Shooting method with adaptive energy bracketing
 * - 6th-order accurate Numerov integration
 * - Automatic matching point optimization
 * - Richardson extrapolation for improved accuracy
 * - Cached computations for performance
 *
 * The TISE is: -hbar^2/(2m) d^2psi/dx^2 + V(x)psi = E*psi
 * Rearranged as: d^2psi/dx^2 = -k^2(x)psi where k^2(x) = 2m(E - V(x))/hbar^2
 */

import QuantumConstants from "./QuantumConstants.js";
import { BoundStateResult, GridConfig, PotentialFunction } from "./PotentialFunction.js";
import qppw from "../../QPPWNamespace.js";

// ============================================================================
// Constants
// ============================================================================

// Numerical method constants
const DEFAULT_MAX_ITERATIONS = 100;
const DEFAULT_CONVERGENCE_TOLERANCE = 1e-10;
const DEFAULT_ENERGY_TOLERANCE = 1e-12;
const DEFAULT_WAVE_FUNCTION_TOLERANCE = 1e-8;

// Integration constants
const NUMEROV_COEFFICIENT = 1.0 / 12.0;  // h^2/12 factor in Numerov method
const MATCHING_POINT_FRACTION = 0.49;    // Matching point fraction for shooting method
const INITIAL_ENERGY_SCALING = 10.0;     // Safety factor for energy bounds
const ENERGY_BRACKET_MULTIPLIER = 2.0;   // Factor for expanding energy search

// Boundary condition constants
const BOUNDARY_PSI_INITIAL = 0.0;        // psi(boundary) = 0
const BOUNDARY_PSI_DERIVATIVE = 1.0;     // Initial derivative (will be normalized)

// Performance constants
const CACHE_SIZE_LIMIT = 100;

// ============================================================================
// Types and Interfaces
// ============================================================================

/**
 * Solver configuration options
 */
export type SolverConfig = {
  reportWarnings?: boolean;
  maxIterations?: number;
  convergenceTolerance?: number;
  energyTolerance?: number;
  waveFunctionTolerance?: number;
  adaptiveMatching?: boolean;
  useRichardsonExtrapolation?: boolean;
  normalizationMethod?: 'max' | 'l2';
  enableCaching?: boolean;
};

/**
 * Detailed solution information
 */
export type QuantumState = {
  energy: number;
  waveFunction: number[];
  nodes: number;
  normalizationConstant: number;
  matchingPoint: number;
  convergenceMetric: number;
};

/**
 * Integration result with detailed metrics
 */
class IntegrationResult {
  constructor(
    public readonly nodeCount: number,
    public readonly logDerivativeMismatch: number,
    public readonly matchingAmplitude: number = 0,
    public readonly convergenceQuality: number = 0
  ) {}

  hasTooManyNodes(desiredNodes: number): boolean {
    return this.nodeCount > desiredNodes ||
           (this.nodeCount === desiredNodes && this.logDerivativeMismatch < 0.0);
  }

  get isConverged(): boolean {
    return Math.abs(this.logDerivativeMismatch) < DEFAULT_CONVERGENCE_TOLERANCE;
  }
}

/**
 * Cache entry for memoization
 */
type CacheEntry<T> = {
  key: string;
  value: T;
  timestamp: number;
};

/**
 * Custom exception for solver errors
 */
export class QuantumSolverException extends Error {
  constructor(
    message: string,
    public readonly code: string = 'QUANTUM_SOLVER_ERROR'
  ) {
    super(message);
    this.name = 'QuantumSolverException';
    Object.setPrototypeOf(this, QuantumSolverException.prototype);
  }
}

// ============================================================================
// Main Solver Class
// ============================================================================

export class QuantumBoundStateSolver {
  // Configuration
  private readonly maxIterations: number;
  private readonly convergenceTolerance: number;
  private readonly energyTolerance: number;
  private readonly waveFunctionTolerance: number;
  private readonly adaptiveMatching: boolean;
  private readonly useRichardsonExtrapolation: boolean;
  private readonly normalizationMethod: 'max' | 'l2';
  private readonly enableCaching: boolean;

  // Grid properties
  private readonly deltaX: number;
  private readonly gridPositions: number[];
  private readonly potentialEnergies: number[];
  private readonly hbarSquaredOver2m: number;
  private gridPoints: number;

  // Cache
  private readonly cache: Map<string, CacheEntry<unknown>>;
  private readonly numerovFactorCache: Map<number, number>;

  // Performance metrics
  private totalIterations = 0;
  private cacheHits = 0;
  private cacheMisses = 0;

  /**
   * Constructor for the Enhanced Quantum Bound State Solver
   *
   * @param mass - Particle mass in kg
   * @param xMin - Minimum x value in meters
   * @param xMax - Maximum x value in meters
   * @param numPoints - Number of grid points
   * @param potentialFunction - Function V(x) returning potential in Joules
   * @param config - Optional solver configuration
   */
  constructor(
    mass: number,
    xMin: number,
    xMax: number,
    numPoints: number,
    private readonly potentialFunction: PotentialFunction,
    config?: SolverConfig
  ) {
    const { HBAR } = QuantumConstants;
    this.hbarSquaredOver2m = (HBAR * HBAR) / (2 * mass);
    this.gridPoints = numPoints;

    // Validate inputs
    this.validateInputs(xMin, xMax, numPoints);

    // Set configuration with defaults
    this.maxIterations = config?.maxIterations ?? DEFAULT_MAX_ITERATIONS;
    this.convergenceTolerance = config?.convergenceTolerance ?? DEFAULT_CONVERGENCE_TOLERANCE;
    this.energyTolerance = config?.energyTolerance ?? DEFAULT_ENERGY_TOLERANCE;
    this.waveFunctionTolerance = config?.waveFunctionTolerance ?? DEFAULT_WAVE_FUNCTION_TOLERANCE;
    this.adaptiveMatching = config?.adaptiveMatching ?? true;
    this.useRichardsonExtrapolation = config?.useRichardsonExtrapolation ?? false;
    this.normalizationMethod = config?.normalizationMethod ?? 'l2';
    this.enableCaching = config?.enableCaching ?? true;

    // Pre-calculate grid properties
    this.deltaX = (xMax - xMin) / (numPoints - 1);
    this.gridPositions = this.calculateGridPositions(xMin, xMax, numPoints);
    this.potentialEnergies = this.calculatePotentialOnGrid();

    // Initialize caches
    this.cache = new Map();
    this.numerovFactorCache = new Map();
  }

  // ========================================================================
  // Public API
  // ========================================================================

  /**
   * Finds the nth eigenstate (energy and wavefunction)
   */
  public findEigenstate(n: number): QuantumState {
    const desiredNodes = n - 1;  // nth state has n-1 nodes

    // Check cache first
    const cacheKey = `eigenstate_${n}`;
    if (this.enableCaching) {
      const cached = this.getFromCache<QuantumState>(cacheKey);
      if (cached) return cached;
    }

    // Find energy
    const energy = this.findEnergyEigenvalue(desiredNodes);

    // Calculate wavefunction
    const waveFunction = this.calculateNormalizedWaveFunction(energy);

    // Determine optimal matching point
    const matchingPoint = this.adaptiveMatching
      ? this.findOptimalMatchingPoint(energy, desiredNodes)
      : Math.floor(this.gridPoints * MATCHING_POINT_FRACTION);

    // Calculate normalization constant
    const normalizationConstant = this.calculateNormalizationConstant(waveFunction);

    // Verify solution quality
    const convergenceMetric = this.calculateConvergenceMetric(energy, waveFunction, desiredNodes);

    const state: QuantumState = {
      energy,
      waveFunction,
      nodes: desiredNodes,
      normalizationConstant,
      matchingPoint,
      convergenceMetric
    };

    // Cache the result
    if (this.enableCaching) {
      this.addToCache(cacheKey, state);
    }

    return state;
  }

  /**
   * Finds multiple eigenstates efficiently
   */
  public findMultipleEigenstates(nStates: number): QuantumState[] {
    const states: QuantumState[] = [];

    for (let n = 1; n <= nStates; n++) {
      states.push(this.findEigenstate(n));
    }

    return states;
  }

  /**
   * Returns the grid positions
   */
  public getGridPositions(): number[] {
    return this.gridPositions.slice();
  }

  /**
   * Gets performance metrics
   */
  public getPerformanceMetrics() {
    return {
      totalIterations: this.totalIterations,
      cacheHits: this.cacheHits,
      cacheMisses: this.cacheMisses,
      cacheHitRate: this.cacheHits / Math.max(1, this.cacheHits + this.cacheMisses),
      averageIterationsPerState: this.totalIterations / Math.max(1, this.cache.size)
    };
  }

  /**
   * Clears all caches
   */
  public clearCache(): void {
    this.cache.clear();
    this.numerovFactorCache.clear();
    this.cacheHits = 0;
    this.cacheMisses = 0;
  }

  // ========================================================================
  // Core Algorithm Methods
  // ========================================================================

  /**
   * Enhanced energy eigenvalue finder with adaptive bracketing
   */
  private findEnergyEigenvalue(desiredNodes: number): number {
    // Get initial energy bounds
    const { lowerBound, upperBound } = this.findEnergyBounds(desiredNodes);

    // Binary search with Richardson extrapolation if enabled
    let energy: number;
    if (this.useRichardsonExtrapolation) {
      energy = this.findEnergyWithRichardson(lowerBound, upperBound, desiredNodes);
    } else {
      energy = this.findEnergyWithBinarySearch(lowerBound, upperBound, desiredNodes);
    }

    // Fine-tune with secant method
    energy = this.refineEnergyWithSecant(energy, desiredNodes);

    return energy;
  }

  /**
   * Finds energy bounds using adaptive bracketing
   */
  private findEnergyBounds(desiredNodes: number): { lowerBound: number; upperBound: number } {
    const boxLength = this.gridPositions[this.gridPoints - 1] - this.gridPositions[0];

    // Get potential range to determine appropriate energy bounds
    const minPotential = Math.min(...this.potentialEnergies);
    const maxPotential = Math.max(...this.potentialEnergies);

    // Particle in box estimate (for positive energy states)
    const particleInBoxEstimate = this.hbarSquaredOver2m * Math.PI * Math.PI *
                                  (desiredNodes + 1) * (desiredNodes + 1) /
                                  (boxLength * boxLength);

    // Determine if we're looking for negative energy states (finite well type)
    // For a finite well with V = -V0 inside and V = 0 outside,
    // bound states have energies between minPotential and maxPotential
    const hasNegativePotential = minPotential < 0;

    let upperBound: number;
    let lowerBound: number;

    if (hasNegativePotential) {
      // For finite wells: search between potential minimum and maximum
      // Bound states have energies between minPotential and maxPotential
      upperBound = maxPotential + particleInBoxEstimate * 0.1;  // Slightly above the well edge
      lowerBound = minPotential * 1.1;  // Slightly below the well bottom
    } else {
      // For positive potentials (like harmonic oscillator): use standard estimate
      upperBound = particleInBoxEstimate * INITIAL_ENERGY_SCALING;
      lowerBound = -Math.abs(upperBound);
    }

    // Refine upper bound (find energy with too many nodes)
    for (let i = 0; i < this.maxIterations; i++) {
      const upperTest = this.integrateSchrodinger(upperBound, desiredNodes);
      if (upperTest.hasTooManyNodes(desiredNodes)) break;
      if (hasNegativePotential) {
        // For negative potential wells, increase energy toward 0
        upperBound += particleInBoxEstimate * 0.5;
        if (upperBound > maxPotential + particleInBoxEstimate) break;
      } else {
        upperBound *= ENERGY_BRACKET_MULTIPLIER;
      }
    }

    // Refine lower bound (find energy with too few nodes)
    for (let i = 0; i < this.maxIterations; i++) {
      const lowerTest = this.integrateSchrodinger(lowerBound, desiredNodes);
      if (!lowerTest.hasTooManyNodes(desiredNodes)) break;
      if (hasNegativePotential) {
        // For negative potential wells, decrease energy further below well bottom
        lowerBound *= 1.2;
      } else {
        lowerBound *= ENERGY_BRACKET_MULTIPLIER;
      }
    }

    return { lowerBound, upperBound };
  }

  /**
   * Binary search for energy
   */
  private findEnergyWithBinarySearch(
    lowerBound: number,
    upperBound: number,
    desiredNodes: number
  ): number {
    let lower = lowerBound;
    let upper = upperBound;
    let energy = 0.5 * (lower + upper); // Initialize to midpoint

    for (let i = 0; i < this.maxIterations; i++) {
      energy = 0.5 * (lower + upper);
      const test = this.integrateSchrodinger(energy, desiredNodes);

      // Check for invalid integration results
      if (!Number.isFinite(test.logDerivativeMismatch)) {
        // Try a different energy if integration fails
        if (test.nodeCount > desiredNodes) {
          upper = energy;
        } else {
          lower = energy;
        }
        continue;
      }

      if (test.isConverged) break;

      if (test.hasTooManyNodes(desiredNodes)) {
        upper = energy;
      } else {
        lower = energy;
      }

      if (Math.abs(upper - lower) < this.energyTolerance) break;

      this.totalIterations++;
    }

    return energy;
  }

  /**
   * Richardson extrapolation for improved accuracy
   */
  private findEnergyWithRichardson(
    lowerBound: number,
    upperBound: number,
    desiredNodes: number
  ): number {
    // Get energy at current grid resolution
    const energy1 = this.findEnergyWithBinarySearch(lowerBound, upperBound, desiredNodes);

    // Temporarily double grid points
    const originalGridPoints = this.gridPoints;
    this.gridPoints = originalGridPoints * 2;
    const energy2 = this.findEnergyWithBinarySearch(lowerBound, upperBound, desiredNodes);
    this.gridPoints = originalGridPoints;

    // Richardson extrapolation formula
    return energy2 + (energy2 - energy1) / 3.0;
  }

  /**
   * Secant method for final energy refinement.
   *
   * The secant method finds roots of f(E) = logDerivativeMismatch by iteratively
   * approximating the derivative using two previous points:
   *   E_new = E_current - f(E_current) * (E_current - E_previous) / (f(E_current) - f(E_previous))
   *
   * This converges faster than bisection (superlinear) while avoiding the need
   * for analytical derivatives.
   */
  private refineEnergyWithSecant(initialEnergy: number, desiredNodes: number): number {
    // Initialize with two energy guesses bracketing the solution
    // These are slightly perturbed from the initial estimate
    let energyPrevious = initialEnergy * 0.999;  // Previous energy iterate
    let energyCurrent = initialEnergy * 1.001;   // Current energy iterate

    // Evaluate the mismatch function at both initial points
    // f(E) = log derivative mismatch, which should be zero at the true eigenvalue
    let mismatchPrevious = this.integrateSchrodinger(energyPrevious, desiredNodes).logDerivativeMismatch;
    let mismatchCurrent = this.integrateSchrodinger(energyCurrent, desiredNodes).logDerivativeMismatch;

    // Secant iteration: refine energy until convergence
    for (let iteration = 0; iteration < 10; iteration++) {
      // Check if the denominator (slope approximation) is too small
      // This would cause numerical instability in the secant formula
      const mismatchDifference = mismatchCurrent - mismatchPrevious;
      if (Math.abs(mismatchDifference) < this.convergenceTolerance) break;

      // Apply the secant formula to compute the next energy estimate
      // This approximates Newton's method using finite differences
      const energyNext = energyCurrent - mismatchCurrent * (energyCurrent - energyPrevious) / mismatchDifference;

      // Evaluate the mismatch at the new energy
      const mismatchNext = this.integrateSchrodinger(energyNext, desiredNodes).logDerivativeMismatch;

      // Check for convergence: mismatch is essentially zero
      if (Math.abs(mismatchNext) < this.convergenceTolerance) return energyNext;

      // Shift the iteration window: current becomes previous, next becomes current
      // This maintains the two-point history needed for the secant approximation
      energyPrevious = energyCurrent;
      mismatchPrevious = mismatchCurrent;
      energyCurrent = energyNext;
      mismatchCurrent = mismatchNext;
    }

    return energyCurrent;
  }

  /**
   * Enhanced Schrodinger equation integration with optimizations
   */
  private integrateSchrodinger(energy: number, desiredNodes: number): IntegrationResult {
    const matchPoint = this.adaptiveMatching
      ? this.findOptimalMatchingPoint(energy, desiredNodes)
      : Math.floor(this.gridPoints * MATCHING_POINT_FRACTION);

    const numerovFactor = this.getNumerovFactor();
    const inverseHbar2Over2m = 1.0 / this.hbarSquaredOver2m;

    // Pre-calculate effective potential
    const effectivePotential = this.potentialEnergies.map(V =>
      inverseHbar2Over2m * (V - energy)
    );

    // Forward integration
    const forwardResult = this.integrateForward(
      effectivePotential,
      numerovFactor,
      matchPoint
    );

    // Backward integration
    const backwardResult = this.integrateBackward(
      effectivePotential,
      numerovFactor,
      matchPoint
    );

    // Calculate total nodes
    const nodeCount = forwardResult.nodes + backwardResult.nodes;

    // Calculate log derivative mismatch
    const logDerivativeMismatch =
      forwardResult.logDerivative + backwardResult.logDerivative;

    // Calculate matching amplitude for quality assessment
    const matchingAmplitude = Math.abs(forwardResult.psiAtMatch);

    // Calculate convergence quality metric
    const convergenceQuality = Math.exp(-Math.abs(logDerivativeMismatch) / this.convergenceTolerance);

    return new IntegrationResult(
      nodeCount,
      logDerivativeMismatch,
      matchingAmplitude,
      convergenceQuality
    );
  }

  /**
   * Forward integration using Numerov method with periodic renormalization
   */
  private integrateForward(
    effectivePotential: number[],
    numerovFactor: number,
    matchPoint: number
  ): { nodes: number; logDerivative: number; psiAtMatch: number } {
    let psi_prev = BOUNDARY_PSI_INITIAL;
    let psi_current = BOUNDARY_PSI_INITIAL;
    let psi_next = BOUNDARY_PSI_DERIVATIVE;

    let nodes = 0;
    const psiValues: number[] = [psi_current, psi_next];
    let scaleFactor = 1.0;

    // Renormalization threshold to prevent overflow
    const RENORM_THRESHOLD = 1e10;

    for (let i = 2; i <= matchPoint + 1; i++) {
      psi_prev = psi_current;
      psi_current = psi_next;

      // Numerov recursion - indices must match wavefunction indices
      // psi_prev = psi[i-2], psi_current = psi[i-1], psi_next = psi[i]
      psi_next = this.numerovStep(
        psi_current,
        psi_prev,
        effectivePotential[i-2],
        effectivePotential[i-1],
        effectivePotential[i],
        numerovFactor
      );

      // Periodic renormalization to prevent overflow
      if (Math.abs(psi_next) > RENORM_THRESHOLD) {
        const norm = Math.abs(psi_next);
        psi_prev /= norm;
        psi_current /= norm;
        psi_next /= norm;
        scaleFactor *= norm;
        // Renormalize stored values
        for (let j = 0; j < psiValues.length; j++) {
          psiValues[j] /= norm;
        }
      }

      psiValues.push(psi_next);

      // Count nodes
      if (i <= matchPoint && this.isNode(psi_current, psi_next)) {
        nodes++;
      }
    }

    // Calculate log derivative using three-point formula
    const logDerivative = this.calculateLogDerivative(
      psiValues[matchPoint - 1],
      psiValues[matchPoint],
      psiValues[matchPoint + 1]
    );

    return {
      nodes,
      logDerivative,
      psiAtMatch: psiValues[matchPoint] * scaleFactor
    };
  }

  /**
   * Backward integration using Numerov method with periodic renormalization
   */
  private integrateBackward(
    effectivePotential: number[],
    numerovFactor: number,
    matchPoint: number
  ): { nodes: number; logDerivative: number; psiAtMatch: number } {
    let psi_prev = BOUNDARY_PSI_INITIAL;
    let psi_current = BOUNDARY_PSI_INITIAL;
    let psi_next = BOUNDARY_PSI_DERIVATIVE;

    let nodes = 0;
    const psiValues: number[] = new Array(this.gridPoints).fill(0);
    psiValues[this.gridPoints - 1] = psi_prev;
    psiValues[this.gridPoints - 2] = psi_current;
    let scaleFactor = 1.0;

    // Renormalization threshold to prevent overflow
    const RENORM_THRESHOLD = 1e10;

    for (let i = this.gridPoints - 3; i >= matchPoint - 1; i--) {
      psi_prev = psi_current;
      psi_current = psi_next;

      // Numerov recursion - indices must match wavefunction indices
      // Going backwards: psi_prev = psi[i+2], psi_current = psi[i+1], psi_next = psi[i]
      psi_next = this.numerovStep(
        psi_current,
        psi_prev,
        effectivePotential[i+2],
        effectivePotential[i+1],
        effectivePotential[i],
        numerovFactor
      );

      // Periodic renormalization to prevent overflow
      if (Math.abs(psi_next) > RENORM_THRESHOLD) {
        const norm = Math.abs(psi_next);
        psi_prev /= norm;
        psi_current /= norm;
        psi_next /= norm;
        scaleFactor *= norm;
        // Renormalize stored values
        for (let j = i + 1; j < this.gridPoints; j++) {
          psiValues[j] /= norm;
        }
      }

      psiValues[i] = psi_next;

      // Count nodes
      if (i >= matchPoint && this.isNode(psi_current, psi_next)) {
        nodes++;
      }
    }

    // Calculate log derivative
    const logDerivative = this.calculateLogDerivative(
      psiValues[matchPoint + 1],
      psiValues[matchPoint],
      psiValues[matchPoint - 1]
    );

    return {
      nodes,
      logDerivative,
      psiAtMatch: psiValues[matchPoint] * scaleFactor
    };
  }

  /**
   * Single Numerov integration step with overflow protection
   */
  private numerovStep(
    psi_current: number,
    psi_prev: number,
    pot_prev: number,
    pot_current: number,
    pot_next: number,
    numerovFactor: number
  ): number {
    const denominator = 1.0 - numerovFactor * pot_next;

    // Prevent division by near-zero values and handle overflow
    if (Math.abs(denominator) < 1e-15) {
      return psi_current * 2.0; // Fallback for numerical stability
    }

    const result = (psi_current * (2.0 + 10.0 * numerovFactor * pot_current) -
            psi_prev * (1.0 - numerovFactor * pot_prev)) / denominator;

    // Cap extreme values to prevent overflow
    const MAX_PSI = 1e100;
    if (!Number.isFinite(result)) {
      return Math.sign(psi_current) * MAX_PSI;
    }
    return Math.max(-MAX_PSI, Math.min(MAX_PSI, result));
  }

  /**
   * Enhanced wavefunction calculation with proper normalization
   */
  private calculateNormalizedWaveFunction(energy: number): number[] {
    const waveFunction = new Array<number>(this.gridPoints);
    const matchPoint = this.adaptiveMatching
      ? this.findOptimalMatchingPoint(energy, 0)
      : Math.floor(this.gridPoints * MATCHING_POINT_FRACTION);

    const numerovFactor = this.getNumerovFactor();
    const inverseHbar2Over2m = 1.0 / this.hbarSquaredOver2m;

    // Calculate effective potential
    const effectivePotential = this.potentialEnergies.map(V =>
      inverseHbar2Over2m * (V - energy)
    );

    // Forward integration
    waveFunction[0] = BOUNDARY_PSI_INITIAL;
    waveFunction[1] = BOUNDARY_PSI_DERIVATIVE;

    for (let i = 2; i <= matchPoint; i++) {
      waveFunction[i] = this.numerovStep(
        waveFunction[i-1],
        waveFunction[i-2],
        effectivePotential[i-2],
        effectivePotential[i-1],
        effectivePotential[i],
        numerovFactor
      );
    }

    const psi_leftAtMatch = waveFunction[matchPoint];

    // Backward integration
    waveFunction[this.gridPoints - 1] = BOUNDARY_PSI_INITIAL;
    waveFunction[this.gridPoints - 2] = BOUNDARY_PSI_DERIVATIVE;

    for (let i = this.gridPoints - 3; i >= matchPoint; i--) {
      waveFunction[i] = this.numerovStep(
        waveFunction[i+1],
        waveFunction[i+2],
        effectivePotential[i+2],
        effectivePotential[i+1],
        effectivePotential[i],
        numerovFactor
      );
    }

    // Scale to match at matching point
    const scaleFactor = psi_leftAtMatch / waveFunction[matchPoint];
    for (let i = matchPoint; i < this.gridPoints; i++) {
      waveFunction[i] *= scaleFactor;
    }

    // Normalize
    return this.normalizeWaveFunction(waveFunction);
  }

  /**
   * Normalizes wavefunction using specified method
   */
  private normalizeWaveFunction(waveFunction: number[]): number[] {
    let norm: number;

    if (this.normalizationMethod === 'l2') {
      // L2 normalization using Simpson's rule
      norm = Math.sqrt(this.calculateL2Norm(waveFunction));
    } else {
      // Max normalization
      norm = Math.max(...waveFunction.map(Math.abs));
    }

    if (norm === 0) {
      return waveFunction;
    }

    return waveFunction.map(psi => psi / norm);
  }

  // ========================================================================
  // Helper Methods
  // ========================================================================

  /**
   * Validates constructor inputs
   */
  private validateInputs(xMin: number, xMax: number, numPoints: number): void {
    if (numPoints < 10) {
      throw new QuantumSolverException(
        'Grid points must be at least 10',
        'INVALID_GRID_POINTS'
      );
    }

    if (xMax <= xMin) {
      throw new QuantumSolverException(
        'xMax must be greater than xMin',
        'INVALID_BOUNDS'
      );
    }

    if (this.hbarSquaredOver2m <= 0) {
      throw new QuantumSolverException(
        'hbarSquaredOver2m must be positive',
        'INVALID_PHYSICS_CONSTANT'
      );
    }
  }

  /**
   * Calculates grid positions
   */
  private calculateGridPositions(xMin: number, xMax: number, numPoints: number): number[] {
    const positions = new Array<number>(numPoints);
    const dx = (xMax - xMin) / (numPoints - 1);
    for (let i = 0; i < numPoints; i++) {
      positions[i] = xMin + i * dx;
    }
    return positions;
  }

  /**
   * Calculates potential on grid
   */
  private calculatePotentialOnGrid(): number[] {
    const potential = new Array<number>(this.gridPoints);

    // Sample all points including boundaries
    for (let i = 0; i < this.gridPoints; i++) {
      potential[i] = this.potentialFunction(this.gridPositions[i]);
    }

    return potential;
  }

  /**
   * Gets cached Numerov factor
   */
  private getNumerovFactor(): number {
    const key = this.deltaX;
    if (!this.numerovFactorCache.has(key)) {
      this.numerovFactorCache.set(key, this.deltaX * this.deltaX * NUMEROV_COEFFICIENT);
    }
    return this.numerovFactorCache.get(key)!;
  }

  /**
   * Checks if there's a node between two points
   */
  private isNode(psi1: number, psi2: number): boolean {
    return (psi1 * psi2 < 0) &&
           Math.abs(psi1) > this.waveFunctionTolerance &&
           Math.abs(psi2) > this.waveFunctionTolerance;
  }

  /**
   * Calculates logarithmic derivative
   */
  private calculateLogDerivative(psiBefore: number, psiAt: number, psiAfter: number): number {
    if (Math.abs(psiAt) < this.waveFunctionTolerance) {
      return 0;
    }
    return (psiAfter - psiBefore) / (2.0 * this.deltaX * psiAt);
  }

  /**
   * Finds optimal matching point for given energy
   */
  private findOptimalMatchingPoint(energy: number, _desiredNodes: number): number {
    // Find classical turning points
    const turningPoints: number[] = [];
    for (let i = 1; i < this.gridPoints - 1; i++) {
      if (this.potentialEnergies[i] <= energy &&
          this.potentialEnergies[i+1] > energy) {
        turningPoints.push(i);
      }
    }

    // Use middle turning point or default
    if (turningPoints.length >= 2) {
      return turningPoints[Math.floor(turningPoints.length / 2)];
    }

    return Math.floor(this.gridPoints * MATCHING_POINT_FRACTION);
  }

  /**
   * Calculates L2 norm using Simpson's rule
   */
  private calculateL2Norm(waveFunction: number[]): number {
    let sum = 0;
    for (let i = 0; i < this.gridPoints - 2; i += 2) {
      sum += waveFunction[i] ** 2 +
             4 * waveFunction[i+1] ** 2 +
             waveFunction[i+2] ** 2;
    }
    return sum * this.deltaX / 3.0;
  }

  /**
   * Calculates normalization constant
   */
  private calculateNormalizationConstant(waveFunction: number[]): number {
    const l2Norm = this.calculateL2Norm(waveFunction);
    return l2Norm > 0 ? 1.0 / Math.sqrt(l2Norm) : 1.0;
  }

  /**
   * Calculates convergence quality metric
   */
  private calculateConvergenceMetric(
    energy: number,
    _waveFunction: number[],
    desiredNodes: number
  ): number {
    const result = this.integrateSchrodinger(energy, desiredNodes);
    const continuity = Math.exp(-Math.abs(result.logDerivativeMismatch));
    const amplitude = result.matchingAmplitude;
    return continuity * amplitude;
  }

  // ========================================================================
  // Cache Management
  // ========================================================================

  /**
   * Gets value from cache
   */
  private getFromCache<T>(key: string): T | null {
    if (!this.enableCaching) return null;

    const entry = this.cache.get(key);
    if (entry) {
      this.cacheHits++;
      return entry.value as T;
    }

    this.cacheMisses++;
    return null;
  }

  /**
   * Adds value to cache with size limit
   */
  private addToCache<T>(key: string, value: T): void {
    if (!this.enableCaching) return;

    // Enforce cache size limit
    if (this.cache.size >= CACHE_SIZE_LIMIT) {
      // Remove oldest entry
      const firstKey = this.cache.keys().next().value;
      if (firstKey !== undefined) {
        this.cache.delete(firstKey);
      }
    }

    this.cache.set(key, {
      key,
      value,
      timestamp: Date.now()
    });
  }
}

// ============================================================================
// Functional API (matches NumerovSolver pattern)
// ============================================================================

/**
 * Solve the 1D Schrodinger equation using the enhanced shooting method.
 *
 * @param potential - Function V(x) that returns potential energy in Joules
 * @param mass - Particle mass in kg
 * @param numStates - Number of bound states to find
 * @param gridConfig - Grid configuration
 * @param energyMin - Minimum energy to search (Joules) - used for validation
 * @param energyMax - Maximum energy to search (Joules) - used for validation
 * @param config - Optional solver configuration
 * @returns Bound state results
 */
export function solveQuantumBound(
  potential: PotentialFunction,
  mass: number,
  numStates: number,
  gridConfig: GridConfig,
  energyMin?: number,
  energyMax?: number,
  config?: SolverConfig
): BoundStateResult {
  const { xMin, xMax, numPoints } = gridConfig;

  // Create solver
  const solver = new QuantumBoundStateSolver(
    mass,
    xMin,
    xMax,
    numPoints,
    potential,
    config
  );

  // Find all eigenstates
  const states = solver.findMultipleEigenstates(numStates);

  // Convert to BoundStateResult format
  const energies: number[] = [];
  const wavefunctions: number[][] = [];

  for (const state of states) {
    // Filter by energy range if provided
    if (energyMin !== undefined && state.energy < energyMin) continue;
    if (energyMax !== undefined && state.energy > energyMax) continue;

    energies.push(state.energy);
    wavefunctions.push(state.waveFunction);
  }

  return {
    energies,
    wavefunctions,
    xGrid: solver.getGridPositions(),
    method: "numerov"  // Compatible with existing method types
  };
}

// ============================================================================
// Utility Functions
// ============================================================================

/**
 * Creates a solver with atomic units (Hartree, Bohr radii)
 */
export function createAtomicUnitsSolver(
  xMin: number,  // Bohr radii
  xMax: number,  // Bohr radii
  numPoints: number,
  potentialFunction: (x: number) => number,  // Hartree
  config?: SolverConfig
): QuantumBoundStateSolver {
  // In atomic units, hbar = m_e = e = 1, so hbar^2/2m = 0.5
  // Convert to SI: use electron mass
  const { ELECTRON_MASS, BOHR_RADIUS, HARTREE_TO_JOULES } = QuantumConstants;

  // Create wrapper that converts from atomic to SI units
  const siPotential = (x: number) => {
    const xBohr = x / BOHR_RADIUS;
    return potentialFunction(xBohr) * HARTREE_TO_JOULES;
  };

  return new QuantumBoundStateSolver(
    ELECTRON_MASS,
    xMin * BOHR_RADIUS,
    xMax * BOHR_RADIUS,
    numPoints,
    siPotential,
    config
  );
}

// ============================================================================
// Registration
// ============================================================================

qppw.register("QuantumBoundStateSolver", {
  QuantumBoundStateSolver,
  solveQuantumBound,
  createAtomicUnitsSolver,
  QuantumSolverException
});

export default QuantumBoundStateSolver;
