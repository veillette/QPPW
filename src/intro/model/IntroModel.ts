/**
 * IntroModel represents the physics model for the intro screen.
 * It extends BaseModel and provides a simplified interface for introductory exploration.
 */

import { NumberProperty, Property } from "scenerystack/axon";
import { Range } from "scenerystack/dot";
import { BaseModel } from "../../common/model/BaseModel.js";
import {
    WellParameters,
    NumericalMethod,
} from "../../common/model/Schrodinger1DSolver.js";
import {
    PotentialType,
    BoundStateResult,
} from "../../common/model/PotentialFunction.js";
import QuantumConstants from "../../common/model/QuantumConstants.js";
import { SuperpositionType, SuperpositionConfig } from "../../common/model/SuperpositionType.js";

export type DisplayMode = "probabilityDensity" | "waveFunction";

export class IntroModel extends BaseModel {
    // Potential type selection
    public readonly potentialTypeProperty: Property<PotentialType>;

    // Well parameters
    public readonly wellWidthProperty: NumberProperty;
    public readonly wellDepthProperty: NumberProperty;
    public readonly wellOffsetProperty: NumberProperty;
    public readonly barrierHeightProperty: NumberProperty;
    public readonly potentialOffsetProperty: NumberProperty;

    // Particle properties
    public readonly particleMassProperty: NumberProperty;

    // Energy level selection
    public readonly selectedEnergyLevelIndexProperty: NumberProperty;

    // Display settings
    public readonly displayModeProperty: Property<DisplayMode>;
    public readonly showRealPartProperty: Property<boolean>;
    public readonly showImaginaryPartProperty: Property<boolean>;
    public readonly showMagnitudeProperty: Property<boolean>;
    public readonly showPhaseProperty: Property<boolean>;

    // Superposition properties (required by WaveFunctionChartNode, but not used in intro screen)
    public readonly superpositionTypeProperty: Property<SuperpositionType>;
    public readonly superpositionConfigProperty: Property<SuperpositionConfig>;

    // Chart visibility
    public readonly showTotalEnergyProperty: Property<boolean>;
    public readonly showPotentialEnergyProperty: Property<boolean>;

    // Cached bound state results
    protected boundStateResult: BoundStateResult | null = null;

    public constructor() {
        super();

        // Initialize potential type
        this.potentialTypeProperty = new Property<PotentialType>(
            PotentialType.INFINITE_WELL,
        );

        // Initialize well parameters
        this.wellWidthProperty = new NumberProperty(1.0, {
            range: new Range(0.1, 6.0),
        });
        this.wellDepthProperty = new NumberProperty(5.0, {
            range: new Range(0.1, 15.0),
        });
        this.wellOffsetProperty = new NumberProperty(0.5, {
            range: new Range(0.0, 1.0),
        });
        this.barrierHeightProperty = new NumberProperty(0.5, {
            range: new Range(0.0, 10.0),
        });
        this.potentialOffsetProperty = new NumberProperty(0.0, {
            range: new Range(-5.0, 15.0),
        });

        // Initialize particle mass
        this.particleMassProperty = new NumberProperty(1.0, {
            range: new Range(0.5, 1.1),
        });

        // Initialize energy level selection
        this.selectedEnergyLevelIndexProperty = new NumberProperty(0, {
            range: new Range(0, 99),
        });

        // Initialize display settings
        this.displayModeProperty = new Property<DisplayMode>("probabilityDensity");
        this.showRealPartProperty = new Property<boolean>(true);
        this.showImaginaryPartProperty = new Property<boolean>(false);
        this.showMagnitudeProperty = new Property<boolean>(false);
        this.showPhaseProperty = new Property<boolean>(false);

        // Initialize superposition properties (not used in intro screen, but required by chart nodes)
        this.superpositionTypeProperty = new Property<SuperpositionType>(SuperpositionType.SINGLE);
        this.superpositionConfigProperty = new Property<SuperpositionConfig>({
            type: SuperpositionType.SINGLE,
            amplitudes: [],
            phases: [],
        });

        // Initialize chart visibility
        this.showTotalEnergyProperty = new Property<boolean>(true);
        this.showPotentialEnergyProperty = new Property<boolean>(true);

        // Recalculate bound states when parameters change
        const invalidateCache = () => {
            this.boundStateResult = null;
        };

        this.potentialTypeProperty.link(invalidateCache);
        this.wellWidthProperty.link(invalidateCache);
        this.wellDepthProperty.link(invalidateCache);
        this.wellOffsetProperty.link(invalidateCache);
        this.barrierHeightProperty.link(invalidateCache);
        this.potentialOffsetProperty.link(invalidateCache);
        this.particleMassProperty.link(invalidateCache);
    }

    /**
     * Called when the solver method or grid points changes.
     */
    protected onSolverMethodChanged(_method: NumericalMethod): void {
        this.boundStateResult = null;
    }

    /**
     * Resets the model to its initial state.
     */
    public reset(): void {
        this.resetAll();
    }

    /**
     * Resets all properties to their initial state.
     */
    public override resetAll(): void {
        super.resetAll();
        this.potentialTypeProperty.reset();
        this.wellWidthProperty.reset();
        this.wellDepthProperty.reset();
        this.wellOffsetProperty.reset();
        this.barrierHeightProperty.reset();
        this.potentialOffsetProperty.reset();
        this.particleMassProperty.reset();
        this.selectedEnergyLevelIndexProperty.reset();
        this.displayModeProperty.reset();
        this.showRealPartProperty.reset();
        this.showImaginaryPartProperty.reset();
        this.showMagnitudeProperty.reset();
        this.showPhaseProperty.reset();
        this.superpositionTypeProperty.reset();
        this.superpositionConfigProperty.reset();
        this.showTotalEnergyProperty.reset();
        this.showPotentialEnergyProperty.reset();
    }

    /**
     * Get all bound state energies and wavefunctions for the current well.
     */
    public getBoundStates(): BoundStateResult | null {
        if (!this.boundStateResult) {
            this.calculateBoundStates();
        }
        return this.boundStateResult;
    }

    /**
     * Calculate bound states using the Schrödinger solver.
     */
    private calculateBoundStates(): void {
        const wellWidth = this.wellWidthProperty.value * QuantumConstants.NM_TO_M;
        const wellDepth =
            this.wellDepthProperty.value * QuantumConstants.EV_TO_JOULES;
        const mass =
            this.particleMassProperty.value * QuantumConstants.ELECTRON_MASS;

        const numStates = 10;
        const CHART_DISPLAY_RANGE_NM = 4;
        const ANALYTICAL_GRID_POINTS = 1000;

        const gridConfig = {
            xMin: -CHART_DISPLAY_RANGE_NM * QuantumConstants.NM_TO_M,
            xMax: CHART_DISPLAY_RANGE_NM * QuantumConstants.NM_TO_M,
            numPoints: ANALYTICAL_GRID_POINTS,
        };

        try {
            const potentialParams: WellParameters = {
                type: this.potentialTypeProperty.value,
                wellWidth: wellWidth,
            };

            // Add type-specific parameters
            switch (this.potentialTypeProperty.value) {
                case PotentialType.FINITE_WELL:
                    potentialParams.wellDepth = wellDepth;
                    break;
                case PotentialType.HARMONIC_OSCILLATOR:
                    potentialParams.springConstant =
                        (8 * wellDepth) / (wellWidth * wellWidth);
                    break;
                case PotentialType.MORSE:
                    potentialParams.dissociationEnergy = wellDepth;
                    potentialParams.equilibriumPosition = 0;
                    potentialParams.wellWidth = wellWidth;
                    break;
                case PotentialType.POSCHL_TELLER:
                    potentialParams.potentialDepth = wellDepth;
                    potentialParams.wellWidth = wellWidth;
                    break;
                case PotentialType.ROSEN_MORSE:
                    potentialParams.potentialDepth = wellDepth;
                    potentialParams.barrierHeight =
                        this.barrierHeightProperty.value * QuantumConstants.EV_TO_JOULES;
                    potentialParams.wellWidth = wellWidth;
                    break;
                case PotentialType.ECKART:
                    potentialParams.potentialDepth = wellDepth;
                    potentialParams.barrierHeight =
                        this.barrierHeightProperty.value * QuantumConstants.EV_TO_JOULES;
                    potentialParams.wellWidth = wellWidth;
                    break;
                case PotentialType.ASYMMETRIC_TRIANGLE:
                    potentialParams.slope = wellDepth / wellWidth;
                    break;
                case PotentialType.TRIANGULAR:
                    potentialParams.potentialDepth = wellDepth;
                    potentialParams.wellWidth = wellWidth;
                    potentialParams.energyOffset =
                        this.potentialOffsetProperty.value * QuantumConstants.EV_TO_JOULES;
                    break;
                case PotentialType.COULOMB_1D:
                case PotentialType.COULOMB_3D: {
                    const coulombConstant = 8.9875517923e9;
                    potentialParams.coulombStrength =
                        coulombConstant *
                        QuantumConstants.ELEMENTARY_CHARGE *
                        QuantumConstants.ELEMENTARY_CHARGE;
                    break;
                }
            }

            this.boundStateResult = this.solver.solveAnalyticalIfPossible(
                potentialParams,
                mass,
                numStates,
                gridConfig,
            );

            // Ensure selected energy level index is within bounds
            if (this.boundStateResult) {
                const maxIndex = this.boundStateResult.energies.length - 1;
                if (this.selectedEnergyLevelIndexProperty.value > maxIndex) {
                    this.selectedEnergyLevelIndexProperty.value = Math.max(0, maxIndex);
                }
            }
        } catch (error) {
            console.error("Error calculating bound states:", error);
            this.boundStateResult = null;
        }
    }

    /**
     * Calculate the classical probability density for a given energy level.
     */
    public getClassicalProbabilityDensity(energyIndex: number): number[] | null {
        if (!this.boundStateResult) {
            this.calculateBoundStates();
        }

        if (
            !this.boundStateResult ||
            energyIndex < 0 ||
            energyIndex >= this.boundStateResult.energies.length
        ) {
            return null;
        }

        const energy = this.boundStateResult.energies[energyIndex];
        const xGrid = this.boundStateResult.xGrid;
        const mass = this.particleMassProperty.value * QuantumConstants.ELECTRON_MASS;

        // Fallback: numerical calculation using potential function
        const potential = this.calculatePotentialEnergy(xGrid);

        // Use BaseModel's common method to calculate classical probability density
        return this.calculateClassicalProbabilityDensity(potential, energy, mass, xGrid);
    }

    /**
     * Calculates the classical turning points for a given energy level.
     */
    public getClassicalTurningPoints(energyLevel: number): {
        left: number;
        right: number;
    } | null {
        if (!this.boundStateResult) {
            this.calculateBoundStates();
        }

        if (
            !this.boundStateResult ||
            energyLevel < 0 ||
            energyLevel >= this.boundStateResult.energies.length
        ) {
            return null;
        }

        const energy = this.boundStateResult.energies[energyLevel];
        const energyEV = energy * QuantumConstants.JOULES_TO_EV;
        const xGrid = this.boundStateResult.xGrid;

        let leftTurningPoint: number | null = null;
        let rightTurningPoint: number | null = null;

        // Handle infinite well directly (classical turning points at well edges)
        if (this.potentialTypeProperty.value === PotentialType.INFINITE_WELL) {
            const halfWidth = this.wellWidthProperty.value / 2;
            return { left: -halfWidth, right: halfWidth };
        }

        // Existing loop to find turning points for other potentials
        for (let i = 0; i < xGrid.length - 1; i++) {
            const x = xGrid[i] * QuantumConstants.M_TO_NM;
            const V = this.getPotentialAtPosition(x);

            const xNext = xGrid[i + 1] * QuantumConstants.M_TO_NM;
            const VNext = this.getPotentialAtPosition(xNext);

            if (
                (V <= energyEV && VNext >= energyEV) ||
                (V >= energyEV && VNext <= energyEV)
            ) {
                // Avoid division by zero when V and VNext are equal
                if (VNext === V) {
                    continue; // skip this segment
                }
                const t = (energyEV - V) / (VNext - V);
                const turningPoint = x + t * (xNext - x);

                if (leftTurningPoint === null) {
                    leftTurningPoint = turningPoint;
                } else {
                    rightTurningPoint = turningPoint;
                }
            }
        }

        if (leftTurningPoint !== null && rightTurningPoint !== null) {
            // Clamp turning points to chart's X-axis range (-X_AXIS_RANGE_NM to X_AXIS_RANGE_NM)
            const X_AXIS_RANGE_NM = 4; // Define X_AXIS_RANGE_NM here
            const minX = -X_AXIS_RANGE_NM;
            const maxX = X_AXIS_RANGE_NM;
            leftTurningPoint = Math.max(minX, Math.min(maxX, leftTurningPoint));
            rightTurningPoint = Math.max(minX, Math.min(maxX, rightTurningPoint));
            return { left: leftTurningPoint, right: rightTurningPoint };
        }

        return null;
    }

    /**
     * Calculates the probability of finding the particle in the classically forbidden region.
     */
    public getClassicallyForbiddenProbability(energyLevel: number): number {
        if (!this.boundStateResult) {
            this.calculateBoundStates();
        }

        if (
            !this.boundStateResult ||
            energyLevel < 0 ||
            energyLevel >= this.boundStateResult.energies.length
        ) {
            return 0;
        }

        const turningPoints = this.getClassicalTurningPoints(energyLevel);
        if (!turningPoints) {
            return 0;
        }

        const wavefunction = this.boundStateResult.wavefunctions[energyLevel];
        const xGrid = this.boundStateResult.xGrid;

        let forbiddenProbability = 0;
        let totalProbability = 0;

        for (let i = 0; i < xGrid.length; i++) {
            const x = xGrid[i] * QuantumConstants.M_TO_NM;
            const psi = wavefunction[i];
            const probabilityDensity = psi * psi;

            let dx = 0;
            if (i === 0) {
                dx = ((xGrid[1] - xGrid[0]) / 2) * QuantumConstants.M_TO_NM;
            } else if (i === xGrid.length - 1) {
                dx = ((xGrid[i] - xGrid[i - 1]) / 2) * QuantumConstants.M_TO_NM;
            } else {
                dx = ((xGrid[i + 1] - xGrid[i - 1]) / 2) * QuantumConstants.M_TO_NM;
            }

            totalProbability += probabilityDensity * dx;

            if (x < turningPoints.left || x > turningPoints.right) {
                forbiddenProbability += probabilityDensity * dx;
            }
        }

        return totalProbability > 0
            ? (forbiddenProbability / totalProbability) * 100
            : 0;
    }

    /**
     * Calculate the potential energy at given positions.
     */
    private calculatePotentialEnergy(xGrid: number[]): number[] {
        const wellWidth = this.wellWidthProperty.value * QuantumConstants.NM_TO_M;
        const wellDepth = this.wellDepthProperty.value * QuantumConstants.EV_TO_JOULES;

        const potential: number[] = [];

        for (let i = 0; i < xGrid.length; i++) {
            const x = xGrid[i];
            let V = 0;

            switch (this.potentialTypeProperty.value) {
                case PotentialType.INFINITE_WELL:
                    V = Math.abs(x) <= wellWidth / 2 ? 0 : Infinity;
                    break;

                case PotentialType.FINITE_WELL:
                    V = Math.abs(x) <= wellWidth / 2 ? 0 : wellDepth;
                    break;

                case PotentialType.HARMONIC_OSCILLATOR: {
                    const springConstant = (8 * wellDepth) / (wellWidth * wellWidth);
                    V = 0.5 * springConstant * x * x;
                    break;
                }

                case PotentialType.MORSE: {
                    const a = 1 / wellWidth;
                    const exponential = Math.exp(-a * x);
                    V = wellDepth * Math.pow(1 - exponential, 2);
                    break;
                }

                case PotentialType.POSCHL_TELLER: {
                    const coshPT = Math.cosh(x / wellWidth);
                    V = -wellDepth / (coshPT * coshPT);
                    break;
                }

                case PotentialType.ROSEN_MORSE: {
                    const barrierHeight =
                        this.barrierHeightProperty.value * QuantumConstants.EV_TO_JOULES;
                    const coshRM = Math.cosh(x / wellWidth);
                    const tanhRM = Math.tanh(x / wellWidth);
                    V = -wellDepth / (coshRM * coshRM) + barrierHeight * tanhRM;
                    break;
                }

                case PotentialType.ECKART: {
                    const barrierHeightE =
                        this.barrierHeightProperty.value * QuantumConstants.EV_TO_JOULES;
                    const expE = Math.exp(x / wellWidth);
                    const denomE = 1 + expE;
                    V = wellDepth / (denomE * denomE) - barrierHeightE / denomE;
                    break;
                }

                case PotentialType.ASYMMETRIC_TRIANGLE: {
                    const slope = wellDepth / wellWidth;
                    V = slope * x;
                    break;
                }

                case PotentialType.TRIANGULAR: {
                    const energyOffset =
                        this.potentialOffsetProperty.value * QuantumConstants.EV_TO_JOULES;
                    if (x < 0) {
                        V = wellDepth + energyOffset;
                    } else if (x < wellWidth) {
                        V = energyOffset + (wellDepth / wellWidth) * x;
                    } else {
                        V = wellDepth + energyOffset;
                    }
                    break;
                }

                case PotentialType.COULOMB_1D:
                case PotentialType.COULOMB_3D: {
                    const coulombConstant = 8.9875517923e9;
                    const coulombStrength =
                        coulombConstant *
                        QuantumConstants.ELEMENTARY_CHARGE *
                        QuantumConstants.ELEMENTARY_CHARGE;
                    const r = Math.abs(x);
                    if (r > 1e-12) {
                        V = -coulombStrength / r;
                    } else {
                        V = -coulombStrength / 1e-12;
                    }
                    break;
                }

                default:
                    V = 0;
                    break;
            }

            potential.push(V);
        }

        return potential;
    }

    /**
     * Get the potential energy at a specific position (in nm).
     */
    private getPotentialAtPosition(xNm: number): number {
        const x = xNm * QuantumConstants.NM_TO_M;
        const potential = this.calculatePotentialEnergy([x]);
        return potential[0] * QuantumConstants.JOULES_TO_EV;
    }
}
