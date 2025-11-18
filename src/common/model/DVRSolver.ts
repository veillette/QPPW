/**
 * Discrete Variable Representation (DVR) method for solving the 1D
 * time-independent Schrödinger equation.
 *
 * This method constructs the Hamiltonian matrix H = T + V and diagonalizes it:
 * - Potential energy V is diagonal: V_ii = V(x_i)
 * - Kinetic energy T uses Colbert-Miller formula for uniform grids
 *
 * Reference: Colbert & Miller, J. Chem. Phys. 96, 1982 (1992)
 */

import QuantumConstants from "./QuantumConstants.js";
import { BoundStateResult, GridConfig, PotentialFunction } from "./PotentialFunction.js";
import { Matrix, diagonalize, normalizeWavefunction, matrixToArray } from "./LinearAlgebraUtils.js";
import qppw from "../../QPPWNamespace.js";

/**
 * Solve the 1D Schrödinger equation using DVR method with Colbert-Miller kinetic energy.
 *
 * @param potential - Function V(x) that returns potential energy in Joules
 * @param mass - Particle mass in kg
 * @param numStates - Number of lowest bound states to return
 * @param gridConfig - Grid configuration
 * @returns Bound state results
 */
export function solveDVR(
  potential: PotentialFunction,
  mass: number,
  numStates: number,
  gridConfig: GridConfig,
): BoundStateResult {
  const { xMin, xMax, numPoints } = gridConfig;
  const dx = (xMax - xMin) / (numPoints - 1);

  // Generate grid
  const xGrid: number[] = [];
  for (let i = 0; i < numPoints; i++) {
    xGrid.push(xMin + i * dx);
  }

  const N = numPoints;

  // Create potential energy values for diagonal matrix
  const potentialValues: number[] = xGrid.map(potential);

  // Create potential energy matrix V (diagonal) using dot's Matrix
  const V = Matrix.diagonalMatrix(potentialValues);

  // Construct kinetic energy matrix T using Colbert-Miller formula
  const T = createKineticEnergyMatrix(N, dx, mass);

  // Construct Hamiltonian H = T + V using dot's Matrix.plus()
  const H = T.plus(V);

  // Convert to array for diagonalization
  const hamiltonianArray = matrixToArray(H);

  // Diagonalize Hamiltonian to get eigenvalues (energies) and eigenvectors (wavefunctions)
  const eigen = diagonalize(hamiltonianArray);

  // Sort eigenvalues and eigenvectors by energy
  const sortedIndices = eigen.eigenvalues
    .map((e, i) => ({ energy: e, index: i }))
    .sort((a, b) => a.energy - b.energy)
    .map((item) => item.index);

  // Extract the lowest numStates bound states
  const energies: number[] = [];
  const wavefunctions: number[][] = [];

  const V_boundary = Math.max(potential(xMin), potential(xMax));

  for (let i = 0; i < Math.min(numStates, N); i++) {
    const idx = sortedIndices[i];
    const energy = eigen.eigenvalues[idx];

    // Only include bound states (energy < V at boundaries)
    if (energy < V_boundary) {
      energies.push(energy);

      // Normalize wavefunction
      const wavefunction = eigen.eigenvectors[idx];
      const normalizedPsi = normalizeWavefunction(wavefunction, dx);
      wavefunctions.push(normalizedPsi);
    }
  }

  return {
    energies,
    wavefunctions,
    xGrid,
    method: "dvr",
  };
}

/**
 * Create kinetic energy matrix using Colbert-Miller formula:
 * T_ij = (ℏ²/2m) * (1/Δx²) * { π²/3 for i=j, 2*(-1)^(i-j)/(i-j)² for i≠j }
 *
 * @param N - Matrix dimension (number of grid points)
 * @param dx - Grid spacing (meters)
 * @param mass - Particle mass (kg)
 * @returns N×N kinetic energy matrix
 */
function createKineticEnergyMatrix(N: number, dx: number, mass: number): Matrix {
  const { HBAR } = QuantumConstants;
  const prefactor = (HBAR * HBAR) / (2 * mass * dx * dx);

  const T = new Matrix(N, N);
  for (let i = 0; i < N; i++) {
    for (let j = 0; j < N; j++) {
      if (i === j) {
        // Diagonal elements
        T.set(i, j, prefactor * (Math.PI * Math.PI) / 3);
      } else {
        // Off-diagonal elements
        const diff = i - j;
        const sign = Math.pow(-1, diff);
        T.set(i, j, prefactor * (2 * sign) / (diff * diff));
      }
    }
  }

  return T;
}

qppw.register("DVRSolver", { solveDVR });
