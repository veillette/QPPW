/**
 * Type definitions for superposition states.
 */

import qppw from "../../QPPWNamespace.js";

/**
 * Type of superposition state
 */
export enum SuperpositionType {
  PSI_I_PSI_J = "psiIPsiJ",
  PSI_K = "psiK",
  LOCALIZED_NARROW = "localizedNarrow",
  LOCALIZED_WIDE = "localizedWide",
  COHERENT = "coherent",
  CUSTOM = "custom",
}

/**
 * Configuration for a superposition state
 */
export interface SuperpositionConfig {
  /** Type of superposition */
  type: SuperpositionType;
  /** Amplitudes for each eigenstate (must sum to 1 when squared) */
  amplitudes: number[];
  /** Phases for each eigenstate (in radians) */
  phases: number[];
  /** Displacement from equilibrium in nm (for coherent states) */
  displacement?: number;
}

qppw.register("SuperpositionType", { SuperpositionType });

export default SuperpositionType;
