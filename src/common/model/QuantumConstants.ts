/**
 * Physical constants for quantum mechanics calculations.
 * All values are in SI units unless otherwise specified.
 */

import qppw from "../../QPPWNamespace.js";

const QuantumConstants = {
  /**
   * Reduced Planck constant (ℏ) in J·s
   */
  HBAR: 1.054571817e-34,

  /**
   * Electron mass in kg
   */
  ELECTRON_MASS: 9.10938356e-31,

  /**
   * Elementary charge in Coulombs
   */
  ELEMENTARY_CHARGE: 1.602176634e-19,

  /**
   * Electron volt in Joules (for unit conversions)
   */
  EV_TO_JOULES: 1.602176634e-19,

  /**
   * Joules to electron volt (for unit conversions)
   */
  JOULES_TO_EV: 1.0 / 1.602176634e-19,

  /**
   * Nanometer to meters conversion
   */
  NM_TO_M: 1e-9,

  /**
   * Meters to nanometer conversion
   */
  M_TO_NM: 1e9,
} as const;

qppw.register("QuantumConstants", QuantumConstants);

export default QuantumConstants;
