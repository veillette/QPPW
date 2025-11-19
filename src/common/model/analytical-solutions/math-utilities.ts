/**
 * Mathematical utility functions for analytical solutions.
 * These are used by various potential solvers for special functions and polynomials.
 */

/**
 * Calculate factorial n!
 */
export function factorial(n: number): number {
  if (n <= 1) return 1;
  let result = 1;
  for (let i = 2; i <= n; i++) {
    result *= i;
  }
  return result;
}

/**
 * Calculate the Hermite polynomial H_n(x) using recurrence relation.
 * H_0(x) = 1
 * H_1(x) = 2x
 * H_(n+1)(x) = 2x*H_n(x) - 2n*H_(n-1)(x)
 */
export function hermitePolynomial(n: number, x: number): number {
  if (n === 0) return 1;
  if (n === 1) return 2 * x;

  let H_prev = 1;
  let H_curr = 2 * x;

  for (let i = 1; i < n; i++) {
    const H_next = 2 * x * H_curr - 2 * i * H_prev;
    H_prev = H_curr;
    H_curr = H_next;
  }

  return H_curr;
}

/**
 * Calculate the associated Laguerre polynomial L_n^α(x) using recurrence relation.
 */
export function associatedLaguerre(n: number, alpha: number, x: number): number {
  if (n === 0) return 1;
  if (n === 1) return 1 + alpha - x;

  let L_prev = 1;
  let L_curr = 1 + alpha - x;

  for (let k = 1; k < n; k++) {
    const L_next =
      ((2 * k + 1 + alpha - x) * L_curr - (k + alpha) * L_prev) / (k + 1);
    L_prev = L_curr;
    L_curr = L_next;
  }

  return L_curr;
}

/**
 * Calculate the Jacobi polynomial P_n^(α,β)(x) using recurrence relation.
 */
export function jacobiPolynomial(n: number, alpha: number, beta: number, x: number): number {
  if (n === 0) return 1;
  if (n === 1) return 0.5 * (alpha - beta + (alpha + beta + 2) * x);

  let P_prev = 1;
  let P_curr = 0.5 * (alpha - beta + (alpha + beta + 2) * x);

  for (let k = 1; k < n; k++) {
    const a1 = 2 * (k + 1) * (k + alpha + beta + 1) * (2 * k + alpha + beta);
    const a2 = (2 * k + alpha + beta + 1) * (alpha * alpha - beta * beta);
    const a3 = (2 * k + alpha + beta) * (2 * k + alpha + beta + 1) * (2 * k + alpha + beta + 2);
    const a4 = 2 * (k + alpha) * (k + beta) * (2 * k + alpha + beta + 2);

    const P_next = ((a2 + a3 * x) * P_curr - a4 * P_prev) / a1;
    P_prev = P_curr;
    P_curr = P_next;
  }

  return P_curr;
}

/**
 * Gamma function approximation using Stirling's formula for large n
 * and direct calculation for small n.
 */
export function gamma(n: number): number {
  // For integer or half-integer values
  if (n === Math.floor(n)) {
    // Integer
    return factorial(n - 1);
  } else if (n - 0.5 === Math.floor(n - 0.5)) {
    // Half-integer: Γ(n+1/2) = sqrt(π) * (2n)! / (4^n * n!)
    const k = n - 0.5;
    return Math.sqrt(Math.PI) * factorial(2 * k) / (Math.pow(4, k) * factorial(k));
  } else {
    // Use Stirling's approximation: Γ(n) ≈ sqrt(2π/n) * (n/e)^n
    return Math.sqrt(2 * Math.PI / n) * Math.pow(n / Math.E, n);
  }
}

/**
 * Log-gamma function for better numerical stability.
 */
export function logGamma(x: number): number {
  if (x <= 0) {
    throw new Error("logGamma undefined for non-positive arguments");
  }

  // For small integers, use log(factorial)
  if (x === Math.floor(x) && x < 20) {
    return Math.log(factorial(x - 1));
  }

  // Stirling's approximation: log(Γ(x)) ≈ (x-0.5)*log(x) - x + 0.5*log(2π)
  return (x - 0.5) * Math.log(x) - x + 0.5 * Math.log(2 * Math.PI);
}

/**
 * Airy function Ai(x) using series expansion for small |x| and asymptotic form for large |x|.
 * Ai(x) is the solution to y'' - xy = 0 that decays exponentially for large positive x.
 */
export function airyAi(x: number): number {
  const ABS_X_THRESHOLD = 5.0;

  if (Math.abs(x) < ABS_X_THRESHOLD) {
    // Series expansion for small |x|
    // Ai(x) = c1 * (1 + x^3/(2*3) + x^6/(2*3*5*6) + ...) - c2 * (x + x^4/(3*4) + x^7/(3*4*6*7) + ...)
    // where c1 = 1/(3^(2/3)*Γ(2/3)), c2 = 1/(3^(1/3)*Γ(1/3))
    const c1 = 0.3550280538878172; // 1/(3^(2/3)*Γ(2/3))
    const c2 = 0.2588194037928068; // 1/(3^(1/3)*Γ(1/3))

    let term1 = 1.0;
    let sum1 = 1.0;
    for (let k = 1; k <= 30; k++) {
      term1 *= (x * x * x) / ((3 * k - 1) * (3 * k));
      sum1 += term1;
      if (Math.abs(term1) < 1e-15) break;
    }

    let term2 = x;
    let sum2 = x;
    for (let k = 1; k <= 30; k++) {
      term2 *= (x * x * x) / ((3 * k) * (3 * k + 1));
      sum2 += term2;
      if (Math.abs(term2) < 1e-15) break;
    }

    return c1 * sum1 - c2 * sum2;
  } else if (x > 0) {
    // Asymptotic expansion for large positive x with higher-order corrections
    // Ai(x) ≈ (1/(2√π)) * x^(-1/4) * exp(-ζ) * (1 - c₁/ζ + c₂/ζ² - ...)
    // where ζ = (2/3) * x^(3/2)
    const zeta = (2.0 / 3.0) * Math.pow(x, 1.5);
    const factor = 0.5 / Math.sqrt(Math.PI) * Math.pow(x, -0.25);

    // Asymptotic series coefficients: cₖ = Γ(3k + 1/2) / (54^k * k! * Γ(k + 1/2))
    // First few: c₁ = 5/72, c₂ = 385/10368
    const invZeta = 1.0 / zeta;
    const correction = 1.0 - (5.0 / 72.0) * invZeta + (385.0 / 10368.0) * invZeta * invZeta;

    return factor * Math.exp(-zeta) * correction;
  } else {
    // Asymptotic expansion for large negative x
    // Ai(x) ≈ (1/√π) * |x|^(-1/4) * sin(ζ + π/4)
    // where ζ = (2/3) * |x|^(3/2)
    const absX = Math.abs(x);
    const zeta = (2.0 / 3.0) * Math.pow(absX, 1.5);
    const factor = 1.0 / Math.sqrt(Math.PI) * Math.pow(absX, -0.25);
    return factor * Math.sin(zeta + Math.PI / 4);
  }
}

/**
 * Airy function Bi(x) using series expansion for small |x| and asymptotic form for large |x|.
 * Bi(x) is the solution to y'' - xy = 0 that grows exponentially for large positive x.
 * Used by the triangular potential solver for finite barrier boundary conditions.
 */
export function airyBi(x: number): number {
  const ABS_X_THRESHOLD = 5.0;

  if (Math.abs(x) < ABS_X_THRESHOLD) {
    // Series expansion for small |x|
    // Bi(x) = c3 * (1 + x^3/(2*3) + x^6/(2*3*5*6) + ...) + c4 * (x + x^4/(3*4) + x^7/(3*4*6*7) + ...)
    // where c3 = 1/(3^(1/6)*Γ(2/3)), c4 = 3^(1/6)/Γ(1/3)
    const c3 = 0.6149266274460007; // 1/(3^(1/6)*Γ(2/3))
    const c4 = 0.4482883573538264; // 3^(1/6)/Γ(1/3)

    let term1 = 1.0;
    let sum1 = 1.0;
    for (let k = 1; k <= 30; k++) {
      term1 *= (x * x * x) / ((3 * k - 1) * (3 * k));
      sum1 += term1;
      if (Math.abs(term1) < 1e-15) break;
    }

    let term2 = x;
    let sum2 = x;
    for (let k = 1; k <= 30; k++) {
      term2 *= (x * x * x) / ((3 * k) * (3 * k + 1));
      sum2 += term2;
      if (Math.abs(term2) < 1e-15) break;
    }

    return c3 * sum1 + c4 * sum2;
  } else if (x > 0) {
    // Asymptotic expansion for large positive x with higher-order corrections
    // Bi(x) ≈ (1/√π) * x^(-1/4) * exp(ζ) * (1 + c₁/ζ + c₂/ζ² + ...)
    // where ζ = (2/3) * x^(3/2)
    const zeta = (2.0 / 3.0) * Math.pow(x, 1.5);
    const factor = 1.0 / Math.sqrt(Math.PI) * Math.pow(x, -0.25);

    // Asymptotic series coefficients (same as Ai but with + signs)
    // First few: c₁ = 5/72, c₂ = 385/10368
    const invZeta = 1.0 / zeta;
    const correction = 1.0 + (5.0 / 72.0) * invZeta + (385.0 / 10368.0) * invZeta * invZeta;

    return factor * Math.exp(zeta) * correction;
  } else {
    // Asymptotic expansion for large negative x
    // Bi(x) ≈ (1/√π) * |x|^(-1/4) * cos(ζ + π/4)
    // where ζ = (2/3) * |x|^(3/2)
    const absX = Math.abs(x);
    const zeta = (2.0 / 3.0) * Math.pow(absX, 1.5);
    const factor = 1.0 / Math.sqrt(Math.PI) * Math.pow(absX, -0.25);
    return factor * Math.cos(zeta + Math.PI / 4);
  }
}

/**
 * Derivative of Airy function Ai'(x) using numerical differentiation.
 * Used by the triangular potential solver for boundary condition matching.
 */
export function airyAiPrime(x: number): number {
  const h = 1e-6;
  return (airyAi(x + h) - airyAi(x - h)) / (2 * h);
}

/**
 * Derivative of Airy function Bi'(x) using numerical differentiation.
 * Used by the triangular potential solver for boundary condition matching.
 */
export function airyBiPrime(x: number): number {
  const h = 1e-6;
  return (airyBi(x + h) - airyBi(x - h)) / (2 * h);
}
