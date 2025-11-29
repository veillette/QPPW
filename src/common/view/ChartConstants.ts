/**
 * ChartConstants defines shared constants for chart layout and appearance.
 * These constants ensure vertical alignment of the x-axes between EnergyChartNode and WaveFunctionChartNode.
 */

/**
 * Chart margins that MUST be identical across all charts to ensure vertical alignment.
 * The left margin determines where the plot area begins, and matching this value ensures
 * that the x-axes of the energy chart and wavefunction chart are perfectly aligned.
 */
export const CHART_MARGINS = {
  left: 60, // Space for y-axis labels (must match across charts for alignment)
  right: 20, // Right padding
  top: 40, // Top padding for energy chart
  bottom: 50, // Bottom padding for energy chart (includes x-axis labels)
} as const;

/**
 * Chart margins for the wavefunction chart.
 * Left and right margins MUST match CHART_MARGINS to ensure x-axis alignment.
 */
export const WAVEFUNCTION_CHART_MARGINS = {
  left: CHART_MARGINS.left, // MUST match CHART_MARGINS.left for alignment
  right: CHART_MARGINS.right, // MUST match CHART_MARGINS.right for consistent width
  top: 10, // Reduced top padding for wavefunction chart
  bottom: 40, // Bottom padding (includes x-axis labels)
} as const;

/**
 * X-axis range in nanometers.
 * The chart x-axis extends from -X_AXIS_RANGE_NM to +X_AXIS_RANGE_NM.
 * This value is shared between EnergyChartNode and WaveFunctionChartNode.
 */
export const X_AXIS_RANGE_NM = 4;

/**
 * Y-axis label offset from the left edge of the chart.
 * This determines how far left the rotated y-axis label text appears.
 * MUST be identical across all charts to ensure visual alignment.
 */
export const Y_AXIS_LABEL_OFFSET = 40;
