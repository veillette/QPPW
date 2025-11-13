/**
 * BaseScreenView is an abstract base class for all screen views in the QPPW simulation.
 * It provides common functionality including standard layout for quantum well simulations.
 */

import { ScreenView, ScreenViewOptions } from "scenerystack/sim";
import { ResetAllButton } from "scenerystack/scenery-phet";
import { Node } from "scenerystack/scenery";
import { OneWellModel } from "../../one-well/model/OneWellModel.js";
import { EnergyChartNode } from "./EnergyChartNode.js";
import { WaveFunctionChartNode } from "./WaveFunctionChartNode.js";
import { ControlPanelNode } from "./ControlPanelNode.js";
import { SimulationControlBar } from "./SimulationControlBar.js";
import { BaseModel } from "../model/BaseModel.js";

export abstract class BaseScreenView extends ScreenView {
  protected readonly resetAllButton: ResetAllButton;
  protected readonly model: BaseModel | OneWellModel;

  // Common components (may be undefined for screens that don't use them)
  protected energyChart?: EnergyChartNode;
  protected waveFunctionChart?: WaveFunctionChartNode;
  protected controlPanel?: ControlPanelNode;
  protected simulationControlBar?: SimulationControlBar;
  protected chartsContainer?: Node;
  protected listBoxParent?: Node;

  protected constructor(
    model: BaseModel | OneWellModel,
    resetCallback: () => void,
    options?: ScreenViewOptions
  ) {
    super(options);

    this.model = model;

    // Create the reset all button in the bottom-right corner
    this.resetAllButton = new ResetAllButton({
      listener: resetCallback,
      right: this.layoutBounds.maxX - 10,
      bottom: this.layoutBounds.maxY - 10,
    });
    this.addChild(this.resetAllButton);
  }

  /**
   * Creates the standard quantum well layout with charts, control panel, and simulation controls.
   * This should be called by subclasses that use the standard layout.
   * @param model - The OneWellModel instance
   */
  protected createStandardLayout(model: OneWellModel): void {

    // Calculate layout dimensions
    const screenWidth = this.layoutBounds.width;
    const screenHeight = this.layoutBounds.height;
    const margin = 20;

    // Fixed chart dimensions (both charts share the same width and have fixed height)
    const chartsWidth = 600; // Fixed width for consistency
    const chartHeight = 300; // Fixed height for both charts

    // Create the energy chart (top plot)
    this.energyChart = new EnergyChartNode(model, {
      width: chartsWidth,
      height: chartHeight,
    });

    // Create the wave function chart (bottom plot)
    this.waveFunctionChart = new WaveFunctionChartNode(model, {
      width: chartsWidth,
      height: chartHeight,
    });

    // Position charts with fixed layout, ensuring x-axes align horizontally
    this.energyChart.left = margin;
    this.energyChart.top = margin;

    this.waveFunctionChart.left = margin; // Same left position to align x-axes
    this.waveFunctionChart.top = margin + chartHeight + margin;

    // Create container for charts
    this.chartsContainer = new Node({
      children: [this.energyChart, this.waveFunctionChart],
    });

    // Create control panel (needs a parent node for ComboBox listbox)
    this.listBoxParent = new Node();
    this.controlPanel = new ControlPanelNode(model, this.listBoxParent);
    this.controlPanel.left = chartsWidth + margin * 2;
    this.controlPanel.top = margin;

    // Create simulation control bar (footer)
    this.simulationControlBar = new SimulationControlBar(model, {
      width: screenWidth,
    });
    this.simulationControlBar.left = 0;
    this.simulationControlBar.bottom = screenHeight;

    // Add all components to the view
    this.addChild(this.chartsContainer);
    this.addChild(this.controlPanel);
    this.addChild(this.simulationControlBar);
    this.addChild(this.listBoxParent); // ListBox parent must be added last for proper z-ordering
  }

  /**
   * Creates the content for the info dialog.
   * Subclasses must implement this method to provide screen-specific information.
   */
  public abstract createInfoDialogContent(): Node;

  /**
   * Creates the screen summary content for accessibility.
   * Subclasses must implement this method to provide screen-specific summary.
   */
  public abstract createScreenSummaryContent(): Node;

  /**
   * Resets the screen view to its initial state.
   * Subclasses should override this method to add screen-specific reset logic.
   */
  public reset(): void {
    // Base implementation - subclasses should call super.reset() and add their own logic
  }

  /**
   * Steps the screen view forward in time.
   * Subclasses should override this method to add screen-specific step logic.
   * @param dt - The time step in seconds
   */
  // eslint-disable-next-line @typescript-eslint/no-unused-vars
  public step(_dt: number): void {
    // Base implementation - subclasses should override
  }
}
