/**
 * BaseScreenView is an abstract base class for all screen views in the QPPW simulation.
 * It provides common functionality including standard layout for quantum well simulations.
 */

import { ScreenView, ScreenViewOptions } from "scenerystack/sim";
import { ResetAllButton } from "scenerystack/scenery-phet";
import { Node } from "scenerystack/scenery";
import { OneWellModel } from "../../one-well/model/OneWellModel.js";
import { TwoWellsModel } from "../../two-wells/model/TwoWellsModel.js";
import { ManyWellsModel } from "../../many-wells/model/ManyWellsModel.js";
import { EnergyChartNode } from "./EnergyChartNode.js";
import { WaveFunctionChartNode } from "./WaveFunctionChartNode.js";
import { ControlPanelNode, ControlPanelNodeOptions } from "./ControlPanelNode.js";
import { SimulationControlBar } from "./SimulationControlBar.js";
import { BaseModel } from "../model/BaseModel.js";

export abstract class BaseScreenView extends ScreenView {
  protected readonly resetAllButton: ResetAllButton;
  protected readonly model: BaseModel | OneWellModel | TwoWellsModel | ManyWellsModel;

  // Common components (may be undefined for screens that don't use them)
  protected energyChart?: EnergyChartNode;
  protected waveFunctionChart?: WaveFunctionChartNode;
  protected controlPanel?: ControlPanelNode;
  protected simulationControlBar?: SimulationControlBar;
  protected chartsContainer?: Node;
  protected listBoxParent?: Node;

  protected constructor(
    model: BaseModel | OneWellModel | TwoWellsModel | ManyWellsModel,
    options?: ScreenViewOptions
  ) {
    super(options);

    this.model = model;

    // Create the reset all button in the bottom-right corner
    this.resetAllButton = new ResetAllButton({
      listener: () => {
        model.resetAll();
        this.reset();
      },
      right: this.layoutBounds.maxX - 10,
      bottom: this.layoutBounds.maxY - 10,
    });
    this.addChild(this.resetAllButton);
  }

  /**
   * Creates the standard quantum well layout with charts, control panel, and simulation controls.
   * This should be called by subclasses that use the standard layout.
   * @param model - The OneWellModel or TwoWellsModel instance
   * @param controlPanelOptions - Optional configuration for the control panel (e.g., hiding mass slider, filtering potential types)
   */
  protected createStandardLayout(model: OneWellModel | TwoWellsModel, controlPanelOptions?: ControlPanelNodeOptions): void {

    // Calculate layout dimensions
    const screenWidth = this.layoutBounds.width;
    const screenHeight = this.layoutBounds.height;
    const margin = 20;

    // Fixed chart dimensions (both charts share the same width and have fixed height)
    const chartsWidth = 600; // Fixed width for consistency
    const energyChartHeight = 300; // Fixed height for energy chart
    const waveFunctionChartHeight = 180; // Fixed height for wavefunction chart (reduced by 40%)

    // Create the energy chart (top plot)
    this.energyChart = new EnergyChartNode(model, {
      width: chartsWidth,
      height: energyChartHeight,
    });

    // Create the wave function chart (bottom plot)
    this.waveFunctionChart = new WaveFunctionChartNode(model, {
      width: chartsWidth,
      height: waveFunctionChartHeight,
    });

    // Position charts with fixed layout, ensuring x-axes align horizontally
    this.energyChart.left = margin;
    this.energyChart.top = 10; // Reduced from margin to move energy chart upward and avoid overlap

    this.waveFunctionChart.left = margin; // Same left position to align x-axes
    this.waveFunctionChart.top = margin + energyChartHeight + 30; // Added extra spacing between charts

    // Create container for charts
    this.chartsContainer = new Node({
      children: [this.energyChart, this.waveFunctionChart],
    });

    // Create listbox parent node for ComboBox popups (if not already created by subclass)
    if (!this.listBoxParent) {
      this.listBoxParent = new Node();
    }

    // Create control panel with optional configuration
    this.controlPanel = new ControlPanelNode(model, this.listBoxParent, controlPanelOptions);
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
    // Update charts to reflect reset model state
    if (this.energyChart) {
      this.energyChart.update();
    }
    if (this.waveFunctionChart) {
      this.waveFunctionChart.update();
    }
    // Base implementation - subclasses should call super.reset() and add their own logic
  }

  /**
   * Steps the screen view forward in time.
   * Subclasses should override this method to add screen-specific step logic.
   * @param dt - The time step in seconds
   */
   
  public step(_dt: number): void {
    // Base implementation - subclasses should override
  }
}
