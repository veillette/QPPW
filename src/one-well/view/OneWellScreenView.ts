/**
 * OneWellScreenView is the main view for the One Well screen.
 * It displays a single quantum potential well and its wave functions.
 */

import { BaseScreenView } from "../../common/view/BaseScreenView.js";
import { OneWellModel } from "../model/OneWellModel.js";
import { ScreenViewOptions } from "scenerystack/sim";
import { Node, VBox, Text } from "scenerystack/scenery";
import { EnergyChartNode } from "./EnergyChartNode.js";
import { WaveFunctionChartNode } from "./WaveFunctionChartNode.js";
import { ControlPanelNode } from "./ControlPanelNode.js";
import { SimulationControlBar } from "./SimulationControlBar.js";
import QPPWColors from "../../QPPWColors.js";

export class OneWellScreenView extends BaseScreenView {
  private readonly model: OneWellModel;
  private readonly energyChart: EnergyChartNode;
  private readonly waveFunctionChart: WaveFunctionChartNode;
  private readonly controlPanel: ControlPanelNode;
  private readonly simulationControlBar: SimulationControlBar;

  public constructor(model: OneWellModel, options?: ScreenViewOptions) {
    super(
      () => {
        model.reset();
        this.reset();
      },
      options,
    );

    this.model = model;

    // Calculate layout dimensions
    const screenWidth = this.layoutBounds.width;
    const screenHeight = this.layoutBounds.height;
    const margin = 20;
    const controlBarHeight = 80;

    // Left side: 70% width for charts
    const chartsWidth = screenWidth * 0.65;
    const chartHeight = (screenHeight - controlBarHeight - margin * 3) / 2;

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

    // Stack charts vertically
    const chartsStack = new VBox({
      spacing: margin,
      align: "left",
      children: [this.energyChart, this.waveFunctionChart],
      left: margin,
      top: margin,
    });

    // Create control panel (needs a parent node for ComboBox listbox)
    const listBoxParent = new Node();
    this.controlPanel = new ControlPanelNode(
      model,
      () => {
        model.reset();
        this.reset();
      },
      listBoxParent,
    );
    this.controlPanel.left = chartsWidth + margin * 2;
    this.controlPanel.top = margin;

    // Create simulation control bar (footer)
    this.simulationControlBar = new SimulationControlBar(model, {
      width: screenWidth,
    });
    this.simulationControlBar.left = 0;
    this.simulationControlBar.bottom = screenHeight;

    // Add all components to the view
    this.addChild(chartsStack);
    this.addChild(this.controlPanel);
    this.addChild(this.simulationControlBar);
    this.addChild(listBoxParent); // ListBox parent must be added last for proper z-ordering
  }

  /**
   * Creates the content for the info dialog.
   */
  public createInfoDialogContent(): Node {
    const text = new Text(
      "Explore quantum mechanics in a single potential well.\n" +
        "Adjust the well parameters to see how energy levels change.",
      {
        font: "14px sans-serif",
        fill: QPPWColors.textFillProperty,
      },
    );
    return new Node({ children: [text] });
  }

  /**
   * Creates the screen summary content for accessibility.
   */
  public createScreenSummaryContent(): Node {
    const text = new Text(
      "One Well screen shows a single quantum potential well with adjustable parameters.",
      {
        font: "14px sans-serif",
        fill: QPPWColors.textFillProperty,
      },
    );
    return new Node({ children: [text] });
  }

  /**
   * Resets the screen view to its initial state.
   */
  public override reset(): void {
    super.reset();
    // Add screen-specific reset logic here
  }

  /**
   * Steps the screen view forward in time.
   * @param dt - The time step in seconds
   */
  public override step(dt: number): void {
    super.step(dt);
    this.model.step(dt);
  }
}
