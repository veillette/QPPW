/**
 * BaseScreenView - Abstract base class for Quantum Bound States (QPPW) screen views
 * Provides common functionality across different simulation screens using the Template Method pattern
 */

import { ScreenView, ScreenViewOptions } from "scenerystack/sim";
import { Node, Rectangle, VBox } from "scenerystack/scenery";
import { ResetAllButton, TimeControlNode } from "scenerystack/scenery-phet";
import { BooleanProperty, Property } from "scenerystack/axon";
import { Panel } from "scenerystack/sun";
import QPPWColors from "../../QPPWColors.js";
import QPPWPreferences from "../../QPPWPreferences.js";

export interface BaseScreenViewOptions extends ScreenViewOptions {
  // Add any common options here
}

/**
 * Abstract base class for screen views in the Quantum Bound States simulation
 */
export abstract class BaseScreenView extends ScreenView {
  // Visualization properties
  protected readonly showGridProperty: Property<boolean>;
  protected readonly showEnergyLevelsProperty: Property<boolean>;
  protected readonly showWaveFunctionProperty: Property<boolean>;
  protected readonly showProbabilityDensityProperty: Property<boolean>;
  protected readonly showMeasuringTapeProperty: Property<boolean>;

  // Time control
  protected readonly isPlayingProperty: Property<boolean>;
  protected readonly timeSpeedProperty: Property<number>;

  // UI Nodes
  protected readonly resetAllButton: ResetAllButton;
  protected readonly timeControlNode?: TimeControlNode;

  protected constructor(options?: BaseScreenViewOptions) {
    super(options);

    // Initialize visualization properties
    this.showGridProperty = new BooleanProperty(false);
    this.showEnergyLevelsProperty = new BooleanProperty(true);
    this.showWaveFunctionProperty = new BooleanProperty(true);
    this.showProbabilityDensityProperty = new BooleanProperty(false);
    this.showMeasuringTapeProperty = new BooleanProperty(false);

    // Initialize time control properties
    this.isPlayingProperty = new BooleanProperty(false);
    this.timeSpeedProperty = new Property<number>(1);

    // Create reset all button
    this.resetAllButton = new ResetAllButton({
      listener: () => {
        this.reset();
      },
      right: this.layoutBounds.maxX - 10,
      bottom: this.layoutBounds.maxY - 10,
    });
    this.addChild(this.resetAllButton);

    // Setup Page Visibility API for auto-pause
    this.setupPageVisibilityHandling();
  }

  /**
   * Setup automatic pause when tab is hidden
   */
  private setupPageVisibilityHandling(): void {
    if (typeof document !== "undefined" && document.addEventListener) {
      document.addEventListener("visibilitychange", () => {
        if (document.hidden && QPPWPreferences.autoPauseWhenTabHiddenProperty.value) {
          this.isPlayingProperty.value = false;
        }
      });
    }
  }

  /**
   * Setup the scene grid
   */
  protected setupGrid(): Node {
    const gridNode = new Rectangle(0, 0, 100, 100, {
      stroke: QPPWColors.sceneGridColorProperty,
      lineWidth: 1,
    });

    this.showGridProperty.link((showGrid) => {
      gridNode.visible = showGrid;
    });

    return gridNode;
  }

  /**
   * Create a standard control panel with consistent styling
   */
  protected createControlPanel(content: Node, options?: { xMargin?: number; yMargin?: number }): Panel {
    return new Panel(content, {
      fill: QPPWColors.controlPanelBackgroundColorProperty,
      stroke: QPPWColors.controlPanelStrokeColorProperty,
      cornerRadius: 5,
      xMargin: options?.xMargin ?? 10,
      yMargin: options?.yMargin ?? 10,
      align: "left",
    });
  }

  /**
   * Create the visualization controls panel
   */
  protected createVisualizationPanel(): Panel {
    // Subclasses can override this to add specific visualization controls
    const content = new VBox({
      spacing: 10,
      align: "left",
    });

    return this.createControlPanel(content);
  }

  /**
   * Create information dialog content - must be implemented by subclasses
   */
  protected abstract createInfoDialogContent(): Node;

  /**
   * Create screen summary content for accessibility - must be implemented by subclasses
   */
  protected abstract createScreenSummaryContent(): Node;

  /**
   * Reset the screen view
   */
  public reset(): void {
    this.showGridProperty.reset();
    this.showEnergyLevelsProperty.reset();
    this.showWaveFunctionProperty.reset();
    this.showProbabilityDensityProperty.reset();
    this.showMeasuringTapeProperty.reset();
    this.isPlayingProperty.reset();
    this.timeSpeedProperty.reset();
  }

  /**
   * Step function called every frame
   */
  public step(dt: number): void {
    // Subclasses should override this to add specific stepping behavior
    if (this.isPlayingProperty.value) {
      const adjustedDt = dt * this.timeSpeedProperty.value;
      this.stepSimulation(adjustedDt);
    }
  }

  /**
   * Step the simulation - subclasses should implement this
   */
  protected stepSimulation(dt: number): void {
    // To be implemented by subclasses
  }

  /**
   * Layout function to position elements when screen size changes
   */
  protected layout(viewBounds: { width: number; height: number }): void {
    // Subclasses can override this to add specific layout behavior
  }
}
