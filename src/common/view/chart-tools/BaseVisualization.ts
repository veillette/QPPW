/**
 * BaseVisualization - Abstract base class for visualization overlays.
 * Provides common functionality for container management, visibility control,
 * and coordinate transformations.
 */

import { Node } from "scenerystack/scenery";

export type BaseVisualizationOptions = {
  dataToViewX: (x: number) => number;
  dataToViewY: (y: number) => number;
};

export abstract class BaseVisualization extends Node {
  protected readonly options: BaseVisualizationOptions;
  protected readonly container: Node;

  constructor(options: BaseVisualizationOptions) {
    super();

    this.options = options;

    // Create container for visual elements
    this.container = new Node({
      visible: false,
    });
    this.addChild(this.container);
  }

  /**
   * Coordinate transformation: convert data x-coordinate to view x-coordinate
   */
  protected dataToViewX(x: number): number {
    return this.options.dataToViewX(x);
  }

  /**
   * Coordinate transformation: convert data y-coordinate to view y-coordinate
   */
  protected dataToViewY(y: number): number {
    return this.options.dataToViewY(y);
  }

  /**
   * Show the visualization
   */
  public show(): void {
    this.container.visible = true;
  }

  /**
   * Hide the visualization
   */
  public hide(): void {
    this.container.visible = false;
  }

  /**
   * Check if the visualization is currently visible
   */
  public isVisible(): boolean {
    return this.container.visible;
  }
}
