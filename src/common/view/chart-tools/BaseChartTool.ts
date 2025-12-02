/**
 * BaseChartTool - Abstract base class for interactive chart tools.
 * Provides common functionality for drag listeners, accessibility, coordinate transformations,
 * and visibility management.
 */

import {
  Node,
  DragListener,
  KeyboardDragListener,
  type NodeOptions,
} from "scenerystack/scenery";
import {
  NumberProperty,
  BooleanProperty,
  DerivedProperty,
} from "scenerystack/axon";
import type { ScreenModel } from "../../model/ScreenModels.js";
import { UtteranceQueue, AriaLiveAnnouncer } from "scenerystack/utterance-queue";

// Create a global utteranceQueue instance for accessibility announcements
// Using AriaLiveAnnouncer for screen reader support via aria-live regions
export const utteranceQueue = new UtteranceQueue(new AriaLiveAnnouncer());

export type BaseChartToolOptions = {
  chartMargins: { left: number; right: number; top: number; bottom: number };
  plotHeight: number;
  plotWidth?: number;
  xMinProperty: NumberProperty;
  xMaxProperty: NumberProperty;
  yMinProperty?: NumberProperty;
  yMaxProperty?: NumberProperty;
  dataToViewX: (x: number) => number;
  dataToViewY: (y: number) => number;
  viewToDataX: (x: number) => number;
  parentNode: Node; // Reference to parent chart node for coordinate conversions
};

export type MarkerConfig = {
  markerNode: Node; // The draggable node (marker line or handle)
  positionProperty: NumberProperty; // Position in data coordinates
  accessibleName: string; // Name for screen readers
  labelContent: string; // Label for PDOM
  helpText: string; // Help text for keyboard users
  constraintFn?: (newX: number) => number; // Optional function to constrain position
  onDragEnd?: (position: number) => void; // Optional callback on drag end
};

export abstract class BaseChartTool extends Node {
  protected readonly model: ScreenModel;
  protected readonly options: BaseChartToolOptions;
  protected readonly container: Node;
  public readonly showProperty: BooleanProperty;

  // Abstract methods that subclasses must implement
  protected abstract setupVisualElements(): void;
  public abstract update(displayMode: string): void;

  constructor(
    model: ScreenModel,
    getEffectiveDisplayMode: () => string,
    options: BaseChartToolOptions,
    toolName: string,
    toolDescription: string,
    nodeOptions?: NodeOptions,
  ) {
    // Initialize showProperty before calling super()
    const showPropertyInternal = new BooleanProperty(false);

    super({
      // pdom - container for the tool
      tagName: "div",
      labelTagName: "h3",
      labelContent: toolName,
      descriptionTagName: "p",
      descriptionContent: new DerivedProperty(
        [showPropertyInternal],
        (enabled) =>
          enabled
            ? toolDescription
            : `${toolName} disabled.`,
      ),
      ...nodeOptions,
    });

    this.model = model;
    this.options = options;
    this.showProperty = showPropertyInternal;

    // Create container
    this.container = new Node({
      visible: false,
    });
    this.addChild(this.container);

    // Subclasses setup their visual elements
    this.setupVisualElements();

    // Store display mode getter
    const getDisplayMode = getEffectiveDisplayMode;

    // Link visibility to property
    this.showProperty.link((show: boolean) => {
      this.container.visible = show;
      if (show) {
        this.update(getDisplayMode());
      }
    });
  }

  /**
   * Setup drag listeners for a marker (both mouse and keyboard).
   * This method handles the common pattern of dragging a marker with coordinate conversion,
   * clamping, and accessibility announcements.
   */
  protected setupMarkerDragListener(config: MarkerConfig): void {
    const { parentNode, viewToDataX, xMinProperty, xMaxProperty } =
      this.options;
    const {
      markerNode,
      positionProperty,
      accessibleName,
      labelContent,
      helpText,
      constraintFn,
      onDragEnd,
    } = config;

    // Setup PDOM attributes for accessibility
    markerNode.tagName = "div";
    markerNode.ariaRole = "slider";
    markerNode.focusable = true;
    markerNode.accessibleName = accessibleName;
    markerNode.labelContent = labelContent;
    markerNode.pdomAttributes = [
      { attribute: "aria-valuemin", value: xMinProperty.value },
      { attribute: "aria-valuemax", value: xMaxProperty.value },
      {
        attribute: "aria-valuenow",
        value: positionProperty.value.toFixed(2),
      },
      {
        attribute: "aria-valuetext",
        value: `Position: ${positionProperty.value.toFixed(2)} nanometers`,
      },
    ];
    markerNode.accessibleHelpText = helpText;

    // Update aria-valuetext when position changes
    positionProperty.link((position) => {
      markerNode.setPDOMAttribute("aria-valuenow", position.toFixed(2));
      markerNode.setPDOMAttribute(
        "aria-valuetext",
        `Position: ${position.toFixed(2)} nanometers`,
      );
    });

    // Mouse drag listener
    const dragListener = new DragListener({
      drag: (event) => {
        // Convert to parent coordinate system
        const parentPoint = parentNode.globalToLocalPoint(event.pointer.point);
        let dataX = viewToDataX(parentPoint.x);

        // Apply constraint function if provided, otherwise just clamp to bounds
        if (constraintFn) {
          dataX = constraintFn(dataX);
        } else {
          dataX = this.clampToRange(
            dataX,
            xMinProperty.value,
            xMaxProperty.value,
          );
        }

        positionProperty.value = dataX;
      },
    });

    markerNode.addInputListener(dragListener);

    // Keyboard drag listener
    const keyboardDragListener = new KeyboardDragListener({
      drag: (_event, listener) => {
        let newX = positionProperty.value + listener.modelDelta.x * 0.1;

        // Apply constraint function if provided, otherwise just clamp to bounds
        if (constraintFn) {
          newX = constraintFn(newX);
        } else {
          newX = this.clampToRange(newX, xMinProperty.value, xMaxProperty.value);
        }

        positionProperty.value = newX;
      },
      dragDelta: 1, // Regular arrow key step
      shiftDragDelta: 0.1, // Fine control with shift
      end: () => {
        // Call optional drag end callback
        if (onDragEnd) {
          onDragEnd(positionProperty.value);
        }
      },
    });

    markerNode.addInputListener(keyboardDragListener);
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
   * Coordinate transformation: convert view x-coordinate to data x-coordinate
   */
  protected viewToDataX(x: number): number {
    return this.options.viewToDataX(x);
  }

  /**
   * Clamps a value to a specified range
   */
  protected clampToRange(value: number, min: number, max: number): number {
    return Math.max(min, Math.min(max, value));
  }

  /**
   * Get the y-coordinate for the top of the chart area
   */
  protected getChartTop(): number {
    return this.options.chartMargins.top;
  }

  /**
   * Get the y-coordinate for the bottom of the chart area
   */
  protected getChartBottom(): number {
    return this.options.chartMargins.top + this.options.plotHeight;
  }
}
