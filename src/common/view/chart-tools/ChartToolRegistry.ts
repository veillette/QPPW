/**
 * ChartToolRegistry manages a collection of chart tools and provides
 * centralized control for updating them based on chart state changes.
 *
 * This class follows the Registry pattern to decouple the chart from
 * individual tool implementations and makes it easier to add, remove,
 * or modify tools without changing the main chart code.
 */

import type { AreaMeasurementTool } from "./AreaMeasurementTool.js";
import type { CurvatureTool } from "./CurvatureTool.js";
import type { DerivativeTool } from "./DerivativeTool.js";
import type { BoundStateResult } from "../../model/PotentialFunction.js";

/**
 * Context object passed to tools during updates.
 * Contains all the information tools need to update themselves.
 */
export interface ToolUpdateContext {
  displayMode: string;
  boundStates: BoundStateResult | null;
  selectedIndex: number;
}

/**
 * Base interface for chart tools that can be registered.
 * All tools must implement this interface to work with the registry.
 */
export interface ChartTool {
  /**
   * Property that controls whether the tool is visible/enabled.
   */
  readonly showProperty: { value: boolean };

  /**
   * Update the tool based on the current chart context.
   */
  update(displayMode: string): void;
}

/**
 * Registry for managing chart tools.
 * Provides centralized control for tool lifecycle and updates.
 */
export class ChartToolRegistry {
  private readonly tools: Map<string, ChartTool> = new Map();

  /**
   * Register a tool with the registry.
   * @param name - Unique identifier for the tool
   * @param tool - The tool instance to register
   */
  public registerTool(name: string, tool: ChartTool): void {
    if (this.tools.has(name)) {
      console.warn(
        `[ChartToolRegistry] Tool "${name}" is already registered. Overwriting.`,
      );
    }
    this.tools.set(name, tool);
  }

  /**
   * Unregister a tool from the registry.
   * @param name - Identifier of the tool to remove
   * @returns true if the tool was removed, false if it wasn't found
   */
  public unregisterTool(name: string): boolean {
    return this.tools.delete(name);
  }

  /**
   * Get a registered tool by name.
   * @param name - Identifier of the tool
   * @returns The tool instance, or undefined if not found
   */
  public getTool(name: string): ChartTool | undefined {
    return this.tools.get(name);
  }

  /**
   * Check if a tool is registered.
   * @param name - Identifier of the tool
   * @returns true if the tool is registered
   */
  public hasTool(name: string): boolean {
    return this.tools.has(name);
  }

  /**
   * Update all enabled tools with the given context.
   * Only tools that have their showProperty set to true will be updated.
   *
   * @param context - The update context containing display mode and chart state
   */
  public updateAllTools(context: ToolUpdateContext): void {
    this.tools.forEach((tool, name) => {
      if (tool.showProperty.value) {
        try {
          tool.update(context.displayMode);
        } catch (error) {
          console.error(
            `[ChartToolRegistry] Error updating tool "${name}":`,
            error,
          );
        }
      }
    });
  }

  /**
   * Update a specific tool by name.
   * @param name - Identifier of the tool to update
   * @param context - The update context
   * @returns true if the tool was found and updated, false otherwise
   */
  public updateTool(name: string, context: ToolUpdateContext): boolean {
    const tool = this.tools.get(name);
    if (tool && tool.showProperty.value) {
      try {
        tool.update(context.displayMode);
        return true;
      } catch (error) {
        console.error(
          `[ChartToolRegistry] Error updating tool "${name}":`,
          error,
        );
      }
    }
    return false;
  }

  /**
   * Get the names of all registered tools.
   * @returns Array of tool names
   */
  public getToolNames(): string[] {
    return Array.from(this.tools.keys());
  }

  /**
   * Get the count of registered tools.
   * @returns Number of registered tools
   */
  public getToolCount(): number {
    return this.tools.size;
  }

  /**
   * Clear all registered tools.
   */
  public clear(): void {
    this.tools.clear();
  }
}
