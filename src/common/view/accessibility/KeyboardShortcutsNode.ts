/**
 * KeyboardShortcutsNode provides the keyboard help dialog for QPPW screens.
 * It describes available keyboard shortcuts and navigation patterns for the
 * quantum physics simulation.
 *
 * This component is displayed when users press the keyboard shortcuts button
 * and helps them understand how to navigate and control the simulation using
 * only the keyboard.
 */

import { Node, Text, VBox, HBox } from "scenerystack/scenery";
import { PhetFont } from "scenerystack/scenery-phet";
import { ReadOnlyProperty } from "scenerystack";
import stringManager from "../../../i18n/StringManager.js";

/**
 * KeyboardShortcutsNode creates an accessible keyboard shortcuts help dialog.
 * The dialog includes:
 * - Simulation control shortcuts (Space for play/pause, R for reset)
 * - Energy level navigation (Arrow keys, Home/End)
 * - Graph interaction instructions
 * - All text is internationalized via StringManager
 */
export class KeyboardShortcutsNode extends Node {
  public constructor() {
    super();

    const keyboardStrings = stringManager.getKeyboardShortcutsStrings();

    // Title
    const titleText = new Text(keyboardStrings.titleStringProperty, {
      font: new PhetFont({ size: 20, weight: "bold" }),
      maxWidth: 600,
    });

    // Simulation Controls Section
    const simulationControlsHeader = new Text(
      keyboardStrings.simulationControlsStringProperty,
      {
        font: new PhetFont({ size: 16, weight: "bold" }),
        maxWidth: 600,
      },
    );

    const playPauseRow = this.createShortcutRow(
      "Space",
      keyboardStrings.playPauseDescriptionStringProperty,
    );
    const resetRow = this.createShortcutRow(
      "R",
      keyboardStrings.resetDescriptionStringProperty,
    );

    const simulationControlsBox = new VBox({
      align: "left",
      spacing: 8,
      children: [simulationControlsHeader, playPauseRow, resetRow],
    });

    // Energy Level Navigation Section
    const energyLevelHeader = new Text(
      keyboardStrings.energyLevelNavigationStringProperty,
      {
        font: new PhetFont({ size: 16, weight: "bold" }),
        maxWidth: 600,
      },
    );

    const arrowUpRow = this.createShortcutRow(
      "↑ or →",
      keyboardStrings.arrowUpRightDescriptionStringProperty,
    );
    const arrowDownRow = this.createShortcutRow(
      "↓ or ←",
      keyboardStrings.arrowDownLeftDescriptionStringProperty,
    );
    const homeRow = this.createShortcutRow(
      "Home",
      keyboardStrings.homeDescriptionStringProperty,
    );
    const endRow = this.createShortcutRow(
      "End",
      keyboardStrings.endDescriptionStringProperty,
    );

    const energyLevelBox = new VBox({
      align: "left",
      spacing: 8,
      children: [energyLevelHeader, arrowUpRow, arrowDownRow, homeRow, endRow],
    });

    // Graph Interactions Section
    const graphInteractionsHeader = new Text(
      keyboardStrings.graphInteractionsStringProperty,
      {
        font: new PhetFont({ size: 16, weight: "bold" }),
        maxWidth: 600,
      },
    );

    const doubleClickRow = this.createShortcutRow(
      "Double-click",
      keyboardStrings.doubleClickResetZoomStringProperty,
    );
    const mouseWheelRow = this.createShortcutRow(
      "Mouse wheel",
      keyboardStrings.mouseWheelZoomStringProperty,
    );
    const dragRow = this.createShortcutRow(
      "Drag",
      keyboardStrings.dragToPanStringProperty,
    );

    const graphInteractionsBox = new VBox({
      align: "left",
      spacing: 8,
      children: [
        graphInteractionsHeader,
        doubleClickRow,
        mouseWheelRow,
        dragRow,
      ],
    });

    // Main content layout
    const content = new VBox({
      align: "left",
      spacing: 20,
      children: [
        titleText,
        simulationControlsBox,
        energyLevelBox,
        graphInteractionsBox,
      ],
    });

    this.addChild(content);

    // PDOM - Make this accessible
    this.tagName = "div";
    this.labelTagName = "h2";
    this.labelContent = "Keyboard Shortcuts";
  }

  /**
   * Creates a row displaying a keyboard shortcut and its description
   */
  private createShortcutRow(
    key: string,
    descriptionProperty: ReadOnlyProperty<string>,
  ): HBox {
    const keyText = new Text(key, {
      font: new PhetFont({ size: 14, weight: "bold" }),
      maxWidth: 150,
    });

    const descriptionText = new Text(descriptionProperty, {
      font: new PhetFont({ size: 14 }),
      maxWidth: 400,
    });

    return new HBox({
      spacing: 15,
      children: [keyText, descriptionText],
      align: "center",
    });
  }
}
