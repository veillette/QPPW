/**
 * IntroScreenView is the view for the Intro screen.
 * It displays quantum wells without play/pause controls or superposition options.
 */

import { ScreenViewOptions } from "scenerystack/sim";
import { Node, Text, VBox, RichText } from "scenerystack/scenery";
import { PhetFont, ResetAllButton } from "scenerystack/scenery-phet";
import { ScreenView } from "scenerystack/sim";
import { IntroModel } from "../model/IntroModel.js";
import { EnergyChartNode } from "../../common/view/EnergyChartNode.js";
import { WaveFunctionChartNode } from "../../common/view/WaveFunctionChartNode.js";
import { IntroControlPanelNode } from "./IntroControlPanelNode.js";
import QPPWColors from "../../QPPWColors.js";
import stringManager from "../../i18n/StringManager.js";

export class IntroScreenView extends ScreenView {
    private readonly model: IntroModel;
    private readonly resetAllButton: ResetAllButton;
    private energyChart: EnergyChartNode;
    private waveFunctionChart: WaveFunctionChartNode;
    private controlPanel: IntroControlPanelNode;
    private listBoxParent: Node;

    public constructor(model: IntroModel, options?: ScreenViewOptions) {
        super(options);

        this.model = model;

        // Calculate layout dimensions
        const margin = 10;

        // Fixed chart dimensions
        const chartsWidth = 600;
        const energyChartHeight = 300;
        const waveFunctionChartHeight = 180;

        // Create the energy chart (top plot)
        // Cast to OneWellModel since IntroModel has all necessary properties
        this.energyChart = new EnergyChartNode(model as unknown as import("../../one-well/model/OneWellModel.js").OneWellModel, {
            width: chartsWidth,
            height: energyChartHeight,
        });

        // Create the wave function chart (bottom plot)
        this.waveFunctionChart = new WaveFunctionChartNode(model as unknown as import("../../one-well/model/OneWellModel.js").OneWellModel, {
            width: chartsWidth,
            height: waveFunctionChartHeight,
        });

        // Position charts
        this.energyChart.left = margin;
        this.energyChart.top = 10 - 200;

        this.waveFunctionChart.left = margin;
        this.waveFunctionChart.top = margin + energyChartHeight + 30;


        // Create listbox parent node for ComboBox popups
        this.listBoxParent = new Node();

        // Create control panel (simplified for intro)
        this.controlPanel = new IntroControlPanelNode(
            model,
            this.listBoxParent,
        );
        this.controlPanel.left = chartsWidth + margin * 2;
        this.controlPanel.top = margin;

        // Create the reset all button in the bottom-right corner
        this.resetAllButton = new ResetAllButton({
            listener: () => {
                model.resetAll();
                this.reset();
            },
            right: this.layoutBounds.maxX - 10,
            bottom: this.layoutBounds.maxY - 10,
        });

        // Add all components to the view
        this.addChild(this.energyChart);
        this.addChild(this.waveFunctionChart);
        this.addChild(this.controlPanel);
        this.addChild(this.resetAllButton);
        this.addChild(this.listBoxParent); // ListBox parent must be added last for proper z-ordering
    }

    /**
     * Creates the content for the info dialog.
     */
    public createInfoDialogContent(): Node {
        const titleText = new Text(stringManager.introStringProperty, {
            font: new PhetFont({ size: 18, weight: "bold" }),
            fill: QPPWColors.textFillProperty,
        });

        const descriptionText = new RichText(
            stringManager.introDescriptionStringProperty,
            {
                font: new PhetFont(14),
                fill: QPPWColors.textFillProperty,
                maxWidth: 500,
            },
        );

        const keyConceptsTitle = new Text(
            stringManager.keyConceptsTitleStringProperty,
            {
                font: new PhetFont({ size: 14, weight: "bold" }),
                fill: QPPWColors.textFillProperty,
            },
        );

        const keyConceptsList = new RichText(
            stringManager.introKeyConceptsStringProperty,
            {
                font: new PhetFont(13),
                fill: QPPWColors.textFillProperty,
                maxWidth: 500,
            },
        );

        const interactionTitle = new Text(
            stringManager.interactionsTitleStringProperty,
            {
                font: new PhetFont({ size: 14, weight: "bold" }),
                fill: QPPWColors.textFillProperty,
            },
        );

        const interactionsList = new RichText(
            stringManager.introInteractionsStringProperty,
            {
                font: new PhetFont(13),
                fill: QPPWColors.textFillProperty,
                maxWidth: 500,
            },
        );

        return new VBox({
            spacing: 12,
            align: "left",
            children: [
                titleText,
                descriptionText,
                keyConceptsTitle,
                keyConceptsList,
                interactionTitle,
                interactionsList,
            ],
        });
    }

    /**
     * Creates the screen summary content for accessibility.
     */
    public createScreenSummaryContent(): Node {
        const summaryText = new RichText(
            stringManager.introEducationalContentStringProperty,
            {
                font: new PhetFont(13),
                fill: QPPWColors.textFillProperty,
                maxWidth: 600,
            },
        );

        return new VBox({
            spacing: 10,
            align: "left",
            children: [summaryText],
        });
    }

    /**
     * Resets the screen view to its initial state.
     */
    public reset(): void {
        // Update charts to reflect reset model state
        if (this.energyChart) {
            this.energyChart.update();
        }
        if (this.waveFunctionChart) {
            this.waveFunctionChart.update();
        }
    }

    /**
     * Steps the screen view forward in time.
     * Note: Intro screen doesn't have play/pause controls, but model still steps
     * @param dt - The time step in seconds
     */
    public override step(dt: number): void {
        super.step(dt);
        if (this.model) {
            this.model.step(dt);
        }
    }
}
