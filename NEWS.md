
# CellChat 1.1.0 (2020-04-26)
## New added functions
* Add `subsetCellChat` to create a object using a portion of cells
* Add `computeAveExpr` to compute the average expression per cell group
* Add `sketchData` to downsample the single cell data for faster calculation
* Add `netAnalysis_signalingChanges_scatter` to identify the signaling changes associated with one cell group

## Important changes
* `computeCommunProb` now supports a faster calculation and displays a ProgressBar
* `computeCommunProbPathway` now returns the significant pathways that are ordered based on the total communication probabilities


# CellChat 1.0.0 (2021-02-17)

* CellChat paper is now officially published [(Jin et al., Nature Communications, 2021)](https://www.nature.com/articles/s41467-021-21246-9). Compared to the preprint, we have now experimentally validated CellChat's predictions on embryonic skin using RNAscope technique, applied CellChat to a human diseased skin dataset and updated many others. 
* We have now developed a [standalone CellChat Shiny App](https://github.com/sqjin/CellChatShiny) for interactive exploration of the cell-cell communication analyzed by CellChat. Want to share your results with your collaborators like biologists for further exploration? Try it out! 


# CellChat 0.5.0 (2021-01-05)

* Slight changes of CellChat object (Please update your previously calculated CellChat object via `updateCellChat()`)
* Enhanced documentation of functions and tutorials (use `help()` to check the documentation, e.g., `help(CellChat)`)
* New features for comparison analysis of multiple datasets
* Support for creating a new CellChat object from Seurat V3 or SingleCellExperiment object


