# CellChat 1.1.2 (2021-07-10)
## New added functions
* Add `netAnalysis_diff_signalingRole_scatter` for 2D visualization of differential signaling roles of each cell group when comparing mutiple datasets.

## Updated functions with minor changes
`netVisual`, `netVisual_aggregate`, `netVisual_individual`, `netVisual_hierarchy1`, `netVisual_hierarchy2`,`netVisual_circle`,`createCellChat`,`netAnalysis_computeCentrality` ,`netAnalysis_signalingRole_scatter`, `netAnalysis_signalingChanges_scatter`

## Important changes
* Add `thresh = 0.05` in `netAnalysis_computeCentrality` to only consider the significant interactions. This will slightly change (very likely quantitative instead of qualitative change) the results computed by previous version of CellChat. 
* In the updated `netAnalysis_computeCentrality`, we now also compute unweighted outdegree (i.e., the total number of outgoing links) and indegree (i.e., the total number of incoming links). 
* `netAnalysis_signalingRole_scatter` and `netAnalysis_signalingChanges_scatter` now support the comparison of the total number of outgoing and incoming links. 
* Change the default setting for visualizing cell-cell communication network: 1) using `circle` plot instead of `hierarchy`; 2) using the same node size instead of different size (setting `vertex.weight = NULL` will give different size as in previous version of CellChat).  


# CellChat 1.1.1 (2021-06-19)
## New added functions
* Add `updateClusterLabels` to update cell cluster labels without re-running the time-consuming function `computeCommunProb`.
## Updated functions with minor changes
`netAnalysis_signalingRole_scatter`, `setIdent`, `netVisual_diffInteraction`, `sketchData`, `identifyOverExpressedGenes`,`netEmbedding`,`runUMAP`

# CellChat 1.1.0 (2021-04-26)
## New added functions
* Add `subsetCellChat` to create a object using a portion of cells
* Add `computeAveExpr` to compute the average expression per cell group
* Add `sketchData` to downsample the single cell data for faster calculation
* Add `netAnalysis_signalingChanges_scatter` to identify the signaling changes associated with one cell group

## Important changes
* `computeCommunProb` now supports a faster calculation and displays a ProgressBar
* `computeCommunProbPathway` now returns the significant pathways that are ordered based on the total communication probabilities

## Publish a new release
* The first release was published as version 1.0.0 before updating to version 1.1.0 

# CellChat 1.0.0 (2021-02-17)

* CellChat paper is now officially published [(Jin et al., Nature Communications, 2021)](https://www.nature.com/articles/s41467-021-21246-9). Compared to the preprint, we have now experimentally validated CellChat's predictions on embryonic skin using RNAscope technique, applied CellChat to a human diseased skin dataset and updated many others. 
* We have now developed a [standalone CellChat Shiny App](https://github.com/sqjin/CellChatShiny) for interactive exploration of the cell-cell communication analyzed by CellChat. Want to share your results with your collaborators like biologists for further exploration? Try it out! 


# CellChat 0.5.0 (2021-01-05)

* Slight changes of CellChat object (Please update your previously calculated CellChat object via `updateCellChat()`)
* Enhanced documentation of functions and tutorials (use `help()` to check the documentation, e.g., `help(CellChat)`)
* New features for comparison analysis of multiple datasets
* Support for creating a new CellChat object from Seurat V3 or SingleCellExperiment object


