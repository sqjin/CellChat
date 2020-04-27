# CellChat: Inference and analysis of cell-cell communication

<p align="center">
  <img width="250"  src="https://github.com/sqjin/CellChat/blob/master/CellChat_Logo.png">
</p>


## Web-based “CellChat Explorer” 

We build a user-friendly web-based “[CellChat Explorer](http://www.cellchat.org/)” that contains two major components:
- **Ligand-Receptor Interaction Explorer** that allows easy exploration of our ligand-receptor interaction database
- **Cell-Cell Communication Atlas Explorer** that allows easy exploration of the cell-cell communication for any given scRNA-seq dataset that has been processed by our R toolkit CellChat. 

## Capabilities
In addition to infer the intercellular communication from any given scRNA-seq data, CellChat provides functionality for further data exploration, analysis, and visualization. 

- It is able to analyze cell-cell communication for continuous states along cellular development trajectories.
- It can quantitatively characterize and compare the inferred cell-cell communication networks using an integrated approach by combining social network analysis, pattern recognition, and manifold learning approaches.
- It provides an easy-to-use tool for extracting and visualizing high-order information of the inferred networks. For example, it allows ready prediction of major signaling inputs and outputs for all cell populations and how these populations and signals coordinate together for functions.
- It provides several visualization outputs to facilitate intuitive user-guided data interpretation.


## Installation
CellChat R package can be easily installed from Github using devtools:  

```
devtools::install_github("sqjin/CellChat")
```

## Vignettes
Please check the vignettes directory of the repo.

- [Basic commands tutorial](https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/vignettes/CellChat-vignette.html)
- [Joint analysis of multiple datasets](https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/vignettes/Joint_analysis_of_multiple_datasets.html)
- [Interface with other single-cell analysis toolkits (e.g., Seurat, Scanpy)](https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/vignettes/Interface_with_other_single-cell_analysis_toolkits.html)

## Help
If you have any problems, comments or suggestions, please contact us at Suoqin Jin (suoqin.jin@uci.edu).


