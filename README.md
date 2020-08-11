# CellChat: Inference and analysis of cell-cell communication

<p align="center">
  <img width="250"  src="https://github.com/sqjin/CellChat/blob/master/CellChat_Logo.png">
</p>


## Web-based “CellChat Explorer” 

We build a user-friendly web-based “[CellChat Explorer](http://www.cellchat.org/)” that contains two major components:
- **Ligand-Receptor Interaction Explorer** that allows easy exploration of our **novel** ligand-receptor interaction database, a comprehensive recapitulation of known molecular compositions including multimeric complexes and co-factors. *Our database CellChatDB is a manually curated database of literature-supported ligand-receptor interactions in both **human and mouse***.
- **Cell-Cell Communication Atlas Explorer** that allows easy exploration of the cell-cell communication for any given scRNA-seq dataset that has been processed by our R toolkit CellChat. 

## Capabilities
In addition to infer the intercellular communication from any given scRNA-seq data, CellChat provides functionality for further data exploration, analysis, and visualization. 

- It is able to analyze cell-cell communication for continuous states along cellular development trajectories.
- It can quantitatively characterize and compare the inferred cell-cell communication networks using an integrated approach by combining social network analysis, pattern recognition, and manifold learning approaches.
- It provides an easy-to-use tool for extracting and visualizing high-order information of the inferred networks. For example, it allows ready prediction of major signaling inputs and outputs for all cell populations and how these populations and signals coordinate together for functions.
- It provides several visualization outputs to facilitate intuitive user-guided data interpretation.

<p align="center">
  <img width="700"  src="https://github.com/sqjin/CellChat/blob/master/overview_CellChat.png">
</p>


Check out [our preprint on bioRxiv](https://www.biorxiv.org/content/10.1101/2020.07.21.214387v1) for the detailed methods and applications, and a [short video](https://youtu.be/lag9UstpYhk) for introducing the key features of CellChat.


## Installation

CellChat R package can be easily installed from Github using devtools:  

```
devtools::install_github("sqjin/CellChat")
```
You might need to manually intall the dependencies [ComplexHeatmap](https://github.com/jokergoo/ComplexHeatmap) using `devtools::install_github("jokergoo/ComplexHeatmap")` if it is not automatically installed. 

## Vignettes
Please check the vignettes directory of the repo.

- [Basic commands tutorial](https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/vignettes/CellChat-vignette.html)
- [Comparison analysis of multiple datasets](https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/vignettes/Joint_analysis_of_multiple_datasets.html)
- [Walkthrough - CellChat analysis of cell-cell communication in mouse skin wounds](https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/vignettes/walkthrough_wound.html)
- [Example streamline for quick analysis and exploration](https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/vignettes/example_streamline_for_quick_analysis_and_exploration.html)
- [Interface with other single-cell analysis toolkits (e.g., Seurat, Scanpy)](https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/vignettes/Interface_with_other_single-cell_analysis_toolkits.html)


## System Requirements
### Hardware requirements
CellChat package requires only a standard computer with enough RAM to support the in-memory operations.

### Software requirements
This package is supported for macOS, Windows and Linux. The package has been tested on macOS: Mojave (10.14.5) and Windows 10. 

Dependencies of CellChat package are indicated in the Description file, and can be automatically installed when installing CellChat pacakge. 

## Help
If you have any problems, comments or suggestions, please contact us at Suoqin Jin (suoqin.jin@uci.edu).

## How should I cite CellChat?
Suoqin Jin, Christian F. Guerrero-Juarez, Lihua Zhang, Ivan Chang, Peggy Myung, Maksim V. Plikus, Qing Nie. **Inference and analysis of cell-cell communication using CellChat**. bioRxiv 2020.07.21.214387; doi: https://doi.org/10.1101/2020.07.21.214387


