
---
title: "CellChat: Inference and analysis of cell-cell communication"
output:
  html_document:
    toc: true
    theme: united
mainfont: Arial
---

<p align="center">
  <img width="250"  src="https://github.com/sqjin/CellChat/blob/master/CellChat_Logo.png">
</p>

## Important update
Janunary 05, 2021 (Version 0.5.0)

* Enhanced documentation of functions and tutorials
* New features for comparison analysis of multiple datasets
* Support for creating a new CellChat object from Seurat V3 or SingleCellExperiment object
* Slight changes of CellChat object (Please update your previously calculated CellChat object via `updateCellChat()`)

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


## Installation

CellChat R package can be easily installed from Github using devtools:  

```
devtools::install_github("sqjin/CellChat")
```
CellChat can be installed on a normal computer within few mins. 

### Installation of other dependencies
- Install [NMF (>= 0.23.0)](http://renozao.github.io/NMF/devel/PAGE-INSTALLATION.html) using `install.packages('NMF')`. Please check [here](https://github.com/sqjin/CellChat/issues/16) for other solutions if you encounter any issue. You might can set `Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS=TRUE)` if it throws R version error. 
- Install [ComplexHeatmap](https://github.com/jokergoo/ComplexHeatmap) using `devtools::install_github("jokergoo/ComplexHeatmap")` if you encounter any issue.
- Install [circlize (>= 0.4.12)](https://github.com/jokergoo/circlize) using `devtools::install_github("jokergoo/circlize")` if you encounter any issue.
- Install UMAP python pacakge for dimension reduction: ```pip install umap-learn```. Please check [here](https://github.com/lmcinnes/umap) if you encounter any issue. 

Some users might have issues when installing CellChat pacakge due to different operating systems and new R version. Please check the following solutions:

- **Installation on Mac OX with R > 3.6**: Please re-install [Xquartz](https://community.rstudio.com/t/imager-package-does-not-work-in-r-3-6-1/38119).
- **Installation on Windows, Linux and Centos**: Please check the solution [here](https://github.com/sqjin/CellChat/issues/5).  



## Tutorials
Please check the tutorial directory of the repo.

- [Full tutorial for CellChat analysis of a single dataset with detailed explanation of each function](https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/CellChat-vignette.html)
- [Full tutorial for comparison analysis of multiple datasets](https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/Comparison_analysis_of_multiple_datasets.html)
- [Comparison analysis of multiple datasets with different cellular compositions](https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/Comparison_analysis_of_multiple_datasets_with_different_cellular_compositions.html)
- [Interface with other single-cell analysis toolkits (e.g., Seurat, SingleCellExperiment, Scanpy)](https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/Interface_with_other_single-cell_analysis_toolkitshtml)
- [Tutorial for updating CellChatDB](https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/Update-CellChatDB.html)



## System Requirements
- Hardware requirements: CellChat package requires only a standard computer with enough RAM to support the in-memory operations.

- Software requirements: This package is supported for macOS, Windows and Linux. The package has been tested on macOS: Mojave (10.14.5) and Windows 10. Dependencies of CellChat package are indicated in the Description file, and can be automatically installed when installing CellChat pacakge. 

## Help
If you have any questions, comments or suggestions, please contact cellchat.package@gmail.com.

## Preprint
Check out [our preprint on bioRxiv](https://www.biorxiv.org/content/10.1101/2020.07.21.214387v1) for the detailed methods and applications.

Suoqin Jin, Christian F. Guerrero-Juarez, Lihua Zhang, Ivan Chang, Peggy Myung, Maksim V. Plikus, Qing Nie. Inference and analysis of cell-cell communication using CellChat. bioRxiv 2020.07.21.214387; doi: https://doi.org/10.1101/2020.07.21.214387


