
# CellChat: Inference and analysis of cell-cell communication

<p align="center">
  <img width="200"  src="https://github.com/sqjin/CellChat/blob/master/CellChat_Logo.png">
</p>

## Update
May 07, 2022 (Version 1.4.0)

We recently utilized CellChat to study the aging-induced signaling changes and thus updated some functions for better interpreting the inferred cell-cell communication. For the version history and detailed important changes, please see the [NEWS file](https://github.com/sqjin/CellChat/blob/master/NEWS.md).


## Capabilities
In addition to infer the intercellular communication from any given scRNA-seq data, CellChat provides functionality for further data exploration, analysis, and visualization. 

- It is able to analyze cell-cell communication for continuous states along cellular development trajectories.
- It can quantitatively characterize and compare the inferred cell-cell communication networks using an integrated approach by combining social network analysis, pattern recognition, and manifold learning approaches.
- It provides an easy-to-use tool for extracting and visualizing high-order information of the inferred networks. For example, it allows ready prediction of major signaling inputs and outputs for all cell populations and how these populations and signals coordinate together for functions.
- It provides several visualization outputs to facilitate intuitive user-guided data interpretation.

Check out [our paper (Jin et al., Nature Communications, 2021)](https://www.nature.com/articles/s41467-021-21246-9) for the detailed methods and applications.

## Installation

CellChat R package can be easily installed from Github using devtools:  

```
devtools::install_github("sqjin/CellChat")
```
**Please make sure you have installed the correct version of `NMF` and `circlize` package**. See instruction below. 

### Installation of other dependencies
- Install [NMF (>= 0.23.0)](http://renozao.github.io/NMF/devel/PAGE-INSTALLATION.html) using `install.packages('NMF')`. Please check [here](https://github.com/sqjin/CellChat/issues/16) for other solutions if you encounter any issue. You might can set `Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS=TRUE)` if it throws R version error. 
- Install [circlize (>= 0.4.12)](https://github.com/jokergoo/circlize) using `devtools::install_github("jokergoo/circlize")` if you encounter any issue.
- Install [ComplexHeatmap](https://github.com/jokergoo/ComplexHeatmap) using `devtools::install_github("jokergoo/ComplexHeatmap")` if you encounter any issue.
- Install UMAP python pacakge for dimension reduction: ```pip install umap-learn```. Please check [here](https://github.com/lmcinnes/umap) if you encounter any issue. 

Some users might have issues when installing CellChat pacakge due to different operating systems and new R version. Please check the following solutions:

- **Installation on Mac OX with R > 3.6**: Please re-install [Xquartz](https://community.rstudio.com/t/imager-package-does-not-work-in-r-3-6-1/38119).
- **Installation on Windows, Linux and Centos**: Please check the solution for [Windows](https://github.com/sqjin/CellChat/issues/5) and [Linux](https://github.com/sqjin/CellChat/issues/131).  



## Tutorials
Please check the tutorial directory of the repo.

- [Full tutorial for CellChat analysis of a single dataset with detailed explanation of each function](https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/CellChat-vignette.html)
- [Full tutorial for comparison analysis of multiple datasets](https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/Comparison_analysis_of_multiple_datasets.html)
- [Comparison analysis of multiple datasets with different cellular compositions](https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/Comparison_analysis_of_multiple_datasets_with_different_cellular_compositions.html)
- [Interface with other single-cell analysis toolkits (e.g., Seurat, SingleCellExperiment, Scanpy)](https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/Interface_with_other_single-cell_analysis_toolkits.html)
- [Tutorial for updating ligand-receptor database CellChatDB](https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/Update-CellChatDB.html)

<p align="center">
  <img width="700"  src="https://github.com/sqjin/CellChat/blob/master/overview_CellChat.png">
</p>


## Web-based “CellChat Explorer” 

We build a user-friendly web-based “[CellChat Explorer](http://www.cellchat.org/)” that contains two major components:
- **Ligand-Receptor Interaction Explorer** that allows easy exploration of our novel ligand-receptor interaction database, a comprehensive recapitulation of known molecular compositions including multimeric complexes and co-factors. *Our database CellChatDB is a manually curated database of literature-supported ligand-receptor interactions in both **human and mouse***.
- **Cell-Cell Communication Atlas Explorer** that allows easy exploration of the cell-cell communication for any given scRNA-seq dataset that has been processed by our R toolkit CellChat.  

We also developed a [standalone CellChat Shiny App](https://github.com/sqjin/CellChatShiny) for our Cell-Cell Communication Atlas Explorer. 


## Help, Suggestion and Contribution
If you have any question, comment or suggestion, please use github issue tracker to report coding related [issues](https://github.com/sqjin/CellChat/issues) of CellChat. 

### Before reporting an issue
- First **check the GitHub [issues](https://github.com/sqjin/CellChat/issues)** to see if the same or a similar issues has been reported and resolved. This relieves the developers from addressing the same issues and helps them focus on adding new features!
- The best way to figure out the issues is **running the sources codes** of the specific functions by yourself. This will also relieve the developers and helps them focus on the common issues! I am sorry, but I have to say I have no idea on many errors except that I can reproduce the issues. 
- Minimal and **reproducible example** are required when filing a GitHub issue. In certain cases, please share your CellChat object and related codes to reproduce the issues. 
- Users are encouraged to discuss issues and bugs using the github [issues](https://github.com/sqjin/CellChat/issues) instead of email exchanges.

### Contribution
CellChat is an open source software package and any contribution is highly appreciated! 

We use GitHub's [Pull Request](https://github.com/sqjin/CellChat/pulls) mechanism for reviewing and accepting submissions of any contribution. Issue a pull request on the GitHub website to request that we merge your branch's changes into CellChat's master branch. Be sure to include a description of your changes in the pull request, as well as any other information that will help the CellChat developers involved in reviewing your code. 

## System Requirements
- Hardware requirements: CellChat package requires only a standard computer with enough RAM to support the in-memory operations.

- Software requirements: This package is supported for macOS, Windows and Linux. The package has been tested on macOS: Mojave (10.14.5) and Windows 10. Dependencies of CellChat package are indicated in the Description file, and can be automatically installed when installing CellChat pacakge. CellChat can be installed on a normal computer within few mins.

## How to cite?
Suoqin Jin et al., Inference and analysis of cell-cell communication using CellChat. Nature Communications, 12:1088 (2021). https://www.nature.com/articles/s41467-021-21246-9 

![](https://api.visitorbadge.io/api/VisitorHit?user=sqjin&repo=CellChat&countColor=%237B1E7A)

