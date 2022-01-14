#' Normalize data using a scaling factor
#'
#' @param data.raw input raw data
#' @param scale.factor the scaling factor used for each cell
#' @param do.log whether do log transformation with pseudocount 1
#' @export
#'
normalizeData <- function(data.raw, scale.factor = 10000, do.log = TRUE) {
  # Scale counts within a sample
  library.size <- Matrix::colSums(data.raw)
  #scale.factor <- median(library.size)
  expr <- Matrix::t(Matrix::t(data.raw) / library.size) * scale.factor
  if (do.log) {
    data.norm <-log1p(expr)
  }
  return(data.norm)
}


#' Scale the data
#'
#' @param data.use input data
#' @param do.center whether center the values
#' @export
#'
scaleData <- function(data.use, do.center = T) {
  data.use <- Matrix::t(scale(Matrix::t(data.use), center = do.center, scale = TRUE))
  return(data.use)
}


#' Scale a data matrix
#'
#' @param x data matrix
#' @param scale the method to scale the data
#' @param na.rm whether remove na
#' @importFrom Matrix rowMeans colMeans rowSums colSums
#' @return
#' @export
#'
#' @examples
scaleMat <- function(x, scale, na.rm=TRUE){

  av <- c("none", "row", "column", 'r1', 'c1')
  i <- pmatch(scale, av)
  if(is.na(i) )
    stop("scale argument shoud take values: 'none', 'row' or 'column'")
  scale <- av[i]

  switch(scale, none = x
         , row = {
           x <- sweep(x, 1L, rowMeans(x, na.rm = na.rm), '-',check.margin = FALSE)
           sx <- apply(x, 1L, sd, na.rm = na.rm)
           sweep(x, 1L, sx, "/", check.margin = FALSE)
         }
         , column = {
           x <- sweep(x, 2L, colMeans(x, na.rm = na.rm), '-',check.margin = FALSE)
           sx <- apply(x, 2L, sd, na.rm = na.rm)
           sweep(x, 2L, sx, "/", check.margin = FALSE)
         }
         , r1 = sweep(x, 1L, rowSums(x, na.rm = na.rm), '/', check.margin = FALSE)
         , c1 = sweep(x, 2L, colSums(x, na.rm = na.rm), '/', check.margin = FALSE)
  )
}

#' Downsampling single cell data using geometric sketching algorithm
#'
#' USERs need to install the python package `pip install geosketch` (https://github.com/brianhie/geosketch)
#'
#' @param object A data matrix (should have row names; samples in rows, features in columns) or a Seurat object.
#'
#' When object is a PCA or UMAP space, please set `do.PCA = FALSE`
#'
#' When object is a data matrix (cells in rows and genes in columns), it is better to use the highly variable genes. PCA will be done on this input data matrix.
#' @param percent the percent of data to sketch
#' @param idents A vector of identity classes to keep for sketching
#' @param do.PCA whether doing PCA on the input data
#' @param dimPC the number of components to use
#' @importFrom reticulate import
#' @return A vector of cell names to use for downsampling
#' @export
#'
sketchData <- function(object, percent, idents = NULL, do.PCA = TRUE, dimPC = 30) {
  # pip install geosketch
  geosketch <- reticulate::import('geosketch')
  if (is(object,"Seurat")) {
    sketch.size <- as.integer(percent*ncol(object))
    if (!is.null(idents)) {
      object <- subset(object, idents = idents)
    }
    object <- object %>% #Seurat::NormalizeData(verbose = FALSE) %>%
      FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
      RunPCA(pc.genes = object@var.genes, npcs = dimPC, verbose = FALSE)

    X.pcs <- object@reductions$pca@cell.embeddings
    cells.all <- Cells(object)

  } else {
    # Get top PCs
    if (do.PCA) {
      X.pcs <- runPCA(object, dimPC = dimPC)
    } else {
      X.pcs <- object
    }

    # Sketch percent of data.
    sketch.size <- as.integer(percent*nrow(X))
    cells.all <- rownames(object)
  }
  sketch.index <- geosketch$gs(X.pcs, sketch.size)
  sketch.index <- unlist(sketch.index) + 1
  sketch.cells <- cells.all[sketch.index]
  return(sketch.cells)
}


#' Add the cell information into meta slot
#'
#' @param object CellChat object
#' @param meta cell information to be added
#' @param meta.name the name of column to be assigned
#'
#' @return
#' @export
#'
#' @examples
addMeta <- function(object, meta, meta.name = NULL) {
  if (is.null(x = meta.name) && is.atomic(x = meta)) {
    stop("'meta.name' must be provided for atomic meta types (eg. vectors)")
  }
  if (inherits(x = meta, what = c("matrix", "Matrix"))) {
    meta <- as.data.frame(x = meta)
  }

  if (is.null(x = meta.name)) {
    meta.name <- names(meta)
  } else {
    names(meta) <- meta.name
  }
  object@meta <- meta
  return(object)
}


#' Set the default identity of cells
#' @param object CellChat object
#' @param ident.use the name of the variable in object.meta;
#' @param levels set the levels of factor
#' @param display.warning whether display the warning message
#' @return
#' @export
#'
#' @examples
setIdent <- function(object, ident.use = NULL, levels = NULL, display.warning = TRUE){
  if (!is.null(ident.use)) {
    object@idents <- as.factor(object@meta[[ident.use]])
  }

  if (!is.null(levels)) {
    object@idents <- factor(object@idents, levels = levels)
  }
  if ("0" %in% as.character(object@idents)) {
    stop("Cell labels cannot contain `0`! ")
  }
  if (length(object@net) > 0) {
    if (all(dimnames(object@net$prob)[[1]] %in% levels(object@idents) )) {
      message("Reorder cell groups! ")
      cat("The cell group order before reordering is ", dimnames(object@net$prob)[[1]],'\n')
      # idx <- match(dimnames(object@net$prob)[[1]], levels(object@idents))
      idx <- match(levels(object@idents), dimnames(object@net$prob)[[1]])
      object@net$prob <- object@net$prob[idx, , ]
      object@net$prob <- object@net$prob[, idx, ]
      object@net$pval <- object@net$pval[idx, , ]
      object@net$pval <- object@net$pval[, idx, ]
      cat("The cell group order after reordering is ", dimnames(object@net$prob)[[1]],'\n')
    } else {
      message("Rename cell groups but do not change the order! ")
      cat("The cell group order before renaming is ", dimnames(object@net$prob)[[1]],'\n')
      dimnames(object@net$prob) <- list(levels(object@idents), levels(object@idents), dimnames(object@net$prob)[[3]])
      dimnames(object@net$pval) <- dimnames(object@net$prob)
      cat("The cell group order after renaming is ", dimnames(object@net$prob)[[1]],'\n')
    }
    if (display.warning) {
      warning("All the calculations after `computeCommunProb` should be re-run!!
    These include but not limited to `computeCommunProbPathway`,`aggregateNet`, and `netAnalysis_computeCentrality`.")
    }


  }
  return(object)
}


#' Update and re-order the cell group names after running `computeCommunProb`
#'
#' @param object CellChat object
#' @param old.cluster.name A vector defining old cell group labels in `object@idents`; Default = NULL, which will use `levels(object@idents)`
#' @param new.cluster.name A vector defining new cell group labels to rename
#' @param new.order reset order of cell group labels
#' @param new.cluster.metaname assign a name of the new labels, which will be the column name of new labels in `object@meta`
#' @return An updated CellChat object
#' @export
#'
updateClusterLabels <- function(object, old.cluster.name = NULL, new.cluster.name = NULL, new.order = NULL, new.cluster.metaname = "new.labels") {
  if (is.null(old.cluster.name)) {
    old.cluster.name <- levels(object@idents)
  }
  if (new.cluster.metaname %in% colnames(object@meta)) {
    stop("Please define another `new.cluster.metaname` as it exists in `colnames(object@meta)`!")
  }
  if (!is.null(new.cluster.name)) {
    labels.new <- plyr::mapvalues(object@idents, from = old.cluster.name, to = new.cluster.name)
    object@meta[[new.cluster.metaname]] <- labels.new
    object <- setIdent(object, ident.use = new.cluster.metaname, display.warning = FALSE)
  } else {
    new.cluster.metaname <- NULL
    cat("Only reorder cell groups but do not rename cell groups!")
  }

  if (!is.null(new.order)) {
    object <- setIdent(object, ident.use = new.cluster.metaname, levels = new.order, display.warning = FALSE)
  }
  message("We now re-run computeCommunProbPathway`,`aggregateNet`, and `netAnalysis_computeCentrality`...")
  object <- computeCommunProbPathway(object)
  ## calculate the aggregated network by counting the number of links or summarizing the communication probability
  object <- aggregateNet(object)
  # network importance analysis
  object <-netAnalysis_computeCentrality(object, slot.name = "netP")
  return(object)
}





#' Subset the expression data of signaling genes for saving computation cost
#'
#' @param object CellChat object
#' @param features default = NULL: subset the expression data of signaling genes in CellChatDB.use
#'
#' @return An updated CellChat object by assigning a subset of the data into the slot `data.signaling`
#' @export
#'
subsetData <- function(object, features = NULL) {
  if (is.null(features)) {
    DB <- object@DB
    gene.use_input <- extractGene(DB)
    gene.use <- intersect(gene.use_input, rownames(object@data))
  } else {
    gene.use <- intersect(features, rownames(object@data))
  }
  object@data.signaling <- object@data[rownames(object@data) %in% gene.use, ]
  return(object)
}



#' Identify over-expressed signaling genes associated with each cell group
#'
#' USERS can use customized gene set as over-expressed signaling genes by setting `object@var.features[[features.name]] <- features.sig`
#' The Bonferroni corrected/adjusted p value can be obtained via `object@var.features[[paste0(features.name, ".info")]]`. Note that by default `features.name = "features"`
#'
#' @param object CellChat object
#' @param data.use a customed data matrix. Default: data.use = NULL and the expression matrix in the slot 'data.signaling' is used
#' @param group.by cell group information; default is `object@idents`; otherwise it should be one of the column names of the meta slot
#' @param idents.use a subset of cell groups used for analysis
#' @param invert whether invert the idents.use
#' @param group.dataset dataset origin information in a merged CellChat object; set it as one of the column names of meta slot when identifying the highly enriched genes in one dataset for each cell group
#' @param pos.dataset the dataset name used for identifying highly enriched genes in this dataset for each cell group
#' @param features.name a char name used for storing the over-expressed signaling genes in `object@var.features[[features.name]]`
#' @param only.pos Only return positive markers
#' @param features features used for identifying Over Expressed genes. default use all features
#' @param return.object whether return the object; otherwise return a data frame consisting of over-expressed signaling genes associated with each cell group
#' @param thresh.pc Threshold of the percent of cells expressed in one cluster
#' @param thresh.fc Threshold of Log Fold Change
#' @param thresh.p Threshold of p-values
#' @importFrom future nbrOfWorkers
#' @importFrom pbapply pbsapply
#' @importFrom future.apply future_sapply
#' @importFrom stats sd wilcox.test
#' @importFrom stats p.adjust
#'
#' @return A CellChat object or a data frame. If returning a CellChat object, two new elements named 'features.name' and paste0(features.name, ".info") will be added into the list `object@var.features`
#' `object@var.features[[features.name]]` is a vector consisting of the identified over-expressed signaling genes;
#' `object@var.features[[paste0(features.name, ".info")]]` is a data frame returned from the differential expression analysis
#' @export
#'
identifyOverExpressedGenes <- function(object, data.use = NULL, group.by = NULL, idents.use = NULL, invert = FALSE, group.dataset = NULL, pos.dataset = NULL, features.name = "features",  only.pos = TRUE, features = NULL, return.object = TRUE,
                                       thresh.pc = 0, thresh.fc = 0, thresh.p = 0.05) {
  if (!is.list(object@var.features)) {
    stop("Please update your CellChat object via `updateCellChat()`")
  }
  if (is.null(data.use)) {
    X <- object@data.signaling
    if (nrow(X) < 3) {stop("Please check `object@data.signaling` and ensure that you have run `subsetData` and that the data matrix `object@data.signaling` looks OK.")}
  } else {
    X <- data.use
  }

  if (is.null(features)) {
    features.use <- row.names(X)
  } else {
    features.use <- intersect(features, row.names(X))
  }
  data.use <- X[features.use,]
  data.use <- as.matrix(data.use)

  if (is.null(group.by)) {
    labels <- object@idents
    if (!is.factor(labels)) {
      message("Use the joint cell labels from the merged CellChat object")
      labels <- object@idents$joint
    }
  } else {
    labels <- object@meta[[group.by]]
  }
  if (!is.factor(labels)) {
    labels <- factor(labels)
  }
  level.use <- levels(labels)[levels(labels) %in% unique(labels)]
  if (!is.null(idents.use)) {
    if (invert) {
      level.use <- level.use[!(level.use %in% idents.use)]
    } else {
      level.use <- level.use[level.use %in% idents.use]
    }
  }
  numCluster <- length(level.use)

  if (!is.null(group.dataset)) {
    labels.dataset <- as.character(object@meta[[group.dataset]])
    if (!(pos.dataset %in% unique(labels.dataset))) {
      cat("Please set pos.dataset to be one of the following dataset names: ", unique(as.character(labels.dataset)))
      stop()
    }
  }

  my.sapply <- ifelse(
    test = future::nbrOfWorkers() == 1,
    yes = pbapply::pbsapply,
    no = future.apply::future_sapply
  )

  mean.fxn <- function(x) {
    return(log(x = mean(x = expm1(x = x)) + 1))
  }
  labels <- as.character(labels)
  genes.de <- vector("list", length = numCluster)
  for (i in 1:numCluster) {
    features <- features.use
    if (is.null(group.dataset)) {
      cell.use1 <- which(labels == level.use[i])
      cell.use2 <- base::setdiff(1:length(labels), cell.use1)
    } else {
      cell.use1 <- which((labels == level.use[i]) & (labels.dataset == pos.dataset))
      cell.use2 <- which((labels == level.use[i]) & (labels.dataset != pos.dataset))
    }

    # feature selection (based on percentages)
    thresh.min <- 0
    pct.1 <- round(
      x = rowSums(data.use[features, cell.use1, drop = FALSE] > thresh.min) /
        length(x = cell.use1),
      digits = 3
    )
    pct.2 <- round(
      x = rowSums(data.use[features, cell.use2, drop = FALSE] > thresh.min) /
        length(x = cell.use2),
      digits = 3
    )
    data.alpha <- cbind(pct.1, pct.2)
    colnames(x = data.alpha) <- c("pct.1", "pct.2")
    alpha.min <- apply(X = data.alpha, MARGIN = 1, FUN = max)
    names(x = alpha.min) <- rownames(x = data.alpha)
    features <- names(x = which(x = alpha.min > thresh.pc))
    if (length(x = features) == 0) {
      #stop("No features pass thresh.pc threshold")
      next
    }

    # feature selection (based on average difference)
    data.1 <- apply(X = data.use[features, cell.use1, drop = FALSE],MARGIN = 1,FUN = mean.fxn)
    data.2 <- apply(X = data.use[features, cell.use2, drop = FALSE],MARGIN = 1,FUN = mean.fxn)
    FC <- (data.1 - data.2)
    if (only.pos) {
      features.diff <- names(which(FC > thresh.fc))
    } else {
      features.diff <- names(which(abs(FC) > thresh.fc))
    }

    features <- intersect(x = features, y = features.diff)
    if (length(x = features) == 0) {
      #  stop("No features pass thresh.fc threshold")
      next
    }

    data1 <- data.use[features, cell.use1, drop = FALSE]
    data2 <- data.use[features, cell.use2, drop = FALSE]

    pvalues <- unlist(
      x = my.sapply(
        X = 1:nrow(x = data1),
        FUN = function(x) {
          # return(wilcox.test(data1[x, ], data2[x, ], alternative = "greater")$p.value)
          return(wilcox.test(data1[x, ], data2[x, ])$p.value)
        }
      )
    )

    pval.adj = stats::p.adjust(
      p = pvalues,
      method = "bonferroni",
      n = nrow(X)
    )
    genes.de[[i]] <- data.frame(clusters = level.use[i], features = as.character(rownames(data1)), pvalues = pvalues, logFC = FC[features], data.alpha[features,, drop = F],pvalues.adj = pval.adj, stringsAsFactors = FALSE)
  }

  markers.all <- data.frame()
  for (i in 1:numCluster) {
    gde <- genes.de[[i]]
    if (!is.null(gde)) {
      gde <- gde[order(gde$pvalues, -gde$logFC), ]
      gde <- subset(gde, subset = pvalues < thresh.p)
      if (nrow(gde) > 0) {
        markers.all <- rbind(markers.all, gde)
      }
    }
  }
  if (only.pos & nrow(markers.all) > 0) {
    markers.all <- subset(markers.all, subset = logFC > 0)
  }
  if (!is.null(group.dataset)) {
    markers.all$datasets[markers.all$logFC > 0] <- pos.dataset
    markers.all$datasets[markers.all$logFC < 0] <- setdiff(unique(labels.dataset), pos.dataset)
    markers.all$datasets <- factor(markers.all$datasets, levels = levels(factor(object@meta[[group.dataset]])))
    markers.all <- markers.all[order(markers.all$datasets, markers.all$pvalues, -markers.all$logFC), ]
  }
  markers.all$features <- as.character(markers.all$features)

  features.sig <- markers.all$features
  object@var.features[[features.name]] <- features.sig
  features.name <- paste0(features.name, ".info")
  object@var.features[[features.name]] <- markers.all

  if (return.object) {
    return(object)
  } else {
    return(markers.all)
  }
}


#' Identify over-expressed ligands and (complex) receptors associated with each cell group
#'
#' This function identifies the over-expressed ligands and (complex) receptors based on the identified signaling genes from 'identifyOverExpressedGenes'.
#'
#' @param object CellChat object
#' @param features.name a char name used for storing the over-expressed ligands and receptors in `object@var.features[[paste0(features.name, ".LR")]]`
#' @param features a vector of features to use. default use all over-expressed genes in `object@var.features[[features.name]]`
#' @param return.object whether returning a CellChat object. If FALSE, it will return a data frame containing over-expressed ligands and (complex) receptors associated with each cell group
#' @importFrom future nbrOfWorkers
#' @importFrom future.apply future_sapply
#' @importFrom pbapply pbsapply
#' @importFrom dplyr select
#'
#' @return A CellChat object or a data frame. If returning a CellChat object, a new element named paste0(features.name, ".LR") will be added into the list `object@var.features`
#' @export
#'
identifyOverExpressedLigandReceptor <- function(object, features.name = "features", features = NULL, return.object = TRUE) {

  features.name.LR <- paste0(features.name, ".LR")
  features.name <- paste0(features.name, ".info")
  DB <- object@DB
  interaction_input <- DB$interaction
  complex_input <- DB$complex
  pairLR <- select(interaction_input, ligand, receptor)
  LR.use <- unique(c(pairLR$ligand, pairLR$receptor))
  if (is.null(features)) {
    if (is.list(object@var.features)) {
      markers.all <- object@var.features[[features.name]] # use the updated CellChat object 12/2020
    } else {
      stop("Please update your CellChat object via `updateCellChat()`")
    }

  } else {
    features.use <- features
    rm(features)
    markers.all <- subset(markers.all, subset = features %in% features.use)
  }

  my.sapply <- ifelse(
    test = future::nbrOfWorkers() == 1,
    yes = pbapply::pbsapply,
    no = future.apply::future_sapply
  )
  complexSubunits <- complex_input[, grepl("subunit" , colnames(complex_input))]

  markers.all.new <- data.frame()
  for (i in 1:nrow(markers.all)) {
    if (markers.all$features[i] %in% LR.use) {
      markers.all.new <- rbind(markers.all.new, markers.all[i, , drop = FALSE])
    } else {
      index.sig <- unlist(
        x = my.sapply(
          X = 1:nrow(complexSubunits),
          FUN = function(x) {
            complexsubunitsV <- unlist(complexSubunits[x,], use.names = F)
            complexsubunitsV <- complexsubunitsV[complexsubunitsV != ""]
            if (markers.all$features[i] %in% complexsubunitsV) {
              return(x)
            }
          }
        )
      )
      complexSubunits.sig <- rownames(complexSubunits[index.sig,])
      markers.all.complex <- data.frame()
      for (j in 1:length(complexSubunits.sig)) {
        markers.all.complex <- rbind(markers.all.complex, markers.all[i, , drop = FALSE])
      }
      markers.all.complex$features <- complexSubunits.sig
      markers.all.new <- rbind(markers.all.new, markers.all.complex)
    }
  }

  object@var.features[[features.name.LR]] <- markers.all.new

  if (return.object) {
    return(object)
  } else {
    return(markers.all.new)
  }
}



#' Identify over-expressed ligand-receptor interactions (pairs) within the used CellChatDB
#'
#' @param object CellChat object
#' @param features.name a char name used for assess the results in `object@var.features[[features.name]]`
#' @param features a vector of features to use. default use all over-expressed genes in `object@var.features[[features.name]]`
#' @param return.object whether returning a CellChat object. If FALSE, it will return a data frame containing the over-expressed ligand-receptor pairs
#' @importFrom future nbrOfWorkers
#' @importFrom future.apply future_sapply
#' @importFrom pbapply pbsapply
#' @importFrom dplyr select
#'
#' @return A CellChat object or a data frame. If returning a CellChat object, a new element named 'LRsig' will be added into the list `object@LR`
#' @export
#'
identifyOverExpressedInteractions <- function(object, features.name = "features", features = NULL, return.object = TRUE) {
  gene.use <- row.names(object@data.signaling)
  DB <- object@DB
  if (is.null(features)) {
    if (is.list(object@var.features)) {
      features.sig <- object@var.features[[features.name]] # use the updated CellChat object 12/2020
    } else {
      stop("Please update your CellChat object via `updateCellChat()`")
    }

  } else {
    features.sig <- features
  }

  interaction_input <- DB$interaction
  complex_input <- DB$complex
  my.sapply <- ifelse(
    test = future::nbrOfWorkers() == 1,
    yes = pbapply::pbsapply,
    no = future.apply::future_sapply
  )
  complexSubunits <- complex_input[, grepl("subunit" , colnames(complex_input))]
  index.sig <- unlist(
    x = my.sapply(
      X = 1:nrow(complexSubunits),
      FUN = function(x) {
        complexsubunitsV <- unlist(complexSubunits[x,], use.names = F)
        complexsubunitsV <- complexsubunitsV[complexsubunitsV != ""]
        if (length(intersect(complexsubunitsV, features.sig)) > 0 & all(complexsubunitsV %in% gene.use)) {
          return(x)
        }
      }
    )
  )
  complexSubunits.sig <- complexSubunits[index.sig,]

  pairLR <- select(interaction_input, ligand, receptor)

  index.sig <- unlist(
    x = my.sapply(
      X = 1:nrow(pairLR),
      FUN = function(x) {
        if (all(unlist(pairLR[x,], use.names = F) %in% c(features.sig, rownames(complexSubunits.sig)))) {
          return(x)
        }
      }
    )
  )
  pairLRsig <- interaction_input[index.sig, ]
  object@LR$LRsig <- pairLRsig
  if (return.object) {
    return(object)
  } else {
    return(pairLRsig)
  }
}


#' Project gene expression data onto a protein-protein interaction network
#'
#' A diffusion process is used to smooth genes’ expression values based on their neighbors’ defined in a high-confidence experimentally validated protein-protein network.
#'
#' This function is useful when analyzing single-cell data with shallow sequencing depth because the projection reduces the dropout effects of signaling genes, in particular for possible zero expression of subunits of ligands/receptors
#'
#' @param object  CellChat object
#' @param adjMatrix adjacency matrix of protein-protein interaction network to use
#' @param alpha numeric in [0,1] alpha = 0: no smoothing; a larger value alpha results in increasing levels of smoothing.
#' @param normalizeAdjMatrix    how to normalize the adjacency matrix
#'                              possible values are 'rows' (in-degree)
#'                              and 'columns' (out-degree)
#' @return a projected gene expression matrix
#' @export
#'
# This function is adapted from https://github.com/BIMSBbioinfo/netSmooth
projectData <- function(object, adjMatrix, alpha=0.5, normalizeAdjMatrix=c('rows','columns')){
  data <- as.matrix(object@data.signaling)
  normalizeAdjMatrix <- match.arg(normalizeAdjMatrix)
  stopifnot(is(adjMatrix, 'matrix') || is(adjMatrix, 'sparseMatrix'))
  stopifnot((is.numeric(alpha) && (alpha > 0 && alpha < 1)))
  if(sum(Matrix::rowSums(adjMatrix)==0)>0) stop("PPI cannot have zero rows/columns")
  if(sum(Matrix::colSums(adjMatrix)==0)>0) stop("PPI cannot have zero rows/columns")
  if(is.numeric(alpha)) {
    if(alpha<0 | alpha > 1) {
      stop('alpha must be between 0 and 1')
    }
    data.projected <- projectAndRecombine(data, adjMatrix, alpha,normalizeAdjMatrix=normalizeAdjMatrix)
  } else stop("unsupported alpha value: ", class(alpha))
  object@data.project <- data.projected
  return(object)
}

#' Perform network projecting on network when the network genes and the
#' experiment genes aren't exactly the same.
#'
#' The gene network might be defined only on a subset of genes that are
#' measured in any experiment. Further, an experiment might not measure all
#' genes that are present in the network. This function projects the experiment
#' data onto the gene space defined by the network prior to projecting. Then,
#' it projects the projected data back into the original dimansions.
#'
#' @param gene_expression  gene expession data to be projected
#'                         [N_genes x M_samples]
#' @param adj_matrix  adjacenty matrix of network to perform projecting over.
#'                    Will be column-normalized.
#'                    Rownames and colnames should be genes.
#' @param alpha  network projecting parameter (1 - restart probability in random
#'                walk model.
#' @param projecting.function  must be a function that takes in data, adjacency
#'                            matrix, and alpha. Will be used to perform the
#'                            actual projecting.
#' @param normalizeAdjMatrix    which dimension (rows or columns) should the
#'                              adjacency matrix be normalized by. rows
#'                              corresponds to in-degree, columns to
#'                              out-degree.
#' @return  matrix with network-projected gene expression data. Genes that are
#'          not present in projecting network will retain original values.
#' @keywords internal
#'
projectAndRecombine <- function(gene_expression, adj_matrix, alpha,
                                projecting.function=randomWalkBySolve,
                                normalizeAdjMatrix=c('rows','columns')) {
  normalizeAdjMatrix <- match.arg(normalizeAdjMatrix)
  gene_expression_in_A_space <- projectOnNetwork(gene_expression,rownames(adj_matrix))
  gene_expression_in_A_space_project <- projecting.function(gene_expression_in_A_space, adj_matrix, alpha, normalizeAdjMatrix)
  gene_expression_project <- projectFromNetworkRecombine(gene_expression, gene_expression_in_A_space_project)
  return(gene_expression_project)
}


#' Project the gene expression matrix onto a lower space
#' of the genes defined in the projecting network
#' @param gene_expression    gene expression matrix
#' @param new_features       the genes in the network, on which to project
#'                           the gene expression matrix
#' @param missing.value      value to assign to genes that are in network,
#'                           but missing from gene expression matrix
#' @return the gene expression matrix projected onto the gene space defined by new_features
#' @keywords internal
projectOnNetwork <- function(gene_expression, new_features, missing.value=0) {
  data_in_new_space = matrix(0, ncol=dim(gene_expression)[2], nrow=length(new_features))
  rownames(data_in_new_space) <- new_features
  colnames(data_in_new_space) <- colnames(gene_expression)
  genes_in_both <- intersect(rownames(data_in_new_space),rownames(gene_expression))
  data_in_new_space[genes_in_both,] <- gene_expression[genes_in_both,]
  genes_only_in_network <- setdiff(new_features, rownames(gene_expression))
  data_in_new_space[genes_only_in_network,] <- missing.value
  return(data_in_new_space)
}

#' project data on graph by solving the linear equation (I - alpha*A) * E_sm = E * (1-alpha)

#' @param E      initial data matrix [NxM]
#' @param A      adjacency matrix of graph to network project on will be column-normalized.
#' @param alpha  projecting coefficient (1 - restart probability of random walk)
#' @return network-projected gene expression
#' @keywords internal
randomWalkBySolve <- function(E, A, alpha, normalizeAjdMatrix=c('rows','columns')) {
  normalizeAjdMatrix <- match.arg(normalizeAjdMatrix)
  if (normalizeAjdMatrix=='rows') {
    Anorm <- l1NormalizeRows(A)
  } else if (normalizeAjdMatrix=='columns') {
    Anorm <- l1NormalizeColumns(A)
  }
  eye <- diag(dim(A)[1])
  AA <- eye - alpha*Anorm
  BB <- (1-alpha) * E
  return(solve(AA, BB))
}

#' Column-normalize a sparse, symmetric matrix (using the l1 norm) so that each
#' column sums to 1.
#'
#' @param A matrix
#' @usage l1NormalizeColumns(A)
#' @return column-normalized sparse matrix object
#' @keywords internal
l1NormalizeColumns <- function(A) {
  return(Matrix::t(Matrix::t(A)/Matrix::colSums(A)))
}

#' Row-normalize a sparse, symmetric matrix (using the l1 norm) so that each
#' row sums to 1.
#'
#' @param A matrix
#' @usage l1NormalizeRows(A)
#' @return row-normalized sparse matrix object
#' @keywords internal
l1NormalizeRows <- function(A) {
  return(A/Matrix::rowSums(A))
}

#' Combine gene expression from projected space (that of the network) with the
#' expression of genes that were not projected (not present in network)
#' @keywords internal
#' @param original_expression    the non-projected expression
#' @param projected_expression    the projected gene expression, in the space
#'                               of the genes defined by the network
#' @return a matrix in the dimensions of original_expression, where values that
#'         are present in projected_expression are copied from there.
projectFromNetworkRecombine <- function(original_expression, projected_expression) {
  data_in_original_space <- original_expression
  genes_in_both <- intersect(rownames(original_expression),rownames(projected_expression))
  data_in_original_space[genes_in_both,] <- as.matrix(projected_expression[genes_in_both,])
  return(data_in_original_space)
}


#' Dimension reduction using PCA
#'
#' @param data.use input data (samples in rows, features in columns)
#' @param do.fast whether do fast PCA
#' @param dimPC the number of components to keep
#' @param seed.use set a seed
#' @param weight.by.var whether use weighted pc.scores
#' @importFrom stats prcomp
#' @importFrom irlba irlba
#' @return
#' @export
#'
#' @examples
runPCA <- function(data.use, do.fast = T, dimPC = 50, seed.use = 42, weight.by.var = T) {
  set.seed(seed = seed.use)
  if (do.fast) {
    dimPC <- min(dimPC, ncol(data.use) - 1)
    pca.res <- irlba::irlba(data.use, nv = dimPC)
    sdev <- pca.res$d/sqrt(max(1, nrow(data.use) - 1))
    if (weight.by.var){
      pc.scores <- pca.res$u %*% diag(pca.res$d)
    } else {
      pc.scores <- pca.res$u
    }
  } else {
    dimPC <- min(dimPC, ncol(data.use) - 1)
    pca.res <- stats::prcomp(x = data.use, rank. = dimPC)
    sdev <- pca.res$sdev
    if (weight.by.var) {
      pc.scores <- pca.res$x %*% diag(pca.res$sdev[1:dimPC]^2)
    } else {
      pc.scores <- pca.res$x
    }
  }
  rownames(pc.scores) <- rownames(data.use)
  colnames(pc.scores) <- paste0('PC', 1:ncol(pc.scores))
  return(pc.scores)
}


#' Run UMAP
#' @param data.use input data matrix
#' @param n_neighbors This determines the number of neighboring points used in
#' local approximations of manifold structure. Larger values will result in more
#' global structure being preserved at the loss of detailed local structure. In general this parameter should often be in the range 5 to 50.
#' @param n_components The dimension of the space to embed into.
#' @param metric This determines the choice of metric used to measure distance in the input space.
#' @param n_epochs the number of training epochs to be used in optimizing the low dimensional embedding. Larger values result in more accurate embeddings. If NULL is specified, a value will be selected based on the size of the input dataset (200 for large datasets, 500 for small).
#' @param learning_rate The initial learning rate for the embedding optimization.
#' @param min_dist This controls how tightly the embedding is allowed compress points together.
#' Larger values ensure embedded points are moreevenly distributed, while smaller values allow the
#' algorithm to optimise more accurately with regard to local structure. Sensible values are in the range 0.001 to 0.5.
#' @param spread he effective scale of embedded points. In combination with min.dist this determines how clustered/clumped the embedded points are.
#' @param set_op_mix_ratio Interpolate between (fuzzy) union and intersection as the set operation used to combine local fuzzy simplicial sets to obtain a global fuzzy simplicial sets.
#' @param local_connectivity The local connectivity required - i.e. the number of nearest neighbors
#' that should be assumed to be connected at a local level. The higher this value the more connected
#' the manifold becomes locally. In practice this should be not more than the local intrinsic dimension of the manifold.
#' @param repulsion_strength Weighting applied to negative samples in low dimensional embedding
#' optimization. Values higher than one will result in greater weight being given to negative samples.
#' @param negative_sample_rate The number of negative samples to select per positive sample in the
#' optimization process. Increasing this value will result in greater repulsive force being applied, greater optimization cost, but slightly more accuracy.
#' @param a More specific parameters controlling the embedding. If NULL, these values are set automatically as determined by min. dist and spread.
#' @param b More specific parameters controlling the embedding. If NULL, these values are set automatically as determined by min. dist and spread.
#' @param seed.use Set a random seed. By default, sets the seed to 42.
#' @param metric_kwds,angular_rp_forest,verbose other parameters used in UMAP
#' @import reticulate
#' @export
#'
runUMAP <- function(
  data.use,
  n_neighbors = 30L,
  n_components = 2L,
  metric = "correlation",
  n_epochs = NULL,
  learning_rate = 1.0,
  min_dist = 0.3,
  spread = 1.0,
  set_op_mix_ratio = 1.0,
  local_connectivity = 1L,
  repulsion_strength = 1,
  negative_sample_rate = 5,
  a = NULL,
  b = NULL,
  seed.use = 42L,
  metric_kwds = NULL,
  angular_rp_forest = FALSE,
  verbose = FALSE){
  if (!reticulate::py_module_available(module = 'umap')) {
    stop("Cannot find UMAP, please install through pip (e.g. pip install umap-learn or reticulate::py_install(packages = 'umap-learn')).")
  }
  set.seed(seed.use)
  reticulate::py_set_seed(seed.use)
  umap_import <- reticulate::import(module = "umap", delay_load = TRUE)
  umap <- umap_import$UMAP(
    n_neighbors = as.integer(n_neighbors),
    n_components = as.integer(n_components),
    metric = metric,
    n_epochs = n_epochs,
    learning_rate = learning_rate,
    min_dist = min_dist,
    spread = spread,
    set_op_mix_ratio = set_op_mix_ratio,
    local_connectivity = local_connectivity,
    repulsion_strength = repulsion_strength,
    negative_sample_rate = negative_sample_rate,
    a = a,
    b = b,
    metric_kwds = metric_kwds,
    angular_rp_forest = angular_rp_forest,
    verbose = verbose
  )
  Rumap <- umap$fit_transform
  umap_output <- Rumap(t(data.use))
  colnames(umap_output) <- paste0('UMAP', 1:ncol(umap_output))
  rownames(umap_output) <- colnames(data.use)
  return(umap_output)
}

.error_if_no_Seurat <- function() {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat installation required for working with Seurat objects")
  }
}


#' Color interpolation
#'
#' This function is modified from https://rdrr.io/cran/circlize/src/R/utils.R
#' Colors are linearly interpolated according to break values and corresponding colors through CIE Lab color space (`colorspace::LAB`) by default.
#' Values exceeding breaks will be assigned with corresponding maximum or minimum colors.
#'
#' @param breaks A vector indicating numeric breaks
#' @param colors A vector of colors which correspond to values in ``breaks``
#' @param transparency A single value in ``[0, 1]``. 0 refers to no transparency and 1 refers to full transparency
#' @param space color space in which colors are interpolated. Value should be one of "RGB", "HSV", "HLS", "LAB", "XYZ", "sRGB", "LUV", see `colorspace::color-class` for detail.
#' @importFrom colorspace coords RGB HSV HLS LAB XYZ sRGB LUV hex
#' @importFrom grDevices col2rgb
#' @return It returns a function which accepts a vector of numeric values and returns interpolated colors.
#' @export
#' @examples
#' \dontrun{
#' col_fun = colorRamp3(c(-1, 0, 1), c("green", "white", "red"))
#' col_fun(c(-2, -1, -0.5, 0, 0.5, 1, 2))
#' }
colorRamp3 = function(breaks, colors, transparency = 0, space = "LAB") {

  if(length(breaks) != length(colors)) {
    stop("Length of `breaks` should be equal to `colors`.\n")
  }

  colors = colors[order(breaks)]
  breaks = sort(breaks)

  l = duplicated(breaks)
  breaks = breaks[!l]
  colors = colors[!l]

  if(length(breaks) == 1) {
    stop("You should have at least two distinct break values.")
  }


  if(! space %in% c("RGB", "HSV", "HLS", "LAB", "XYZ", "sRGB", "LUV")) {
    stop("`space` should be in 'RGB', 'HSV', 'HLS', 'LAB', 'XYZ', 'sRGB', 'LUV'")
  }

  colors = t(grDevices::col2rgb(colors)/255)

  attr = list(breaks = breaks, colors = colors, transparency = transparency, space = space)

  if(space == "LUV") {
    i = which(apply(colors, 1, function(x) all(x == 0)))
    colors[i, ] = 1e-5
  }

  transparency = 1-ifelse(transparency > 1, 1, ifelse(transparency < 0, 0, transparency))[1]
  transparency_str = sprintf("%X", round(transparency*255))
  if(nchar(transparency_str) == 1) transparency_str = paste0("0", transparency_str)

  fun = function(x = NULL, return_rgb = FALSE, max_value = 1) {
    if(is.null(x)) {
      stop("Please specify `x`\n")
    }

    att = attributes(x)
    if(is.data.frame(x)) x = as.matrix(x)

    l_na = is.na(x)
    if(all(l_na)) {
      return(rep(NA, length(l_na)))
    }

    x2 = x[!l_na]

    x2 = ifelse(x2 < breaks[1], breaks[1],
                ifelse(x2 > breaks[length(breaks)], breaks[length(breaks)],
                       x2
                ))
    ibin = .bincode(x2, breaks, right = TRUE, include.lowest = TRUE)
    res_col = character(length(x2))
    for(i in unique(ibin)) {
      l = ibin == i
      res_col[l] = .get_color(x2[l], breaks[i], breaks[i+1], colors[i, ], colors[i+1, ], space = space)
    }
    res_col = paste(res_col, transparency_str[1], sep = "")

    if(return_rgb) {
      res_col = t(grDevices::col2rgb(as.vector(res_col), alpha = TRUE)/255)
      return(res_col)
    } else {
      res_col2 = character(length(x))
      res_col2[l_na] = NA
      res_col2[!l_na] = res_col

      attributes(res_col2) = att
      return(res_col2)
    }
  }

  attributes(fun) = attr
  return(fun)
}

.restrict_in = function(x, lower, upper) {
  x[x > upper] = upper
  x[x < lower] = lower
  x
}

# x: vector
# break1 single value
# break2 single value
# rgb1 vector with 3 elements
# rgb2 vector with 3 elements
.get_color = function(x, break1, break2, col1, col2, space) {

  col1 = colorspace::coords(as(colorspace::sRGB(col1[1], col1[2], col1[3]), space))
  col2 = colorspace::coords(as(colorspace::sRGB(col2[1], col2[2], col2[3]), space))

  res_col = matrix(ncol = 3, nrow = length(x))
  for(j in 1:3) {
    xx = (x - break2)*(col2[j] - col1[j]) / (break2 - break1) + col2[j]
    res_col[, j] = xx
  }

  res_col = get(space)(res_col)
  res_col = colorspace::coords(as(res_col, "sRGB"))
  res_col[, 1] = .restrict_in(res_col[,1], 0, 1)
  res_col[, 2] = .restrict_in(res_col[,2], 0, 1)
  res_col[, 3] = .restrict_in(res_col[,3], 0, 1)
  colorspace::hex(colorspace::sRGB(res_col))
}
