
#' The CellChat Class
#'
#' The CellChat object is created from a single-cell transcriptomic data matrix, Seurat V3 or SingleCellExperiment object.
#' When inputting an data matrix, it takes a digital data matrices as input. Genes should be in rows and cells in columns. rownames and colnames should be included.
#' The class provides functions for data preprocessing, intercellular communication network inference, communication network analysis, and visualization.
#'
#'
#'# Class definitions
#' @importFrom methods setClassUnion
#' @importClassesFrom Matrix dgCMatrix
setClassUnion(name = 'AnyMatrix', members = c("matrix", "dgCMatrix"))
setClassUnion(name = 'AnyFactor', members = c("factor", "list"))

#' The key slots used in the CellChat object are described below.
#'
#' @slot data.raw raw count data matrix
#' @slot data normalized data matrix for CellChat analysis (Genes should be in rows and cells in columns)
#' @slot data.signaling a subset of normalized matrix only containing signaling genes
#' @slot data.scale scaled data matrix
#' @slot data.project projected data
#' @slot images a list of spatial image objects
#' @slot net a three-dimensional array P (K×K×N), where K is the number of cell groups and N is the number of ligand-receptor pairs. Each row of P indicates the communication probability originating from the sender cell group to other cell groups.
#' @slot netP a three-dimensional array representing cel-cell communication networks on a signaling pathway level
#' @slot DB ligand-receptor interaction database used in the analysis (a subset of CellChatDB)
#' @slot LR a list of information related with ligand-receptor pairs
#' @slot meta data frame storing the information associated with each cell
#' @slot idents a factor defining the cell identity used for all analysis. It becomes a list for a merged CellChat object
#' @slot var.features A list: one element is a vector consisting of the identified over-expressed signaling genes; one element is a data frame returned from the differential expression analysis
#' @slot dr List of the reduced 2D coordinates, one per method, e.g., umap/tsne/dm
#' @slot options List of miscellaneous data, such as parameters used throughout analysis, and a indicator whether the CellChat object is a single or merged
#'
#' @exportClass CellChat
#' @importFrom Rcpp evalCpp
#' @importFrom methods setClass
# #' @useDynLib CellChat
CellChat <- methods::setClass("CellChat",
                                 slots = c(data.raw = 'AnyMatrix',
                                           data = 'AnyMatrix',
                                           data.signaling = "AnyMatrix",
                                           data.scale = "matrix",
                                           data.project = "AnyMatrix",
                                           images = "list",
                                           net = "list",
                                           netP = "list",
                                           meta = "data.frame",
                                           idents = "AnyFactor",
                                           DB = "list",
                                           LR = "list",
                                           var.features = "list",
                                           dr = "list",
                                           options = "list")
)
#' show method for CellChat
#'
#' @param CellChat object
#' @param show show the object
#' @param object object
#' @docType methods
#'
setMethod(f = "show", signature = "CellChat", definition = function(object) {
  if (object@options$mode == "single") {
    cat("An object of class", class(object), "created from a single dataset", "\n", nrow(object@data), "genes.\n",  ncol(object@data), "cells. \n")
  } else if (object@options$mode == "merged") {
    cat("An object of class", class(object), "created from a merged object with multiple datasets", "\n", nrow(object@data.signaling), "signaling genes.\n",  ncol(object@data.signaling), "cells. \n")
  }
  if (object@options$datatype == "RNA") {
    cat("CellChat analysis of single cell RNA-seq data! \n")
  } else {
    cat("CellChat analysis of", object@options$datatype, "data! The input spatial locations are \n")
    print(head(object@images$coordinates))
  }


  invisible(x = NULL)
})



#' Create a new CellChat object from a data matrix, Seurat or SingleCellExperiment object
#'
#' @param object a normalized (NOT count) data matrix (genes by cells), Seurat or SingleCellExperiment object
#' @param meta a data frame (rows are cells with rownames) consisting of cell information, which will be used for defining cell groups.
#' If input is a Seurat or SingleCellExperiment object, the meta data in the object will be used
#' @param group.by a char name of the variable in meta data, defining cell groups.
#' If input is a data matrix and group.by is NULL, the input `meta` should contain a column named 'labels',
#' If input is a Seurat or SingleCellExperiment object, USER must provide `group.by` to define the cell groups. e.g, group.by = "ident" for Seurat object
#' @param datatype By default datatype = "RNA"; when running CellChat on spatial imaging data, set datatype = "spatial" and input `scale.factors`
#'
#' @param coordinates a data matrix in which each row gives the spatial locations/coordinates of each cell/spot
#' @param scale.factors a list containing the scale factors and spot diameter for the full/high/low resolution images.
#'
#' USER must input this list when datatype = "spatial". scale.factors must contain an element named `spot.diameter`, which is the theoretical spot size; e.g., 10x Visium (spot.size = 65 microns), and another element named `spot`, which is the number of pixels that span the diameter of a theoretical spot size in the original, full-resolution image.
#'
#' For 10X visium, scale.factors are in the `scalefactors_json.json`. scale.factors$spot is the `spot.size.fullres `
#'
#' @param assay Assay to use when the input is a Seurat object. NB: The data in the `integrated` assay is not suitable for CellChat analysis because it contains negative values.
#' @param do.sparse whether use sparse format
#'
#' @return
#' @export
#' @importFrom methods as new
#' @examples
#' \dontrun{
#' Create a CellChat object from single-cell transcriptomics data
#' # Input is a data matrix
#' ## create a dataframe consisting of the cell labels
#' meta = data.frame(labels = cell.labels, row.names = names(cell.labels))
#' cellChat <- createCellChat(object = data.input, meta = meta, group.by = "labels")
#'
#' # input is a Seurat object
#' ## use the default cell identities of Seurat object
#' cellChat <- createCellChat(object = seurat.obj, group.by = "ident", assay = "RNA")
#' ## use other meta information as cell groups
#' cellChat <- createCellChat(object = seurat.obj, group.by = "seurat.clusters")
#'
#' # input is a SingleCellExperiment object
#' cellChat <- createCellChat(object = sce.obj, group.by = "sce.clusters")
#'
#' Create a CellChat object from spatial imaging data
#' # Input is a data matrix
#' cellChat <- createCellChat(object = data.input, meta = meta, group.by = "labels",
#'                            datatype = "spatial", coordinates = coordinates, scale.factors = scale.factors)
#'
#' # input is a Seurat object
#' cellChat <- createCellChat(object = seurat.obj, group.by = "ident", assay = "SCT",
#'                            datatype = "spatial", scale.factors = scale.factors)
#'
#' }
createCellChat <- function(object, meta = NULL, group.by = NULL,
                           datatype = c("RNA", "spatial"), coordinates = NULL, scale.factors = NULL,
                           assay = NULL, do.sparse = T) {
  datatype <- match.arg(datatype)
  # data matrix as input
  if (inherits(x = object, what = c("matrix", "Matrix", "dgCMatrix"))) {
    print("Create a CellChat object from a data matrix")
    data <- object
    if (is.null(group.by)) {
      group.by <- "labels"
    }
  }
  # Seurat object as input
  if (is(object,"Seurat")) {
    .error_if_no_Seurat()
    print("Create a CellChat object from a Seurat object")
    if (is.null(assay)) {
      assay = DefaultAssay(object)
      if (assay == "integrated") {
        warning("The data in the `integrated` assay is not suitable for CellChat analysis! Please use the `RNA` or `SCT` assay! ")
      }
      cat(paste0("The `data` slot in the default assay is used. The default assay is ", assay),'\n')
    }

    data <- Seurat::GetAssayData(object, assay = assay, slot = "data") # normalized data matrix
    if (min(data) < 0) {
      stop("The data matrix contains negative values. Please ensure the normalized data matrix is used.")
    }
    if (is.null(meta)) {
      cat("The `meta.data` slot in the Seurat object is used as cell meta information",'\n')
      meta <- object@meta.data
      meta$ident <- Seurat::Idents(object)
    }
    if (is.null(group.by)) {
      group.by <- "ident"
    }
    if (datatype %in% c("spatial")) {
      if (is.null(coordinates)) {
        coordinates <- GetTissueCoordinates(object, scale = NULL, cols = c("imagerow", "imagecol"))
        # scale.factors <- object@images[["slice1"]]@scale.factors
      }
    }


  }
  # SingleCellExperiment object as input
  if (is(object,"SingleCellExperiment")) {
    print("Create a CellChat object from a SingleCellExperiment object")
    if ("logcounts" %in% SummarizedExperiment::assayNames(object)) {
      cat("The `logcounts` assay is used",'\n')
      data <- SingleCellExperiment::logcounts(object)
    } else {
      stop("SingleCellExperiment object must contain an assay named `logcounts`")
    }
    if (is.null(meta)) {
      cat("The `colData` assay in the SingleCellExperiment object is used as cell meta information",'\n')
      meta <- as.data.frame(SingleCellExperiment::colData(object))
    }
    if (is.null(group.by)) {
      stop("`group.by` should be defined!")
    }
  }

  if (!inherits(x = data, what = c("dgCMatrix")) & do.sparse) {
    data <- as(data, "dgCMatrix")
  }

  if (!is.null(meta)) {
    if (inherits(x = meta, what = c("matrix", "Matrix"))) {
      meta <- as.data.frame(x = meta)
    }
    if (!is.data.frame(meta)) {
      stop("The input `meta` should be a data frame")
    }
    if (!identical(rownames(meta), colnames(data))) {
      cat("The cell barcodes in 'meta' is ", head(rownames(meta)),'\n')
      warning("The cell barcodes in 'meta' is different from those in the used data matrix.
              We now simply assign the colnames in the data matrix to the rownames of 'mata'!")
      rownames(meta) <- colnames(data)
    }
  } else {
    meta <- data.frame()
  }
  if (datatype %in% c("spatial")) {
    if (ncol(coordinates) == 2) {
      colnames(coordinates) <- c("x_cent","y_cent")
    } else {
      stop("Please check the input 'coordinates' and make sure it is a two column matrix.")
    }
    if (is.null(scale.factors) | !("spot.diameter" %in% names(scale.factors)) | !("spot" %in% names(scale.factors))) {
      stop("scale.factors with elements named `spot.diameter` and `spot` should be provided!")
    } else {
      images = list("coordinates" = coordinates,
                    "scale.factors" = scale.factors)
    }
    cat("Create a CellChat object from spatial imaging data...",'\n')
  } else {
    images <- list()
  }

  object <- methods::new(Class = "CellChat",
                         data = data,
                         images = images,
                         meta = meta)

  if (!is.null(meta) & nrow(meta) > 0) {
    cat("Set cell identities for the new CellChat object", '\n')
    if (!(group.by %in% colnames(meta))) {
      stop("The 'group.by' is not a column name in the `meta`, which will be used for cell grouping.")
    }
    object <- setIdent(object, ident.use = group.by) # set "labels" as default cell identity
    cat("The cell groups used for CellChat analysis are ", levels(object@idents), '\n')
  }
  object@options$mode <- "single"
  object@options$datatype <- datatype
  return(object)
}


#' Merge CellChat objects
#'
#' @param object.list  A list of multiple CellChat objects
#' @param add.names A vector containing the name of each dataset
#' @param merge.data whether merging the data for ALL genes. Default only merges the data of signaling genes
#' @param cell.prefix whether prefix cell names
#' @importFrom methods slot new
#'
#' @return
#' @export
#'
#' @examples
mergeCellChat <- function(object.list, add.names = NULL, merge.data = FALSE, cell.prefix = FALSE) {
  if (is.null(add.names)) {
    add.names <- paste("Dataset",1:length(object.list),sep = "_")
  }
  slot.name <- c("net", "netP", "idents" ,"LR", "var.features", "images")
  slot.combined <- vector("list", length(slot.name))
  names(slot.combined) <- slot.name
  for (i in 1:length(slot.name)) {
    object.slot <- vector("list", length(object.list))
    for (j in 1:length(object.list)) {
      object.slot[[j]] <- slot(object.list[[j]], slot.name[i])
    }
    slot.combined[[i]] <- object.slot
    names(slot.combined[[i]]) <- add.names
  }

  if (cell.prefix) {
    warning("Prefix cell names!")
    for (i in 1:length(object.list)) {colnames(object.list[[i]]@data) <- paste(colnames(object.list[[i]]@data), add.names[i], sep = "_")}
  } else {
    cell.names <- c()
    for (i in 1:length(object.list)) {
      cell.names <- c(cell.names, colnames(object.list[[i]]@data))
    }
    if (sum(duplicated(cell.names))) {
      stop("Duplicated cell names were detected across datasets!! Please set cell.prefix = TRUE")
    }
  }

  meta.use <- colnames(object.list[[1]]@meta)
  for (i in 2:length(object.list)) {
    meta.use <- meta.use[meta.use %in% colnames(object.list[[i]]@meta)]
  }

  dataset.name <- c()
  cell.names <- c()
  meta.joint <- data.frame()
  for (i in 1:length(object.list)) {
    dataset.name <- c(dataset.name, rep(add.names[i], length(colnames(object.list[[i]]@data))))
    cell.names <- c(cell.names, colnames(object.list[[i]]@data))
    meta.joint <- rbind(meta.joint, object.list[[i]]@meta[ , meta.use, drop = FALSE])
  }
  if (!identical(rownames(meta.joint), cell.names)) {
    cat("The cell barcodes in merged 'meta' is ", head(rownames(meta.joint)),'\n')
    warning("The cell barcodes in merged 'meta' is different from those in the used data matrix.
              We now simply assign the colnames in the data matrix to the rownames of merged 'mata'!")
    rownames(meta.joint) <- cell.names
  }

  #dataset.name <- data.frame(dataset.name = dataset.name, row.names = cell.names)
  meta.joint$datasets <- factor(dataset.name, levels = add.names)

  genes.use <- rownames(object.list[[1]]@data)
  for (i in 2:length(object.list)) {
    genes.use <- genes.use[genes.use %in% rownames(object.list[[i]]@data)]
  }
  data.joint <- c()
  for (i in 1:length(object.list)) {
    data.joint <- cbind(data.joint, object.list[[i]]@data[genes.use, ])
  }
  gene.signaling.joint = unique(unlist(lapply(object.list, function(x) rownames(x@data.signaling))))
  data.signaling.joint <- data.joint[rownames(data.joint) %in% gene.signaling.joint, ]

  idents.joint <- c()
  idents.levels <- c()
  for (i in 1:length(object.list)) {
    idents.joint <- c(idents.joint, as.character(object.list[[i]]@idents))
    idents.levels <- union(idents.levels, levels(object.list[[i]]@idents))
  }
  names(idents.joint) <- cell.names
  idents.joint <- factor(idents.joint, levels = idents.levels)
  slot.combined$idents$joint <- idents.joint

  if (merge.data) {
    message("Merge the following slots: 'data','data.signaling','images','net', 'netP','meta', 'idents', 'var.features', 'DB', and 'LR'.")
    merged.object <- methods::new(
      Class = "CellChat",
      data = data.joint,
      data.signaling = data.signaling.joint,
      images = slot.combined$images,
      net = slot.combined$net,
      netP = slot.combined$netP,
      meta = meta.joint,
      idents = slot.combined$idents,
      var.features = slot.combined$var.features,
      LR = slot.combined$LR,
      DB = object.list[[1]]@DB)
  } else {
    message("Merge the following slots: 'data.signaling','images','net', 'netP','meta', 'idents', 'var.features' , 'DB', and 'LR'.")
    merged.object <- methods::new(
      Class = "CellChat",
      data.signaling = data.signaling.joint,
      images = slot.combined$images,
      net = slot.combined$net,
      netP = slot.combined$netP,
      meta = meta.joint,
      idents = slot.combined$idents,
      var.features = slot.combined$var.features,
      LR = slot.combined$LR,
      DB = object.list[[1]]@DB)
  }
  merged.object@options$mode <- "merged"

  datatype.joint <- c()
  for (j in 1:length(object.list)) {
    datatype.joint <- union(datatype.joint, slot(object.list[[j]], "options")$datatype)
  }
  if (length(datatype.joint) == 1){
    merged.object@options$datatype <- datatype.joint
  } else {
    message("The data types in these objects are  ", datatype.joint,'\n')
    stop("Comparison analysis is not suggested for different types of data.")
  }
  return(merged.object)
}



#' Update a single CellChat object
#'
#' Update a single previously calculated CellChat object (version < 1.6.0)
#'
#' version < 0.5.0: `object@var.features` is now `object@var.features$features`; `object@net$sum` is now `object@net$weight` if `aggregateNet` has been run.
#'
#' version 1.6.0: a `object@images` slot is added and `datatype` is added in `object@options$datatype`
#'
#' @param object CellChat object
#'
#' @return a updated CellChat object
#' @export
#'
updateCellChat <- function(object) {

  if (is.character(object@var.features)) {
    message("Update slot 'var.features' from a vector to a list")
    var.features.new <- list(features = object@var.features)
  } else {
    var.features.new <- object@var.features
  }
  if ("sum" %in% names(object@net)) {
    net <- object@net
    net$weight <- net$sum
  } else {
    net <- object@net
  }
  if (!("mode" %in% names(object@options))) {
    object@options$mode <- "single"
  }
  if (!("datatype" %in% names(object@options))) {
    object@options$datatype <- "RNA"
    images = list()
  } else {
    images = object@images
  }
  object.new <- methods::new(
    Class = "CellChat",
    data.raw = object@data.raw,
    data = object@data,
    data.signaling = object@data.signaling,
    data.scale = object@data.scale,
    data.project = object@data.project,
    images = images,
    net = net,
    netP = object@netP,
    meta = object@meta,
    idents = object@idents,
    DB = object@DB,
    LR = object@LR,
    var.features = var.features.new,
    dr = object@dr,
    options = object@options
  )
  return(object.new)
}

#' Update a CellChat object by lifting up the cell groups to the same cell labels across all datasets
#'
#' This function is useful when comparing inferred communications across different datasets with different cellular compositions
#'
#' @param object A single or merged CellChat object
#' @param group.new A char vector giving the cell labels to lift up. The order of cell labels in the vector will be used for setting the new cell identity.
#'
#'  If the input is a merged CellChat object and group.new = NULL, it will use the cell labels from one dataset with the maximum number of cell groups
#'
#'  If the input is a single CellChat object, `group.new` must be defined.
#'
#' @return a updated CellChat object
#'
#' @export
#'
liftCellChat <- function(object, group.new = NULL) {
  if (object@options$mode == "merged") {
    idents <- object@idents[1:(length(object@idents)-1)]
    if (is.null(group.new)) {
      group.max.all <- unique(unlist(sapply(idents, levels)))
      group.num <- sapply(idents, nlevels)
      group.num.max <- max(group.num)
      group.max <- levels(idents[[which(group.num == group.num.max)]])
      if (length(group.max) != length(group.max.all)) {
        stop("CellChat object cannot lift up due to the missing cell groups in any dataset. Please define the parameter `group.new`!")
      }
    } else {
      group.max <- group.new
      group.num.max <- length(group.new)
    }
    message(paste0("The CellChat object will be lifted up using the cell labels ", paste(group.max, collapse=", ")))
    for (i in 1:length(idents)) {
      cat("Update slots object@net, object@netP, object@idents in dataset ", names(object@idents)[i],'\n')
      # cat("Update slot object@net...", '\n')
      net <- object@net[[i]]
      group.i <- levels(idents[[i]])
      # group.existing <- group.max[group.max %in% group.i]
      group.existing <- group.i[group.i %in% group.max]
      group.existing.index <- which(group.max %in% group.existing)
      for (net.j in names(net)) {
        values <- net[[net.j]]
        if (net.j %in% c("prob","pval")) {
          values.new <- array(data = 0, dim = c(group.num.max, group.num.max, dim(values)[3]),
                              dimnames = list(group.max, group.max, dimnames(values)[[3]]))
          values.new[group.existing.index, group.existing.index, ] <- values
          net[[net.j]] <- values.new
        }
        if (net.j %in% c("count","sum","weight")) {
          values.new <- array(data = 0, dim = c(group.num.max, group.num.max),
                              dimnames = list(group.max, group.max))
          values.new[group.existing.index, group.existing.index] <- values
          net[[net.j]] <- values.new
        }
        if (net.j %in% c("pairwiseRank")) {
          for (k in 1:length(values)) {
            values.new1 <- vector("list", group.num.max)
            values.new1[group.existing.index] <- values[[k]]
            temp <- values[[k]][[1]]
            temp$prob  <- 0; temp$pval <- 1
            for (kk in setdiff(1:group.num.max, group.existing.index)) {
              values.new1[[kk]] <- temp
            }
            names(values.new1) <- group.max
            values[[k]] <- values.new1
          }
          values.new <- vector("list", group.num.max)
          values.new[group.existing.index] <- values
          temp <- lapply(values.new1, function(x) {
            x$prob  <- 0; x$pval <- 1
            return(x)
          })
          for (kk in setdiff(1:group.num.max, group.existing.index)) {
            values.new[[kk]] <- temp
          }
          names(values.new) <- group.max
        }
        net[[net.j]] <- values.new
      }
      object@net[[i]] <- net

      # cat("Update slot object@netP...", '\n')
      netP <- object@netP[[i]]
      for (netP.j in names(netP)) {
        values <- netP[[netP.j]]
        if (netP.j %in% c("pathways")) {
          values.new <- values
          netP[[netP.j]] <- values.new
        }
        if (netP.j %in% c("prob")) {
          values.new <- array(data = 0, dim = c(group.num.max, group.num.max, dim(values)[3]),
                              dimnames = list(group.max, group.max, dimnames(values)[[3]]))
          values.new[group.existing.index, group.existing.index, ] <- values
          netP[[netP.j]] <- values.new
        }
        if (netP.j %in% c("centr")) {
          for (k in 1:length(values)) {
            values.new <- lapply(values, function(x) {
              values.new2 <- lapply(x, function(x) {
                values.new1 = as.vector(matrix(0, nrow = 1, ncol = group.num.max))
                values.new1[group.existing.index] <- x
                names(values.new1) <- group.max
                return(values.new1)
              })
              names(values.new2) <- names(x)
              return(values.new2)
            })
            names(values.new) <- names(values)
          }
          netP[[netP.j]] <- values.new
        }

      }
      object@netP[[i]] <- netP
      # cat("Update slot object@idents...", '\n')
      # idents[[i]] <- factor(group.max, levels = group.max)
      idents[[i]] <- factor(idents[[i]], levels = group.max)
    }
    object@idents[1:(length(object@idents)-1)] <- idents
  } else {
    if (is.null(group.new)) {
      stop("Please define the parameter `group.new`!")
    } else {
      group.max <- as.character(group.new)
      group.num.max <- length(group.new)
      message(paste0("The CellChat object will be lifted up using the cell labels ", paste(group.max, collapse=", ")))
    }
    cat("Update slots object@net, object@netP, object@idents in a single dataset...", '\n')
    # cat("Update slot object@net...", '\n')
    net <- object@net
    idents <- object@idents
    group.i <- levels(idents)
    # group.existing <- group.max[group.max %in% group.i]
    group.existing <- group.i[group.i %in% group.max]
    group.existing.index <- which(group.max %in% group.existing)
    for (net.j in names(net)) {
      values <- net[[net.j]]
      if (net.j %in% c("prob","pval")) {
        values.new <- array(data = 0, dim = c(group.num.max, group.num.max, dim(values)[3]),
                            dimnames = list(group.max, group.max, dimnames(values)[[3]]))
        values.new[group.existing.index, group.existing.index, ] <- values
        net[[net.j]] <- values.new
      }
      if (net.j %in% c("count","sum","weight")) {
        values.new <- array(data = 0, dim = c(group.num.max, group.num.max),
                            dimnames = list(group.max, group.max))
        values.new[group.existing.index, group.existing.index] <- values
        net[[net.j]] <- values.new
      }
      if (net.j %in% c("pairwiseRank")) {
        for (k in 1:length(values)) {
          values.new1 <- vector("list", group.num.max)
          values.new1[group.existing.index] <- values[[k]]
          temp <- values[[k]][[1]]
          temp$prob  <- 0; temp$pval <- 1
          for (kk in setdiff(1:group.num.max, group.existing.index)) {
            values.new1[[kk]] <- temp
          }
          names(values.new1) <- group.max
          values[[k]] <- values.new1
        }
        values.new <- vector("list", group.num.max)
        values.new[group.existing.index] <- values
        temp <- lapply(values.new1, function(x) {
          x$prob  <- 0; x$pval <- 1
          return(x)
        })
        for (kk in setdiff(1:group.num.max, group.existing.index)) {
          values.new[[kk]] <- temp
        }
        names(values.new) <- group.max
      }
      net[[net.j]] <- values.new
    }
    object@net <- net


    # cat("Update slot object@netP...", '\n')
    netP <- object@netP
    for (netP.j in names(netP)) {
      values <- netP[[netP.j]]
      if (netP.j %in% c("pathways")) {
        values.new <- values
        netP[[netP.j]] <- values.new
      }
      if (netP.j %in% c("prob")) {
        values.new <- array(data = 0, dim = c(group.num.max, group.num.max, dim(values)[3]),
                            dimnames = list(group.max, group.max, dimnames(values)[[3]]))
        values.new[group.existing.index, group.existing.index, ] <- values
        netP[[netP.j]] <- values.new
      }
      if (netP.j %in% c("centr")) {
        for (k in 1:length(values)) {
          values.new <- lapply(values, function(x) {
            values.new2 <- lapply(x, function(x) {
              values.new1 = as.vector(matrix(0, nrow = 1, ncol = group.num.max))
              values.new1[group.existing.index] <- x
              names(values.new1) <- group.max
              return(values.new1)
            })
            names(values.new2) <- names(x)
            return(values.new2)
          })
          names(values.new) <- names(values)
        }
      }
      netP[[netP.j]] <- values.new
    }
    object@netP <- netP

    # cat("Update slot object@idents...", '\n')
    idents <- factor(idents, levels = group.max)
    object@idents <- idents
  }

  return(object)
}


#' Subset CellChat object using a portion of cells
#'
#' @param object  A CellChat object (either an object from a single dataset or a merged objects from multiple datasets)
#' @param cells.use a char vector giving the cell barcodes to subset. If cells.use = NULL, USER must define `idents.use`
#' @param idents.use a subset of cell groups used for analysis
#' @param group.by cell group information; default is `object@idents`; otherwise it should be one of the column names of the meta slot
#' @param invert whether invert the idents.use
#' @param thresh threshold of the p-value for determining significant interaction. A parameter as an input of the function `computeCommunProbPathway`
#' @importFrom methods slot new
#'
#' @return
#' @export
#'
subsetCellChat <- function(object, cells.use = NULL, idents.use = NULL, group.by = NULL, invert = FALSE, thresh = 0.05) {
  if (!is.null(idents.use)) {
    if (is.null(group.by)) {
      labels <- object@idents
      if (object@options$mode == "merged") {
        message("Use the joint cell labels from the merged CellChat object")
        labels <- object@idents$joint
      }
    } else {
      labels <- object@meta[[group.by]]
    }
    if (!is.factor(labels)) {
      labels <- factor(labels)
    }
    level.use0 <- levels(labels)
    level.use <- levels(labels)[levels(labels) %in% unique(labels)]

    if (invert) {
      level.use <- level.use[!(level.use %in% idents.use)]
    } else {
      level.use <- level.use[level.use %in% idents.use]
    }
    cells.use.index <- which(as.character(labels) %in% level.use)
    cells.use <- names(labels)[cells.use.index]
  } else if (!is.null(cells.use)) {
    labels <- object@idents
    if (object@options$mode == "merged") {
      message("Use the joint cell labels from the merged CellChat object")
      labels <- object@idents$joint
    }
    level.use0 <- levels(labels)
    level.use <- levels(labels)[levels(labels) %in% unique(as.character(labels[cells.use]))]
    cells.use.index <- which(names(labels) %in% cells.use)
  } else {
    stop("USER should define either `cells.use` or `idents.use`!")
  }
  cat("The subset of cell groups used for CellChat analysis are ", level.use, '\n')

  if (nrow(object@data) > 0) {
    data.subset <- object@data[, cells.use.index]
  } else {
    data.subset <- matrix(0, nrow = 0, ncol = 0)
  }
  if (nrow(object@data.project) > 0) {
    data.project.subset <- object@data.project[, cells.use.index]
  } else {
    data.project.subset <- matrix(0, nrow = 0, ncol = 0)
  }
  data.signaling.subset <- object@data.signaling[, cells.use.index]

  meta.subset <- object@meta[cells.use.index, , drop = FALSE]


  if (object@options$mode == "merged") {
    idents <- object@idents[1:(length(object@idents)-1)]
    group.existing <- level.use0[level.use0 %in% level.use]
    group.existing.index <- which(level.use0 %in% level.use)
    net.subset <- vector("list", length = length(object@net))
    netP.subset <- vector("list", length = length(object@netP))
    idents.subset <- vector("list", length = length(idents))
    names(net.subset) <- names(object@net)
    names(netP.subset) <- names(object@netP)
    names(idents.subset) <- names(object@idents[1:(length(object@idents)-1)])
    images.subset <- vector("list", length = length(idents))
    names(images.subset) <- names(object@idents[1:(length(object@idents)-1)])

    for (i in 1:length(idents)) {
      cat("Update slots object@images, object@net, object@netP, object@idents in dataset ", names(object@idents)[i],'\n')
      images <- object@images[[i]]
      for (images.j in names(images)) {
        values <- images[[images.j]]
        if (images.j %in% c("coordinates")) {
          values.new <- values[cells.use.index, ]
          images[[images.j]] <- values.new
        }
        if (images.j %in% c("distance")) {
          values.new <- values[group.existing.index, group.existing.index, ,drop = FALSE]
          images[[images.j]] <- values.new
        }
      }
      images.subset[[i]] <- images

      # cat("Update slot object@net...", '\n')
      net <- object@net[[i]]
      for (net.j in names(net)) {
        values <- net[[net.j]]
        if (net.j %in% c("prob","pval")) {
          values.new <- values[group.existing.index, group.existing.index, ]
          net[[net.j]] <- values.new
        }
        if (net.j %in% c("count","sum","weight")) {
          values.new <- values[group.existing.index, group.existing.index]
          net[[net.j]] <- values.new
        }
       # net[[net.j]] <- values.new
      }
      net.subset[[i]] <- net

      # cat("Update slot object@netP...", '\n')
      # netP <- object@netP[[i]]
      # for (netP.j in names(netP)) {
      #   values <- netP[[netP.j]]
      #   if (netP.j %in% c("pathways")) {
      #     values.new <- values
      #     netP[[netP.j]] <- values.new
      #   }
      #   if (netP.j %in% c("prob")) {
      #     values.new <- values[group.existing.index, group.existing.index, ]
      #     netP[[netP.j]] <- values.new
      #   }
      #   if (netP.j %in% c("centr")) {
      #     for (k in 1:length(values)) {
      #       values.new <- lapply(values, function(x) {
      #         values.new2 <- lapply(x, function(x) {
      #           values.new1 <- x[group.existing.index]
      #           names(values.new1) <- group.existing
      #           return(values.new1)
      #         })
      #         names(values.new2) <- names(x)
      #         return(values.new2)
      #       })
      #       names(values.new) <- names(values)
      #     }
      #   }
      #   netP[[netP.j]] <- values.new
      # }
      netP = computeCommunProbPathway(net = net.subset[[i]], pairLR.use = object@LR[[i]]$LRsig, thresh = thresh)
      netP$centr = netAnalysis_computeCentrality(net =  net.subset[[i]]$prob)
      netP.subset[[i]] <- netP
      idents.subset[[i]] <- idents[[i]][names(idents[[i]]) %in% cells.use]
      idents.subset[[i]] <- factor(idents.subset[[i]], levels = levels(idents[[i]])[levels(idents[[i]]) %in% level.use])
    }
    idents.subset$joint <- factor(object@idents$joint[cells.use.index], levels = level.use)

  } else {
    cat("Update slots object@images, object@net, object@netP in a single dataset...", '\n')

    group.existing <- level.use0[level.use0 %in% level.use]
    group.existing.index <- which(level.use0 %in% level.use)

    images <- object@images
    for (images.j in names(images)) {
      values <- images[[images.j]]
      if (images.j %in% c("coordinates")) {
        values.new <- values[cells.use.index, ]
        images[[images.j]] <- values.new
      }
      if (images.j %in% c("distance")) {
        values.new <- values[group.existing.index, group.existing.index, ,drop = FALSE]
        images[[images.j]] <- values.new
      }
    }
    images.subset <- images


    # cat("Update slot object@net...", '\n')
    net <- object@net
    for (net.j in names(net)) {
      values <- net[[net.j]]
      if (net.j %in% c("prob","pval")) {
        values.new <- values[group.existing.index, group.existing.index, ,drop = FALSE]
        net[[net.j]] <- values.new
      }
      if (net.j %in% c("count","sum","weight")) {
        values.new <- values[group.existing.index, group.existing.index, ,drop = FALSE]
        net[[net.j]] <- values.new
      }
    }
    net.subset <- net

    # cat("Update slot object@netP...", '\n')
    # netP <- object@netP
    # for (netP.j in names(netP)) {
    #   values <- netP[[netP.j]]
    #   if (netP.j %in% c("pathways")) {
    #     values.new <- values
    #     netP[[netP.j]] <- values.new
    #   }
    #   if (netP.j %in% c("prob")) {
    #     values.new <- values[group.existing.index, group.existing.index, ]
    #     netP[[netP.j]] <- values.new
    #   }
    #   if (netP.j %in% c("centr")) {
    #     for (k in 1:length(values)) {
    #       values.new <- lapply(values, function(x) {
    #         values.new2 <- lapply(x, function(x) {
    #           values.new1 <- x[group.existing.index]
    #           names(values.new1) <- group.existing
    #           return(values.new1)
    #         })
    #         names(values.new2) <- names(x)
    #         return(values.new2)
    #       })
    #       names(values.new) <- names(values)
    #     }
    #   }
    #   netP[[netP.j]] <- values.new
    # }
    netP = computeCommunProbPathway(net = net.subset, pairLR.use = object@LR$LRsig, thresh = thresh)
    netP$centr = netAnalysis_computeCentrality(net = net.subset$prob)
    netP.subset <- netP
    idents.subset <- object@idents[cells.use.index]
    idents.subset <- factor(idents.subset, levels = level.use)
  }


  object.subset <- methods::new(
    Class = "CellChat",
    data = data.subset,
    data.signaling = data.signaling.subset,
    data.project = data.project.subset,
    images = images.subset,
    net = net.subset,
    netP = netP.subset,
    meta = meta.subset,
    idents = idents.subset,
    var.features = object@var.features,
    LR = object@LR,
    DB = object@DB,
    options = object@options
  )
  return(object.subset)
}


