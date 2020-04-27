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

#' Downsampling single cell data
#'
#' @param X input data (samples in rows, features in columns)
#' @param percent the percent of data to sketch
#' @param dimPC the number of components to use
#' @importFrom reticulate import
#' @return
#' @export
#'
sketchData <- function(X, percent, dimPC = 50) {
  # pip install geosketch
  geosketch <- reticulate::import('geosketch')
  # Get top PCs
  X.pcs <- runPCA(X, dimPC = dimPC)
  # Sketch percent of data.
  sketch.size <- as.integer(percent*nrow(X))
  sketch.index <- geosketch$gs(X.pcs, sketch.size)
  sketch.index <- unlist(sketch.index) + 1
  return(sketch.index)
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
#' @return
#' @export
#'
#' @examples
setIdent <- function(object, ident.use = NULL, levels = NULL){
  object@idents <- as.factor(object@meta[[ident.use]])
  if (!is.null(levels)) {
    object@idents <- factor(object@idents, levels = levels)
  }
  return(object)
}



#' Subset the expression data of signaling genes in CellChatDB
#'
#' @param object CellChat object
#' @param features default = NULL: subset the expression data of signaling genes in CellChatDB
#'
#' @return
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

#' Identify over-expressed genes in each cell group
#'
#' @param object CellChat object
#' @param features a vector of features
#' @param thresh.p threshold of p-values
#' @importFrom stats wilcox.test
#' @importFrom future nbrOfWorkers
#' @importFrom future.apply future_sapply
#' @importFrom pbapply pbsapply
#' @importFrom stats p.adjust
#'
#' @return
#' @export
#'
#' @examples
identifyOverExpressedGenes <- function(object, features = NULL, thresh.p = 0.05) {
  data.use <- object@data.signaling
  if (is.null(features)) {
    features <- row.names(data.use)
  } else {
    features <- intersect(features, row.names(data.use))
  }

  data.use <- as.matrix(data.use[features, ])
  labels <- object@idents
  level.use <- levels(labels)
  level.use <- level.use[level.use %in% unique(labels)]
  numCluster <- length(level.use)

  my.sapply <- ifelse(
    test = future::nbrOfWorkers() == 1,
    yes = pbapply::pbsapply,
    no = future.apply::future_sapply
  )

  Pvalues <- matrix(nrow = length(features), ncol = numCluster)
  for (i in 1:numCluster) {
    cell.use1 <- which(labels %in% level.use[i])
    cell.use2 <- base::setdiff(1:length(labels), cell.use1)
    data1 <- data.use[, cell.use1, drop = FALSE]
    data2 <- data.use[, cell.use2, drop = FALSE]
    pvalues <- unlist(
      x = my.sapply(
        X = 1:nrow(data1),
        FUN = function(x) {
          return(wilcox.test(data1[x, ], data2[x, ], alternative = "greater")$p.value)
        }
      )
    )
    Pvalues[, i] <- pvalues
  }
  features.sig <- features[rowSums(Pvalues < thresh.p) > 0]
  object@var.features <- features.sig
  return(object)
}


#' Identify over-expressed ligand-receptor interactions
#'
#' @param object CellChat object
#' @param features a vector of features
#' @importFrom future nbrOfWorkers
#' @importFrom future.apply future_sapply
#' @importFrom pbapply pbsapply
#' @importFrom dplyr select
#'
#' @return
#' @export
#'
#' @examples
identifyOverExpressedInteractions <- function(object, features = NULL) {
  if (is.null(features)) {
    gene.use <- row.names(object@data.signaling)
  } else {
    gene.use <- intersect(features, row.names(object@data.signaling))
  }
  DB <- object@DB
  features.sig <- object@var.features
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
  return(object)
}


#' Project gene expression data onto a protein-protein interaction network
#' @param object   CellChat object
#' @param adjMatrix    adjacency matrix of protein-protein interaction network to use
#' @param alpha    numeric in [0,1]
#' @param normalizeAdjMatrix    how to normalize the adjacency matrix
#'                              possible values are 'rows' (in-degree)
#'                              and 'columns' (out-degree)
#' @return projected gene expression matrix
#' @examples
#' @export
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
  data_in_new_space = matrix(rep(0, length(new_features)*dim(gene_expression)[2]),nrow=length(new_features))
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
#' @param n.neighbors This determines the number of neighboring points used in
#' local approximations of manifold structure. Larger values will result in more
#' global structure being preserved at the loss of detailed local structure. In general this parameter should often be in the range 5 to 50.
#' @param n.components The dimension of the space to embed into.
#' @param distance This determines the choice of metric used to measure distance in the input space.
#' @param n.epochs the number of training epochs to be used in optimizing the low dimensional embedding. Larger values result in more accurate embeddings. If NULL is specified, a value will be selected based on the size of the input dataset (200 for large datasets, 500 for small).
#' @param learning.rate The initial learning rate for the embedding optimization.
#' @param min.dist This controls how tightly the embedding is allowed compress points together.
#' Larger values ensure embedded points are moreevenly distributed, while smaller values allow the
#' algorithm to optimise more accurately with regard to local structure. Sensible values are in the range 0.001 to 0.5.
#' @param spread he effective scale of embedded points. In combination with min.dist this determines how clustered/clumped the embedded points are.
#' @param set.op.mix.ratio Interpolate between (fuzzy) union and intersection as the set operation used to combine local fuzzy simplicial sets to obtain a global fuzzy simplicial sets.
#' @param local.connectivity The local connectivity required - i.e. the number of nearest neighbors
#' that should be assumed to be connected at a local level. The higher this value the more connected
#' the manifold becomes locally. In practice this should be not more than the local intrinsic dimension of the manifold.
#' @param repulsion.strength Weighting applied to negative samples in low dimensional embedding
#' optimization. Values higher than one will result in greater weight being given to negative samples.
#' @param negative.sample.rate The number of negative samples to select per positive sample in the
#' optimization process. Increasing this value will result in greater repulsive force being applied, greater optimization cost, but slightly more accuracy.
#' @param a More specific parameters controlling the embedding. If NULL, these values are set automatically as determined by min. dist and spread.
#' @param b More specific parameters controlling the embedding. If NULL, these values are set automatically as determined by min. dist and spread.
#' @param seed.use Set a random seed. By default, sets the seed to 42.
#' @param metric.kwds,angular.rp.forest,verbose other parameters used in UMAP
#' @import reticulate
#' @export
#'
runUMAP <- function(
  data.use,
  n.neighbors = 30L,
  n.components = 2L,
  distance = "correlation",
  n.epochs = NULL,
  learning.rate = 1.0,
  min.dist = 0.3,
  spread = 1.0,
  set.op.mix.ratio = 1.0,
  local.connectivity = 1L,
  repulsion.strength = 1,
  negative.sample.rate = 5,
  a = NULL,
  b = NULL,
  seed.use = 42L,
  metric.kwds = NULL,
  angular.rp.forest = FALSE,
  verbose = FALSE){
  if (!reticulate::py_module_available(module = 'umap')) {
    stop("Cannot find UMAP, please install through pip (e.g. pip install umap-learn or reticulate::py_install(packages = 'umap-learn')).")
  }
  set.seed(seed.use)
  reticulate::py_set_seed(seed.use)
  umap_import <- reticulate::import(module = "umap", delay_load = TRUE)
  umap <- umap_import$UMAP(
    n_neighbors = as.integer(n.neighbors),
    n_components = as.integer(n.components),
    metric = distance,
    n_epochs = n.epochs,
    learning_rate = learning.rate,
    min_dist = min.dist,
    spread = spread,
    set_op_mix_ratio = set.op.mix.ratio,
    local_connectivity = local.connectivity,
    repulsion_strength = repulsion.strength,
    negative_sample_rate = negative.sample.rate,
    a = a,
    b = b,
    metric_kwds = metric.kwds,
    angular_rp_forest = angular.rp.forest,
    verbose = verbose
  )
  Rumap <- umap$fit_transform
  umap_output <- Rumap(t(data.use))
  colnames(umap_output) <- paste0('UMAP', 1:ncol(umap_output))
  rownames(umap_output) <- colnames(data.use)
  return(umap_output)
}

