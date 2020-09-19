
#' Compute and visualize the contribution of each ligand-receptor pair in the overall signaling pathways
#'
#' @param object CellChat object
#' @param signaling a signaling pathway name
#' @param signaling.name alternative signaling pathway name to show on the plot
#' @param width the width of individual bar
#' @param vertex.receiver a numeric vector giving the index of the cell groups as targets in the first hierarchy plot
#' @param thresh threshold of the p-value for determining significant interaction
#' @param x.rotation rotation of x-label
#' @param title the title of the plot
#' @param font.size font size of the text
#' @param font.size.title font size of the title
#' @importFrom dplyr select
#' @importFrom ggplot2 ggplot geom_bar aes coord_flip scale_x_discrete element_text theme ggtitle
#' @importFrom cowplot ggdraw draw_label plot_grid
#'
#' @return
#' @export
#'
#' @examples
netAnalysis_contribution <- function(object, signaling, signaling.name = NULL, width = 0.1, vertex.receiver = NULL, thresh = 0.05,
                                     x.rotation = 0, title = "Contribution of each L-R pair",
                                     font.size = 10, font.size.title = 10) {
  pairLR <- searchPair(signaling = signaling, pairLR.use = object@LR$LRsig, key = "pathway_name", matching.exact = T, pair.only = T)
  pair.name.use = select(object@DB$interaction[rownames(pairLR),],"interaction_name_2")
  if (is.null(signaling.name)) {
    signaling.name <- signaling
  }
  net <- object@net
  pairLR.use.name <- dimnames(net$prob)[[3]]
  pairLR.name <- intersect(rownames(pairLR), pairLR.use.name)
  pairLR <- pairLR[pairLR.name, ]
  prob <- net$prob
  pval <- net$pval

  prob[pval > thresh] <- 0
  if (length(pairLR.name) > 1) {
    pairLR.name.use <- pairLR.name[apply(prob[,,pairLR.name], 3, sum) != 0]
  } else {
    pairLR.name.use <- pairLR.name[sum(prob[,,pairLR.name]) != 0]
  }


  if (length(pairLR.name.use) == 0) {
    stop(paste0('There is no significant communication of ', signaling.name))
  } else {
    pairLR <- pairLR[pairLR.name.use,]
  }

  prob <- prob[,,pairLR.name.use]

  if (length(dim(prob)) == 2) {
    prob <- replicate(1, prob, simplify="array")
    dimnames(prob)[3] <- pairLR.name.use
  }
  prob <-(prob-min(prob))/(max(prob)-min(prob))

  if (is.null(vertex.receiver)) {
    pSum <- apply(prob, 3, sum)
    pSum.max <- sum(prob)
    pSum <- pSum/pSum.max
    pSum[is.na(pSum)] <- 0
    y.lim <- max(pSum)

    pair.name <- unlist(dimnames(prob)[3])
    pair.name <- factor(pair.name, levels = unique(pair.name))
    if (!is.null(pair.name.use)) {
      pair.name <- pair.name.use[as.character(pair.name),1]
      pair.name <- factor(pair.name, levels = unique(pair.name))
    }
    mat <- pSum
    df1 <- data.frame(name = pair.name, contribution = mat)
    if(nrow(df1) < 10) {
      df2 <- data.frame(name = as.character(1:(10-nrow(df1))), contribution = rep(0, 10-nrow(df1)))
      df <- rbind(df1, df2)
    } else {
      df <- df1
    }
    gg <- ggplot(df, aes(x=name, y=contribution)) + geom_bar(stat="identity", width = 0.7) +
      theme_classic() + theme(axis.text.y = element_text(angle = x.rotation, hjust = 1,size=font.size, colour = 'black'), axis.text=element_text(size=font.size),
                              axis.title.y = element_text(size= font.size), axis.text.x = element_blank(), axis.ticks = element_blank()) +
      xlab("") + ylab("Relative contribution") + ylim(0,y.lim) + coord_flip() + theme(legend.position="none") +
      scale_x_discrete(limits = rev(levels(df$name)), labels = c(rep("", max(0, 10-nlevels(df1$name))),rev(levels(df1$name))))
    if (!is.null(title)) {
      gg <- gg + ggtitle(title)+ theme(plot.title = element_text(hjust = 0.5, size = font.size.title))
    }
    gg

  } else {
    pair.name <- factor(unlist(dimnames(prob)[3]), levels = unique(unlist(dimnames(prob)[3])))
    # show all the communications
    pSum <- apply(prob, 3, sum)
    pSum.max <- sum(prob)
    pSum <- pSum/pSum.max
    pSum[is.na(pSum)] <- 0
    y.lim <- max(pSum)

    df<- data.frame(name = pair.name, contribution = pSum)
    gg <- ggplot(df, aes(x=name, y=contribution)) + geom_bar(stat="identity",width = 0.2) +
      theme_classic() + theme(axis.text=element_text(size=10),axis.text.x = element_text(angle = x.rotation, hjust = 1,size=8),
                              axis.title.y = element_text(size=10)) +
      xlab("") + ylab("Relative contribution") + ylim(0,y.lim)+ ggtitle("All")+ theme(plot.title = element_text(hjust = 0.5))#+

    # show the communications in Hierarchy1
    if (dim(prob)[3] > 1) {
      pSum <- apply(prob[,vertex.receiver,], 3, sum)
    } else {
      pSum <- sum(prob[,vertex.receiver,])
    }

    pSum <- pSum/pSum.max
    pSum[is.na(pSum)] <- 0

    df<- data.frame(name = pair.name, contribution = pSum)
    gg1 <- ggplot(df, aes(x=name, y=contribution)) + geom_bar(stat="identity",width = 0.2) +
      theme_classic() + theme(axis.text=element_text(size=10),axis.text.x = element_text(angle = x.rotation, hjust = 1,size=8), axis.title.y = element_text(size=10)) +
      xlab("") + ylab("Relative contribution") + ylim(0,y.lim)+ ggtitle("Hierarchy1") + theme(plot.title = element_text(hjust = 0.5))#+
    #scale_x_discrete(limits = c(0,1))

    # show the communications in Hierarchy2

    if (dim(prob)[3] > 1) {
      pSum <- apply(prob[,setdiff(1:dim(prob)[1],vertex.receiver),], 3, sum)
    } else {
      pSum <- sum(prob[,setdiff(1:dim(prob)[1],vertex.receiver),])
    }
    pSum <- pSum/pSum.max
    pSum[is.na(pSum)] <- 0

    df<- data.frame(name = pair.name, contribution = pSum)
    gg2 <- ggplot(df, aes(x=name, y=contribution)) + geom_bar(stat="identity", width=0.9) +
      theme_classic() + theme(axis.text=element_text(size=10),axis.text.x = element_text(angle = x.rotation, hjust = 1,size=8), axis.title.y = element_text(size=10)) +
      xlab("") + ylab("Relative contribution") + ylim(0,y.lim)+ ggtitle("Hierarchy2")+ theme(plot.title = element_text(hjust = 0.5))#+
    #scale_x_discrete(limits = c(0,1))
    title <- cowplot::ggdraw() + cowplot::draw_label(paste0("Contribution of each signaling in ", signaling.name, " pathway"), fontface='bold', size = 10)
    gg.combined <- cowplot::plot_grid(gg, gg1, gg2, nrow = 1)
    gg.combined <- cowplot::plot_grid(title, gg.combined, ncol = 1, rel_heights=c(0.1, 1))
    gg <- gg.combined
    gg
  }
    return(list(LRpair = as.character(pair.name), gg.obj = gg))
}


#' Identification of dominant senders, receivers, mediators and influencers in the intercellular communication network
#'
#' @param object CellChat object
#' @param slot.name the slot name of object that is used to compute centrality measures of signaling networks
#' @param net compute the centrality measures on a specific signaling network given by a 2 or 3 dimemsional array net
#' @param net.name a character vector giving the name of signaling networks
#' @importFrom future nbrOfWorkers
#' @importFrom methods slot
#' @importFrom future.apply future_sapply
#' @importFrom pbapply pbsapply
#'
#' @return
#' @export
#'
#' @examples
netAnalysis_signalingRole <- function(object, slot.name = "netP", net = NULL, net.name = NULL) {
  if (is.null(net)) {
    net = methods::slot(object, slot.name)$prob
  }
  if (is.null(net.name)) {
    net.name <- dimnames(net)[[3]]
  }
  if (length(dim(net)) == 3) {
    nrun <- dim(net)[3]
    my.sapply <- ifelse(
      test = future::nbrOfWorkers() == 1,
      yes = pbapply::pbsapply,
      no = future.apply::future_sapply
    )
    centr.all = my.sapply(
      X = 1:nrun,
      FUN = function(x) {
        net0 <- net[ , , x]
        return(computeCentralityLocal(net0))
      },
      simplify = FALSE
    )
  } else {
    centr.all <- as.list(computeCentralityLocal(net))
  }
  names(centr.all) <- net.name
  slot(object, slot.name)$centr <- centr.all
  return(object)
}



#' Compute Centrality measures for a signaling network
#'
#' @param net compute the centrality measures on a specific signaling network given by a 2 or 3 dimemsional array net
#' @importFrom igraph graph_from_adjacency_matrix strength hub_score authority_score eigen_centrality page_rank betweenness
#' @importFrom sna flowbet infocent
#'
#' @return
computeCentralityLocal <- function(net) {
  centr <- vector("list")
  G <- igraph::graph_from_adjacency_matrix(net, mode = "directed", weighted = T)
  centr$outdeg <- igraph::strength(G, mode="out")
  centr$indeg <- igraph::strength(G, mode="in")
  centr$hub <- igraph::hub_score(G)$vector
  centr$authority <- igraph::authority_score(G)$vector # A node has high authority when it is linked by many other nodes that are linking many other nodes.
  centr$eigen <- igraph::eigen_centrality(G)$vector # A measure of influence in the network that takes into account second-order connections
  centr$page_rank <- igraph::page_rank(G)$vector
  E(G)$weight <- 1/E(G)$weight
  centr$betweenness <- igraph::betweenness(G)
  #centr$flowbet <- try(sna::flowbet(net)) # a measure of its role as a gatekeeper for the flow of communication between any two cells; the total maximum flow (aggregated across all pairs of third parties) mediated by v.
  #centr$info <- try(sna::infocent(net)) # actors with higher information centrality are predicted to have greater control over the flow of information within a network; highly information-central individuals tend to have a large number of short paths to many others within the social structure.
  centr$flowbet <- tryCatch({
    sna::flowbet(net)
  }, error = function(e) {
    as.vector(matrix(0, nrow = nrow(net), ncol = 1))
  })
  centr$info <- tryCatch({
    sna::infocent(net, diag = T, rescale = T, cmode = "lower")
  }, error = function(e) {
    as.vector(matrix(0, nrow = nrow(net), ncol = 1))
  })
  return(centr)
}



#' Identification of major signals for specific cell groups and general communication patterns
#'
#' @param object CellChat object
#' @param slot.name the slot name of object that is used to compute centrality measures of signaling networks
#' @param pattern "outgoing" or "incoming"
#' @param k the number of patterns
#' @param k.range a range of the number of patterns
#' @param heatmap.show whether showing heatmap
#' @param color.use the character vector defining the color of each cell group
#' @param color.heatmap a color name in brewer.pal
#' @param title.legend the title of legend in heatmap
#' @param width width of heatmap
#' @param height height of heatmap
#' @param font.size fontsize in heatmap
#' @importFrom methods slot
#' @importFrom NMF nmfEstimateRank nmf
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#' @importFrom ComplexHeatmap Heatmap HeatmapAnnotation draw
#' @importFrom stats setNames
#' @importFrom grid grid.grabExpr grid.newpage pushViewport grid.draw unit gpar viewport popViewport
#'
#' @return
#' @export
#'
#' @examples

identifyCommunicationPatterns <- function(object, slot.name = "netP", pattern = c("outgoing","incoming"), k = NULL, k.range = seq(2,10), heatmap.show = TRUE,
                                          color.use = NULL, color.heatmap = "Spectral", title.legend = "Contributions",
                                          width = 4, height = 6, font.size = 8) {
  pattern <- match.arg(pattern)
  prob <- methods::slot(object, slot.name)$prob
  if (pattern == "outgoing") {
    data_sender <- apply(prob, c(1,3), sum)
    data_sender = sweep(data_sender, 2L, apply(data_sender, 2, function(x) max(x, na.rm = TRUE)), '/', check.margin = FALSE)
    data0 = as.matrix(data_sender)
  } else if (pattern == "incoming") {
    data_receiver <- apply(prob, c(2,3), sum)
    data_receiver = sweep(data_receiver, 2L, apply(data_receiver, 2, function(x) max(x, na.rm = TRUE)), '/', check.margin = FALSE)
    data0 = as.matrix(data_receiver)
  }
  options(warn = -1)
  data <- data0
  data <- data[rowSums(data)!=0,]
  if (is.null(k)) {
    res <- NMF::nmfEstimateRank(data, range = k.range, method = 'lee', nrun=10)
    # It is possible to specify that the currently registered backend should be used, by setting argument .pbackend=NULL.
    idx <- min(which(abs(diff(res$measures$cophenetic)) > 1*sd(res$measures$cophenetic)))
    k = res$measures$rank[idx]
    cat("Plot the cophenetic measure over a given range of rank (i.e., number of patterns): ",'\n')
    pdf(file = paste0("cophenetic measure_",k,"_.pdf"), width = 4, height = 3)
    plot(res, 'cophenetic',xlab = 'Number of patterns', ylab = 'Cophenetic', main = " ")
    dev.off()
    cat("Based on the cophenetic measure, the optimal number of patterns is", k,'\n')
  }

  outs_NMF <- NMF::nmf(data, rank = k, method = 'lee', seed = 'nndsvd')
  W <- scaleMat(outs_NMF@fit@W, 'r1')
  H <- scaleMat(outs_NMF@fit@H, 'c1')
  colnames(W) <- paste0("Pattern ", seq(1,ncol(W))); rownames(H) <- paste0("Pattern ", seq(1,nrow(H)));
  if (heatmap.show) {
    net <- W
    if (is.null(color.use)) {
      color.use <- scPalette(length(rownames(net)))
    }
    color.heatmap = grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n = 9, name = color.heatmap)))(100)

    df<- data.frame(group = rownames(net)); rownames(df) <- rownames(net)
    cell.cols.assigned <- setNames(color.use, unique(as.character(df$group)))
    row_annotation <- HeatmapAnnotation(df = df, col = list(group = cell.cols.assigned),which = "row",
                                        show_legend = FALSE, show_annotation_name = FALSE,
                                        simple_anno_size = grid::unit(0.2, "cm"))

    ht1 = Heatmap(net, col = color.heatmap, na_col = "white", name = "Contribution",
                  left_annotation = row_annotation,
                  cluster_rows = T,cluster_columns = F,clustering_method_rows = "average",
                  row_names_side = "left",row_names_rot = 0,row_names_gp = gpar(fontsize = font.size),column_names_gp = gpar(fontsize = font.size),
                  width = unit(width, "cm"), height = unit(height, "cm"),
                  show_heatmap_legend = F,
                  column_title = "Cell patterns",column_title_gp = gpar(fontsize = 10)
    )

    df<- data.frame(group = rownames(net)); rownames(df) <- rownames(net)
    cell.cols.assigned <- setNames(color.use, unique(as.character(df$group)))
    row_annotation <- HeatmapAnnotation(df = df, col = list(group = cell.cols.assigned),which = "row",
                                        show_legend = FALSE, show_annotation_name = FALSE,
                                        simple_anno_size = grid::unit(0.2, "cm"))

    net <- t(H)

    ht2 = Heatmap(net, col = color.heatmap, na_col = "white", name = "Contribution",
                  cluster_rows = T,cluster_columns = F,clustering_method_rows = "average",
                  row_names_side = "left",row_names_rot = 0,row_names_gp = gpar(fontsize = font.size),column_names_gp = gpar(fontsize = font.size),
                  width = unit(width, "cm"), height = unit(height, "cm"),
                  column_title = "Communication patterns",column_title_gp = gpar(fontsize = 10),
                  heatmap_legend_param = list(title = title.legend, title_gp = gpar(fontsize = 8, fontface = "plain"),title_position = "leftcenter-rot",
                                              border = NA, at = c(round(min(net, na.rm = T), digits = 1), round(max(net, na.rm = T), digits = 1)),
                                              legend_height = unit(20, "mm"),labels_gp = gpar(fontsize = 6),grid_width = unit(2, "mm"))
    )

    gb_ht1 = grid.grabExpr(draw(ht1))
    gb_ht2 = grid.grabExpr(draw(ht2))
    #grid.newpage()
    pushViewport(viewport(x = 0.1, y = 0.1, width = 0.2, height = 0.5, just = c("left", "bottom")))
    grid.draw(gb_ht1)
    popViewport()

    pushViewport(viewport(x = 0.6, y = 0.1, width = 0.2, height = 0.5, just = c("left", "bottom")))
    grid.draw(gb_ht2)
    popViewport()

  }

  data_W <- as.data.frame(as.table(W)); colnames(data_W) <- c("CellGroup","Pattern","Contribution")
  data_H <- as.data.frame(as.table(H)); colnames(data_H) <- c("Pattern","Signaling","Contribution")

  res.pattern = list("cell" = data_W, "signaling" = data_H)
  methods::slot(object, slot.name)$pattern[[pattern]] <- list(data = data0, pattern = res.pattern)
  return(object)
}


#' Compute signaling network similarity for any pair of signaling networks
#'
#' @param object CellChat object
#' @param slot.name the slot name of object that is used to compute centrality measures of signaling networks
#' @param type "functional","structural"
#' @param k the number of nearest neighbors
#' @param thresh the fraction (0 to 0.25) of interactions to be trimmed before computing network similarity
#' @importFrom methods slot

#'
#' @return
#' @export
#'
#' @examples
computeNetSimilarity <- function(object, slot.name = "netP", type = c("functional","structural"), k = NULL, thresh = NULL) {
  type <- match.arg(type)
  prob = methods::slot(object, slot.name)$prob
  if (is.null(k)) {
    if (dim(prob)[3] <= 25) {
      k <- ceiling(sqrt(dim(prob)[3]))
    } else {
      k <- ceiling(sqrt(dim(prob)[3])) + 1
    }

  }
  if (!is.null(thresh)) {
    prob[prob < quantile(c(prob[prob != 0]), thresh)] <- 0
  }
  if (type == "functional") {
    # compute the functional similarity
    D_signalings <- matrix(0, nrow = dim(prob)[3], ncol = dim(prob)[3])
    S2 <- D_signalings; S3 <- D_signalings;
    for (i in 1:(dim(prob)[3]-1)) {
      for (j in (i+1):dim(prob)[3]) {
        Gi <- (prob[ , ,i] > 0)*1
        Gj <- (prob[ , ,j] > 0)*1
        S3[i,j] <- sum(Gi * Gj)/sum(Gi+Gj-Gi*Gj,na.rm=TRUE)
      }
    }
    # define the similarity matrix
    S3[is.na(S3)] <- 0; S3 <- S3 + t(S3); diag(S3) <- 1
    # S_signalings <- S1 *S2
    S_signalings <- S3
  } else if (type == "structural") {
    # compute the structure distance
    D_signalings <- matrix(0, nrow = dim(prob)[3], ncol = dim(prob)[3])
    for (i in 1:(dim(prob)[3]-1)) {
      for (j in (i+1):dim(prob)[3]) {
        Gi <- (prob[ , ,i] > 0)*1
        Gj <- (prob[ , ,j] > 0)*1
        D_signalings[i,j] <- computeNetD_structure(Gi,Gj)
      }
    }
    # define the structure similarity matrix
    D_signalings[is.infinite(D_signalings)] <- 0
    D_signalings[is.na(D_signalings)] <- 0
    D_signalings <- D_signalings + t(D_signalings)
    S_signalings <- 1-D_signalings
  }

  # smooth the similarity matrix using SNN
  SNN <- buildSNN(S_signalings, k = k, prune.SNN = 1/15)
  Similarity <- as.matrix(S_signalings*SNN)
  rownames(Similarity) <- dimnames(prob)[[3]]
  colnames(Similarity) <- dimnames(prob)[[3]]

  methods::slot(object, slot.name)$similarity[[type]]$matrix <- Similarity
  return(object)
}



#' Compute signaling network similarity for a pair of datasets
#'
#' @param object A merged CellChat object
#' @param slot.name the slot name of object that is used to compute centrality measures of signaling networks
#' @param type "functional","structural"
#' @param k the number of nearest neighbors
#' @param thresh the fraction (0 to 0.25) of interactions to be trimmed before computing network similarity
#' @importFrom methods slot
#'
#' @return
#' @export
#'
#' @examples
computeNetSimilarityPairwise <- function(object, slot.name = "netP", type = c("functional","structural"), k = NULL, thresh = NULL) {
  type <- match.arg(type)
  net <- list()
  signalingAll <- c()
  for (i in 1:length(setdiff(names(methods::slot(object, slot.name)), "similarity"))) {
    object.net <- methods::slot(object, slot.name)[[i]]
    net[[i]] = object.net$prob
    signalingAll <- c(signalingAll, dimnames(net[[i]])[[3]])
  }
  net.dim <- sapply(net, dim)[3,]
  nnet <- sum(net.dim)
  position <- cumsum(net.dim); position <- c(0,position)

  if (is.null(k)) {
    if (nnet <= 25) {
      k <- ceiling(sqrt(nnet))
    } else {
      k <- ceiling(sqrt(nnet)) + 1
    }

  }
  if (!is.null(thresh)) {
    for (i in 1:length(net)) {
      neti <- net[[i]]
      neti[neti < quantile(c(neti[neti != 0]), thresh)] <- 0
      net[[i]] <- neti
    }
  }
  if (type == "functional") {
    # compute the functional similarity
    S3 <- matrix(0, nrow = nnet, ncol = nnet)
    for (i in 1:nnet) {
      for (j in 1:nnet) {
        idx.i <- which(position - i >= 0)[1]
        idx.j <- which(position - j >= 0)[1]
        net.i <- net[[idx.i-1]]
        net.j <- net[[idx.j-1]]
        Gi <- (net.i[ , ,i-position[idx.i-1]] > 0)*1
        Gj <- (net.j[ , ,j-position[idx.j-1]] > 0)*1
        S3[i,j] <- sum(Gi * Gj)/sum(Gi+Gj-Gi*Gj,na.rm=TRUE)
      }
    }

    # define the similarity matrix
    S3[is.na(S3)] <- 0;  diag(S3) <- 1
    S_signalings <- S3
  } else if (type == "structural") {
    # compute the structure distance
    D_signalings <- matrix(0, nrow = nnet, ncol = nnet)
    for (i in 1:nnet) {
      for (j in 1:nnet) {
        idx.i <- which(position - i >= 0)[1]
        idx.j <- which(position - j >= 0)[1]
        net.i <- net[[idx.i-1]]
        net.j <- net[[idx.j-1]]
        Gi <- (net.i[ , ,i-position[idx.i-1]] > 0)*1
        Gj <- (net.j[ , ,j-position[idx.j-1]] > 0)*1
        D_signalings[i,j] <- computeNetD_structure(Gi,Gj)
      }
    }
    # define the structure similarity matrix
    D_signalings[is.infinite(D_signalings)] <- 0
    D_signalings[is.na(D_signalings)] <- 0
    S_signalings <- 1-D_signalings
  }
  # smooth the similarity matrix using SNN
  SNN <- buildSNN(S_signalings, k = k, prune.SNN = 1/15)
  Similarity <- as.matrix(S_signalings*SNN)
  rownames(Similarity) <- signalingAll
  colnames(Similarity) <- rownames(Similarity)

  methods::slot(object, slot.name)$similarity[[type]]$matrix <- Similarity
  return(object)
}


#' Manifold learning of the signaling networks based on their similarity
#'
#' @param object CellChat object
#' @param slot.name the slot name of object that is used to compute centrality measures of signaling networks
#' @param type "functional","structural"
#' @param k the number of nearest neighbors in running umap
#' @param pathway.remove a range of the number of patterns
#' @importFrom methods slot
#' @return
#' @export
#'
#' @examples
netEmbedding <- function(object, slot.name = "netP", type = c("functional","structural"), pathway.remove = NULL, k = NULL) {
  Similarity <- methods::slot(object, slot.name)$similarity[[type]]$matrix
  if (is.null(pathway.remove)) {
    pathway.remove <- rownames(Similarity)[which(colSums(Similarity) == 1)]
  }
  if (length(pathway.remove) > 0) {
    pathway.remove.idx <- which(rownames(Similarity) %in% pathway.remove)
    Similarity <- Similarity[-pathway.remove.idx, -pathway.remove.idx]
  }
  if (is.null(k)) {
    k <- ceiling(sqrt(dim(Similarity)[1])) + 1
  }
  options(warn = -1)
  # dimension reduction
  Y <- runUMAP(Similarity, min.dist = 0.3, n.neighbors = k)
  methods::slot(object, slot.name)$similarity[[type]]$dr <- Y
  return(object)
}


#' Classification learning of the signaling networks
#'
#' @param object CellChat object
#' @param slot.name the slot name of object that is used to compute centrality measures of signaling networks
#' @param type "functional","structural"
#' @param k the number of signaling groups when running kmeans
#' @param methods the methods for clustering: "kmeans" or "spectral"
#' @param do.parallel whether doing parallel when inferring the number of signaling groups when running kmeans
#' @param nCores number of workers when doing parallel
#' @param k.eigen the number of eigenvalues used when doing spectral clustering
#' @importFrom methods slot
#' @importFrom future nbrOfWorkers plan
#' @importFrom future.apply future_sapply
#' @importFrom pbapply pbsapply
#' @return
#' @export
#'
#' @examples
netClustering <- function(object, slot.name = "netP", type = c("functional","structural"), k = NULL, methods = "kmeans", do.parallel = TRUE, nCores = 4, k.eigen = NULL) {
  type <- match.arg(type)
  Y <- methods::slot(object, slot.name)$similarity[[type]]$dr
  data.use <- Y
  if (methods == "kmeans") {
    if (!is.null(k)) {
      clusters = kmeans(data.use,k,nstart=10)$cluster
    } else {
      N <- nrow(data.use)
      kRange <- seq(2,min(N-1, 10),by = 1)
      if (do.parallel) {
        future::plan("multiprocess", workers = nCores)
        options(future.globals.maxSize = 1000 * 1024^2)
      }
      my.sapply <- ifelse(
        test = future::nbrOfWorkers() == 1,
        yes = pbapply::pbsapply,
        no = future.apply::future_sapply
      )
      results = my.sapply(
        X = 1:length(kRange),
        FUN = function(x) {
          idents <- kmeans(data.use,kRange[x],nstart=10)$cluster
          clusIndex <- idents
          #adjMat0 <- as.numeric(outer(clusIndex, clusIndex, FUN = "==")) - outer(1:N, 1:N, "==")
          adjMat0 <- Matrix::Matrix(as.numeric(outer(clusIndex, clusIndex, FUN = "==")), nrow = N, ncol = N)
          return(list(adjMat = adjMat0, ncluster = length(unique(idents))))
        },
        simplify = FALSE
      )
      adjMat <- lapply(results, "[[", 1)
      CM <- Reduce('+', adjMat)/length(kRange)

      numCluster <- computeEigengap(as.matrix(CM))$upper_bound
      clusters = kmeans(data.use,numCluster,nstart=10)$cluster
    }

  } else if (methods == "spectral") {
    A <- as.matrix(data.use)
    D <- apply(A, 1, sum)
    L <- diag(D)-A                       # unnormalized version
    L <- diag(D^-0.5)%*%L%*% diag(D^-0.5) # normalized version
    evL <- eigen(L,symmetric=TRUE)  # evL$values is decreasing sorted when symmetric=TRUE
    # pick the first k first k eigenvectors (corresponding k smallest) as data points in spectral space
    plot(rev(evL$values)[1:30])
    Z <- evL$vectors[,(ncol(evL$vectors)-k.eigen+1):ncol(evL$vectors)]
    clusters = kmeans(Z,k,nstart=20)$cluster
  }
  methods::slot(object, slot.name)$similarity[[type]]$group <- clusters
  return(object)
}


#' Build SNN matrix
# #' Adapted from swne (https://github.com/yanwu2014/swne)
#' @param data.use Features x samples matrix to use to build the SNN
#' @param k Defines k for the k-nearest neighbor algorithm
#' @param k.scale Granularity option for k.param
#' @param prune.SNN Sets the cutoff for acceptable Jaccard distances when
#'                  computing the neighborhood overlap for the SNN construction.
#'
#' @return Returns similarity matrix in sparse matrix format
#'
#' @importFrom FNN get.knn
#' @importFrom Matrix sparseMatrix
#' @export
#'
buildSNN <- function(data.use, k = 10, k.scale = 10, prune.SNN = 1/15) {
  n.cells <- ncol(data.use)
  if (n.cells < k) {
    stop("k cannot be greater than the number of samples")
  }

  ## find the k-nearest neighbors for each single cell
  my.knn <- FNN::get.knn(t(as.matrix(data.use)), k = min(k.scale * k, n.cells - 1))
  nn.ranked <- cbind(1:n.cells, my.knn$nn.index[, 1:(k - 1)])
  nn.large <- my.knn$nn.index

  w <- ComputeSNN(nn.ranked, prune.SNN)
  colnames(w) <- rownames(w) <- colnames(data.use)

  Matrix::diag(w) <- 1
  return(w)
}



#' Compute the eigengap of a given matrix for inferring the number of clusters
#'
#' @param CM consensus matrix
#' @param tau truncated consensus matrix
#' @param tol tolerance
#' @param do.plot whether showing the eigenspectrum; Default will save the plot
#' @return
#' @import ggplot2
#' @export
computeEigengap <- function(CM, tau = NULL, tol = 0.01, do.plot = TRUE){
  # compute the drop tolerance, enforcing parsimony of components
  K.init <- computeLaplacian(CM, tol = tol)$n_zeros
  if (is.null(tau)) {
    if (K.init <= 5) {
      tau = 0.3
    } else if (K.init <= 10){
      tau = 0.4
    } else {
      tau = 0.5
    }
  }

  # truncate the ensemble consensus matrix
  CM[CM <= tau] <- 0;
  # normalize and make symmetric
  CM <- (CM + t(CM))/2
  eigs <- computeLaplacian(CM, tol = tol)

  # compute the largest eigengap
  gaps <- diff(eigs$val)
  upper_bound <- which(gaps == max(gaps))

  # compute the number of zero eigenvalues
  lower_bound <- eigs$n_zeros

  if (do.plot) {
    df <- data.frame(nCluster = 1:min(c(30,length(eigs$val))), eigenVal = eigs$val[1:min(c(30,length(eigs$val)))])
    g <- ggplot(df, aes(x = nCluster, y = eigenVal)) + geom_point(size = 1) +
      geom_point(aes(x= upper_bound, y= eigs$val[upper_bound]), colour="red", size = 3, pch = 1) + theme(legend.position="none")
    title.name <- paste0('Inferred number of clusters: ', upper_bound,'; Min number: ', lower_bound)
    g <- g + labs(title = title.name) + theme_bw() + scale_x_continuous(breaks=seq(0,30,5)) +
      theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
      theme(text = element_text(size = 10)) + labs(x = 'Number of clusters', y = 'Eigenvalue of graph Laplacian')+
      theme(axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8))
    ggsave(filename= paste0("estimationNumCluster_eigenspectrum",sample.int(100,1),".pdf"), plot=g, width = 3.5, height = 3, units = 'in', dpi = 300)
  }
  return(list(upper_bound = upper_bound,
              lower_bound = lower_bound,
              eigs = eigs))
}


#' Compute eigenvalues of associated Laplacian matrix of a given matrix
#'
#' @param CM consensus matrix
#' @param tol tolerance
#' @return
#' @importFrom RSpectra eigs_sym
#' @importFrom Matrix colSums
#' @export
computeLaplacian <- function(CM, tol = 0.01) {
  # Normalized Laplacian:
  Dsq <- sqrt(Matrix::colSums(CM))
  L <- -Matrix::t(CM / Dsq) / Dsq
  Matrix::diag(L) <- 1 + Matrix::diag(L)

  numEigs <- min(100,nrow(CM))
  res <- RSpectra::eigs_sym(L, k = numEigs, which = "SM", opt = list(tol = 1e-4))
  eigs <- abs(Re(res$values))
  n_zeros <- sum(eigs <= tol)
  return(list(val = sort(eigs), n_zeros = n_zeros))
}


#' Rank the similarity of the shared signaling pathways based on their joint manifold learning
#'
#' @param object CellChat object
#' @param slot.name the slot name of object that is used to compute centrality measures of signaling networks
#' @param type "functional","structural"
#' @param comparison a numerical vector giving the datasets for comparison
#' @param pathway.remove a character vector defining the signaling to remove
#' @param x.rotation rotation of x-labels
#' @param title main title of the plot
#' @param bar.w the width of bar plot
#' @param color.use defining the color for each cell group
#' @param font.size font size
#' @import ggplot2
#' @importFrom methods slot
#' @return
#' @export
#'
#' @examples
rankSimilarity <- function(object, slot.name = "netP", type = c("functional","structural"), comparison = c(1,2),  pathway.remove = NULL,
                           x.rotation = 90, title = NULL, color.use = NULL, bar.w = NULL, font.size = 8) {
  type <- match.arg(type)
  net <- list()
  for (i in 1:length(setdiff(names(methods::slot(object, slot.name)), "similarity"))) {
    net[[i]] = methods::slot(object, slot.name)[[i]]$prob
  }
  net.dim <- sapply(net, dim)[3,]
  position <- cumsum(net.dim); position <- c(0,position)
  if (is.null(pathway.remove)) {
    similarity <- methods::slot(object, slot.name)$similarity[[type]]$matrix
    pathway.remove <- rownames(similarity)[which(colSums(similarity) == 1)]
    pathway.remove.idx <- which(rownames(similarity) %in% pathway.remove)
  }
  if (length(pathway.remove.idx) > 0) {
    for (i in 1:length(pathway.remove.idx)) {
      idx <- which(position - pathway.remove.idx[i] > 0)
      position[idx[1]] <- position[idx[1]] - 1
      if (idx[1] == 2) {
        position[3] <- position[3] - 1
      }
    }
  }

  Y <- methods::slot(object, slot.name)$similarity[[type]]$dr
  data1 <- Y[(position[comparison[1]]+1):position[comparison[1]+1], ]
  data2 <- Y[(position[comparison[2]]+1):position[comparison[2]+1], ]

  pathway.show = as.character(intersect(rownames(data1), rownames(data2)))
  data1 <- data1[pathway.show, ]
  data2 <- data2[pathway.show, ]
  euc.dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))
  dist <- NULL
  for(i in 1:nrow(data1)) dist[i] <- euc.dist(data1[i,],data2[i,])
  df <- data.frame(name = pathway.show, dist = dist, row.names = pathway.show)
  df <- df[order(df$dist), , drop = F]
  df$name <- factor(df$name, levels = as.character(df$name))

  gg <- ggplot(df, aes(x=name, y=dist)) + geom_bar(stat="identity",width = bar.w) +
    theme_classic() + theme(text=element_text(size=font.size),axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.y = element_text(size=font.size)) +
    xlab("") + ylab("Pathway distance") + coord_flip()#+
  if (!is.null(title)) {
    gg <- gg + ggtitle(title)+ theme(plot.title = element_text(hjust = 0.5))
  }
  if (!is.null(color.use)) {
    gg <- gg + scale_fill_manual(values = ggplot2::alpha(color.use, alpha = 1), drop = FALSE, na.value = "white")
    gg <- gg + scale_colour_manual(values = color.use, drop = FALSE, na.value = "white")
  }
  return(gg)
}



#' Rank signaling networks based on the information flow
#'
#' @param object CellChat object
#' @param slot.name the slot name of object that is used to compute centrality measures of signaling networks
#' @param mode "single","comparison"
#' @param comparison a numerical vector giving the datasets for comparison; a single value means ranking for only one dataset and two values means ranking comparison for two datasets
#' @param x.rotation rotation of x-labels
#' @param title main title of the plot
#' @param bar.w the width of bar plot
#' @param color.use defining the color for each cell group
#' @param font.size font size

#' @import ggplot2
#' @importFrom methods slot
#' @return
#' @export
#'
#' @examples
rankNet <- function(object, slot.name = "netP", mode = c("single","comparison"), comparison = c(1,2), x.rotation = 90, title = NULL, color.use = NULL, bar.w = 0.75, font.size = 8) {
  mode <- match.arg(mode)
  object.names <- names(methods::slot(object, slot.name))

  if (mode == "single") {
    object1 <- methods::slot(object, slot.name)[[comparison[1]]]
    prob1 = object1$prob
    pSum <- apply(prob1, 3, sum)
    pSum <- pSum
    pSum <- -1/log(pSum)
    pSum[is.infinite(pSum) | pSum < 0] <- max(pSum+1)
    pSum[is.na(pSum)] <- 0
    pair.name <- names(pSum)

    df<- data.frame(name = pair.name, contribution = pSum)
    idx <- with(df, order(df$contribution))
    df <- df[idx, ]
    df$name <- factor(df$name, levels = as.character(df$name))
    gg <- ggplot(df, aes(x=name, y=contribution)) + geom_bar(stat="identity",width = bar.w) +
      theme_classic() + theme(axis.text=element_text(size=10),axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.y = element_text(size=10)) +
      xlab("") + ylab("Information flow") + coord_flip()#+
    if (!is.null(title)) {
      gg <- gg + ggtitle(title)+ theme(plot.title = element_text(hjust = 0.5))
    }

    return(gg)

  } else {
    if (is.null(title)) {
      title <- object.names
    }
    object1 <- methods::slot(object, slot.name)[[comparison[1]]]
    object2 <- methods::slot(object, slot.name)[[comparison[2]]]
    prob1 = object1$prob
    prob2 = object2$prob

    pSum <- apply(prob1, 3, sum)
    pSum2 <- apply(prob2, 3, sum)

    pSum1 <- pSum
    pSum1 <- -1/log(pSum1)
    pSum1[is.infinite(pSum1) | pSum1 < 0] <- max(pSum1+1)

    pSum1[is.na(pSum1)] <- 0
    pair.name1 <- names(pSum1)
    pSum2 <- pSum2
    pSum2 <- -1/log(pSum2)
    pSum2[is.infinite(pSum2) | pSum2 < 0] <- max(pSum2+1)

    pSum2[is.na(pSum2)] <- 0
    pair.name2 <- names(pSum2)

    pair.name <- union(pair.name1, pair.name2)
    df1 <- data.frame(name = pair.name, contribution = 0, group = title[comparison[1]], row.names = pair.name)
    df1[pair.name1,2] <- -pSum1
    df2 <- data.frame(name = pair.name, contribution = 0, group = title[comparison[2]], row.names = pair.name)
    df2[pair.name2,2] <- pSum2

    contribution.relative <- as.numeric(format(df2$contribution/abs(df1$contribution), digits=1))
    df1$contribution.relative <- contribution.relative
    df2$contribution.relative <- contribution.relative
    idx <- with(df1, order(-contribution.relative, -contribution))
    df1 <- df1[idx, ]
    df2 <- df2[idx, ]
    df1$name <- factor(df1$name, levels = as.character(df1$name))
    df2$name <- factor(df2$name, levels = as.character(df2$name))
    df <- rbind(df1, df2)

    if (is.null(color.use)) {
      color.use =  ggPalette(2)
    }

    if (!is.null(color.use)) {
      colors.text <- ifelse(df$contribution.relative < 1, color.use[1], ifelse(df$contribution.relative > 1, color.use[2], "black"))
      gg <- ggplot(df, aes(x=name, y=contribution, fill = group)) + geom_bar(stat="identity",width = bar.w) +
        theme_classic() + theme(axis.text=element_text(size=font.size), axis.text.x = element_blank(),axis.ticks.x = element_blank(), axis.title.y = element_text(size=font.size)) +
        xlab("") + ylab("Information flow") + coord_flip()#+
      gg <- gg + scale_fill_manual(name = "", values = color.use) + theme(axis.text.y = element_text(colour = colors.text))
    } else {
      colors.text <- ifelse(df$contribution.relative < 1, "black", ifelse(df$contribution.relative > 1, "#696969", "#C0C0C0"))
      gg <- ggplot(df, aes(x=name, y=contribution)) + geom_bar(data = df[df$group == title[comparison[1]], ],fill = "grey50", color = "grey50",stat="identity",width = bar.w)
      gg <- gg + geom_bar(data = df[df$group == title[comparison[2]], ], fill = "white", color = "grey50",stat="identity",width = bar.w)
      gg <- gg + theme_classic() + theme(axis.text=element_text(size=font.size),axis.text.x = element_blank(),axis.ticks.x = element_blank(), axis.title.y = element_text(size=font.size)) +
        xlab("") + ylab("Information flow") + coord_flip() #+
      gg <- gg + theme(axis.text.y = element_text(colour = colors.text))
    }
    return(gg)
  }
}

#' Rank ligand-receptor interactions for any pair of two cell groups
#'
#' @param object CellChat object
#' @param LR.use ligand-receptor interactions used in inferring communication network
#' @return
#' @export
#'
rankNetPairwise <- function(object, LR.use = NULL) {
  if (is.null(LR.use)) {
    pairLR.use <- object@LR$LRsig
  } else {
    pairLR.use <- LR.use
  }
  net <- object@net
  prob <- net$prob
  pval <- net$pval
  numCluster <- dim(prob)[1]
  pairwiseLR <- list()
  for (i in 1:numCluster) {
    temp <- list()
    for (j in 1:numCluster) {
      pvalij <- pval[i,j,]; pvalij <- as.vector(pvalij)
      probij <- prob[i,j,]; probij <- as.vector(probij)
      index <- 1:length(pvalij)
      data <- data.frame(pathway_index = index, interaction_name = pairLR.use$interaction_name, interaction_name_2 = pairLR.use$interaction_name_2, pathway_name = pairLR.use$pathway_name, ligand = pairLR.use$ligand, receptor = pairLR.use$receptor,
                         prob = probij, pval = pvalij, row.names = rownames(pairLR.use))
      temp[[j]] <- data[with(data, order(pval, -prob)), ]
    }
    names(temp) <- colnames(prob)
    pairwiseLR[[i]] <- temp
  }
  names(pairwiseLR) <- rownames(prob)
  object@net$pairwiseRank <- pairwiseLR
  return(object)
}


#' compute the Shannon entropy
#'
#' @param a a numeric vector
#' @return
entropia<-function(a){
  a<-a[which(a>0)]
  return(-sum(a*log(a)))
}


#' compute the node distance matrix
#'
#' @param g a graph objecct
#' @return
node_distance<-function(g){
  n<-length(V(g))
  if(n==1){
    retorno=1
  }

  if(n>1){
    a<-Matrix::Matrix(0,nrow=n,ncol=n,sparse=TRUE)
    m<-shortest.paths(g,algorithm=c("unweighted"))
    m[which(m=="Inf")]<-n
    quem<-setdiff(intersect(m,m),0)
    for(j in (1:length(quem))){

      l<-which(m==quem[j])/n

      linhas<-floor(l)+1

      posicoesm1<-which(l==floor(l))

      if(length(posicoesm1)>0){
        linhas[posicoesm1]<-linhas[posicoesm1]-1
      }
      a[1:n,quem[j]]<-hist(linhas,plot=FALSE,breaks=(0:n))$counts

    }
    retorno=(a/(n-1))
  }
  return(retorno)
}


#' compute nnd
#'
#' @param g a graph objecct
#' @return
nnd<-function(g){

  N<-length(V(g))

  nd<-node_distance(g)

  pdfm<-Matrix::colMeans(nd)

  norm<-log(max(c(2,length(which(pdfm[1:(N-1)]>0))+1)))

  return(c(pdfm,max(c(0,entropia(pdfm)-entropia(as.matrix(nd))/N))/norm))
}

#' compute alpha centrality
#'
#' @param g a graph objecct
#' @return
alpha_centrality<-function(g){

  N<-length(V(g))

  r<-sort(igraph::alpha.centrality(g,exo=degree(g)/(N-1),alpha=1/N))/((N^2))

  return(c(r,max(c(0,1-sum(r)))))

}

#' Compute the structural distance between two signaling networks
#'
#' @param g a graph object of one signaling network
#' @param h a graph object of another signaling network
#' @param w1 parameter
#' @param w2 parameter
#' @param w3 parameter
#' @import igraph
#' @return
#' @export
#'
#' @examples
computeNetD_structure <- function(g, h, w1 = 0.45, w2 = 0.45, w3 = 0.1){

  first<-0

  second<-0

  third<-0

  # g<-read.graph(g,format=c("edgelist"),directed=FALSE)
  #
  # h<-read.graph(h,format=c("edgelist"),directed=FALSE)

  g <- graph_from_adjacency_matrix(g,mode="directed")
  h <- graph_from_adjacency_matrix(h,mode="directed")

  N<-length(V(g))

  M<-length(V(h))

  PM<-matrix(0,ncol=max(c(M,N)))

  if(w1+w2>0){

    pg = nnd(g)

    PM[1:(N-1)]=pg[1:(N-1)]

    PM[length(PM)]<-pg[N]

    ph=nnd(h)

    PM[1:(M-1)]=PM[1:(M-1)]+ph[1:(M-1)]

    PM[length(PM)]<-PM[length(PM)]+ph[M]

    PM<-PM/2

    first<-sqrt(max(c((entropia(PM)-(entropia(pg[1:N])+entropia(ph[1:M]))/2)/log(2),0)))

    second<-abs(sqrt(pg[N+1])-sqrt(ph[M+1]))


  }

  if(w3>0){

    pg<-alpha_centrality(g)

    ph<-alpha_centrality(h)

    m<-max(c(length(pg),length(ph)))

    Pg<-matrix(0,ncol=m)

    Ph<-matrix(0,ncol=m)

    Pg[(m-length(pg)+1):m]<-pg

    Ph[(m-length(ph)+1):m]<-ph

    third<-third+sqrt((entropia((Pg+Ph)/2)-(entropia(pg)+entropia(ph))/2)/log(2))/2

    g<-graph.complementer(g)

    h<-graph.complementer(h)


    pg<-alpha_centrality(g)

    ph<-alpha_centrality(h)

    m<-max(c(length(pg),length(ph)))

    Pg<-matrix(0,ncol=m)

    Ph<-matrix(0,ncol=m)

    Pg[(m-length(pg)+1):m]<-pg

    Ph[(m-length(ph)+1):m]<-ph

    third<-third+sqrt((entropia((Pg+Ph)/2)-(entropia(pg)+entropia(ph))/2)/log(2))/2
  }
  return(w1*first+w2*second+w3*third)
}
