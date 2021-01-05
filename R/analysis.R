
#' Compute and visualize the contribution of each ligand-receptor pair in the overall signaling pathways
#'
#' @param object CellChat object
#' @param signaling a signaling pathway name
#' @param signaling.name alternative signaling pathway name to show on the plot
#' @param width the width of individual bar
#' @param vertex.receiver a numeric vector giving the index of the cell groups as targets in the first hierarchy plot
#' @param thresh threshold of the p-value for determining significant interaction
#' @param return.data whether return the data.frame consisting of the predicted L-R pairs and their contribution
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
netAnalysis_contribution <- function(object, signaling, signaling.name = NULL, width = 0.1, vertex.receiver = NULL, thresh = 0.05, return.data = FALSE,
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
    if (!is.null(pairLR.name.use)) {
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
  if (return.data) {
    return(list(LR.contribution = df, gg.obj = gg))
  } else {
    return(gg)
  }
}


#' Compute the network centrality scores allowing identification of dominant senders, receivers, mediators and influencers in all inferred communication networks
#'
#' NB: This function was previously named as `netAnalysis_signalingRole`.  The previous function `netVisual_signalingRole` is now named as `netAnalysis_signalingRole_network`.
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
netAnalysis_computeCentrality <- function(object, slot.name = "netP", net = NULL, net.name = NULL) {
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
#' @importFrom igraph graph_from_adjacency_matrix strength hub_score authority_score eigen_centrality page_rank betweenness E
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
  igraph::E(G)$weight <- 1/igraph::E(G)$weight
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


#' Select the number of the patterns for running `identifyCommunicationPatterns`
#'
#' We infer the number of patterns based on two metrics that have been implemented in the NMF R package, including Cophenetic and Silhouette. Both metrics measure the stability for a particular number of patterns based on a hierarchical clustering of the consensus matrix. For a range of the number of patterns, a suitable number of patterns is the one at which Cophenetic and Silhouette values begin to drop suddenly.
#'
#' @param object CellChat object
#' @param slot.name the slot name of object that is used to compute centrality measures of signaling networks
#' @param pattern "outgoing" or "incoming"
#' @param k.range a range of the number of patterns
#' @param title.name title of plot
#' @param do.facet whether use facet plot showing the two measures
#' @param nrun number of runs when performing NMF
#' @param seed.use seed when performing NMF
#' @importFrom methods slot
# #' @importFrom NMF nmfEstimateRank
#' @import NMF
# #' @importFrom ggplot2 scale_color_brewer
#' @import ggplot2
#' @return a ggplot object
#' @export
#'
#' @examples
selectK <- function(object, slot.name = "netP", pattern = c("outgoing","incoming"), title.name = NULL, do.facet = TRUE, k.range = seq(2,10), nrun = 30, seed.use = 10) {
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

  if (is.null(title.name)) {
    title.name <- paste0(pattern, " signaling \n")
    # title.name <- paste0(pattern, " signaling \n (nrun = ", nrun, ", seed = ", seed.use, ")")
  }

  res <- NMF::nmfEstimateRank(data, range = k.range, method = 'lee', nrun=nrun, seed = seed.use)
  df1 <- data.frame(k = res$measures$rank, score = res$measures$cophenetic, Measure = "Cophenetic")
  df2 <- data.frame(k = res$measures$rank, score = res$measures$silhouette.consensus, Measure = "Silhouette")
  # df3 <- data.frame(k = res$measures$rank, score = res$measures$dispersion, Measure = "Dispersion")
  df <- rbind(df1, df2)
  #df <- rbind(df1, df2, df3)
  gg <- ggplot(df, aes(x = k, y = score, group = Measure, color = Measure)) + geom_line(size=1) +
    geom_point() +
    theme_classic() + labs(x = 'Number of patterns', y='Measure score') +
    labs(title = title.name) +  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
    theme(legend.position = "right") + theme(text = element_text(size = 10)) + scale_x_discrete(limits = (unique(df$k))) +
    scale_color_brewer(palette="Set2") + guides(color=guide_legend("Measure type"))
  if (do.facet) {
    gg <- gg + facet_wrap(~ Measure, scales='free')
  }
  gg
  return(gg)
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
    stop("Please run the function `selectK` for selecting a suitable k!")
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
    color.heatmap = grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n = 9, name = color.heatmap)))(255)

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
#' @param stacked whether plot the stacked bar plot
#' @param axis.gap whetehr making gaps in y-axes
#' @param ylim,segments,tick_width,rel_heights parameters in the function gg.gap when making gaps in y-axes
#' e.g., ylim = c(0, 35), segments = list(c(11, 14),c(16, 28)), tick_width = c(5,2,5), rel_heights = c(0.8,0,0.1,0,0.1)
#' https://tobiasbusch.xyz/an-r-package-for-everything-ep2-gaps
#' @param tol a tolerance when considering the relative contribution being equal between two datasets. contribution.relative between 1-tol and 1+tol will be considered as equal contribution
#' @param return.data whether return the data.frame consisting of the calculated information flow of each signaling pathway or L-R pair
#' @param show.raw whether show the raw information flow. Default = FALSE, showing the scaled information flow to provide compariable data scale
#' @param x.rotation rotation of x-labels
#' @param title main title of the plot
#' @param bar.w the width of bar plot
#' @param color.use defining the color for each cell group
#' @param font.size font size

#' @import ggplot2
#' @importFrom methods slot
#' @importFrom gg.gap gg.gap
#' @return
#' @export
#'
#' @examples
rankNet <- function(object, slot.name = "netP", mode = c("single","comparison"), comparison = c(1,2), stacked = FALSE, axis.gap = FALSE, ylim = NULL, segments = NULL, tick_width = NULL, rel_heights = c(0.9,0,0.1), tol = 0.05, return.data = FALSE, show.raw = FALSE, x.rotation = 90, title = NULL, color.use = NULL, bar.w = 0.75, font.size = 8) {
  mode <- match.arg(mode)
  options(warn = -1)
  object.names <- names(methods::slot(object, slot.name))
  if (mode == "single") {
    object1 <- methods::slot(object, slot.name)
    prob1 = object1$prob
    pSum <- apply(prob1, 3, sum)
    pSum.original <- pSum
    pSum <- -1/log(pSum)
    pSum[is.na(pSum)] <- 0
    idx1 <- which(is.infinite(pSum) | pSum < 0)
    values.assign <- seq(max(pSum)*1.1, max(pSum)*1.5, length.out = length(idx1))
    position <- sort(pSum.original[idx1], index.return = TRUE)$ix
    pSum[idx1] <- values.assign[match(1:length(idx1), position)]
    pair.name <- names(pSum)

    df<- data.frame(name = pair.name, contribution = pSum.original, contribution.scaled = pSum, group = object.names[comparison[1]])
    idx <- with(df, order(df$contribution.scaled))
    df <- df[idx, ]
    df$name <- factor(df$name, levels = as.character(df$name))
    gg <- ggplot(df, aes(x=name, y=contribution.scaled)) + geom_bar(stat="identity",width = bar.w) +
      theme_classic() + theme(axis.text=element_text(size=10),axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.y = element_text(size=10)) +
      xlab("") + ylab("Information flow") + coord_flip()#+
    if (!is.null(title)) {
      gg <- gg + ggtitle(title)+ theme(plot.title = element_text(hjust = 0.5))
    }

  } else {
    object1 <- methods::slot(object, slot.name)[[comparison[1]]]
    object2 <- methods::slot(object, slot.name)[[comparison[2]]]
    prob1 = object1$prob
    prob2 = object2$prob

    pSum <- apply(prob1, 3, sum)
    pSum2 <- apply(prob2, 3, sum)

    pSum1 <- pSum
    pSum1.original <- pSum1
    pSum1 <- -1/log(pSum1)
    pair.name1 <- names(pSum1)
    pSum2 <- pSum2
    pSum2.original <- pSum2
    pSum2 <- -1/log(pSum2)
    pair.name2 <- names(pSum2)

    pSum1[is.na(pSum1)] <- 0
    pSum2[is.na(pSum2)] <- 0
    idx1 <- which(is.infinite(pSum1) | pSum1 < 0)
    idx2 <- which(is.infinite(pSum2) | pSum2 < 0)
    values.assign <- seq(max(c(pSum1, pSum2))*1.1, max(c(pSum1, pSum2))*1.5, length.out = length(c(idx1, idx2)))
    position <- sort(c(pSum1.original[idx1], pSum2.original[idx2]), index.return = TRUE)$ix
    pSum1[idx1] <- values.assign[match(1:length(idx1), position)]
    pSum2[idx2] <- values.assign[match((length(idx1)+1):length(position), position)]

    pair.name <- union(pair.name1, pair.name2)
    df1 <- data.frame(name = pair.name, contribution = 0, contribution.scaled = 0, group = object.names[comparison[1]], row.names = pair.name)
    df1[pair.name1,3] <- pSum1
    df1[pair.name1,2] <- pSum1.original
    df2 <- data.frame(name = pair.name, contribution = 0, contribution.scaled = 0, group = object.names[comparison[2]], row.names = pair.name)
    df2[pair.name2,3] <- pSum2
    df2[pair.name2,2] <- pSum2.original

    contribution.relative <- as.numeric(format(df2$contribution/abs(df1$contribution), digits=1))
    df1$contribution.relative <- contribution.relative
    df2$contribution.relative <- contribution.relative
    df1$contribution.data2 <- df2$contribution
    #df1$contribution.data2[!is.infinite(df1$contribution.data2)] <- 0
    idx <- with(df1, order(-contribution.relative, contribution, -contribution.data2))
    df1 <- df1[idx, ]
    df2 <- df2[idx, ]
    df1$contribution.data2 <- NULL
    df1$name <- factor(df1$name, levels = as.character(df1$name))
    df2$name <- factor(df2$name, levels = as.character(df2$name))
    df <- rbind(df1, df2)
    df$group <- factor(df$group, levels = c(object.names[comparison[1]], object.names[comparison[2]]))

    if (is.null(color.use)) {
      color.use =  ggPalette(2)
    }
    if (!show.raw) {
      if (stacked) {
        colors.text <- ifelse(df$contribution.relative < 1-tol, color.use[1], ifelse(df$contribution.relative > 1+tol, color.use[2], "black"))
        gg <- ggplot(df, aes(x=name, y=contribution.scaled, fill = group)) + geom_bar(stat="identity",width = bar.w, position ="fill") +
          theme_classic() + theme(axis.text=element_text(size=font.size), axis.title.y = element_text(size=font.size)) +
          xlab("") + ylab("Relative information flow") + coord_flip()#+ theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())
        gg <- gg + scale_fill_manual(name = "", values = color.use) + theme(axis.text.y = element_text(colour = colors.text)) +
          # scale_y_continuous(labels = c(0,0.5,1)) +
          geom_hline(yintercept = 0.5, linetype="dashed", color = "grey50", size=0.5)
        if (!is.null(title)) {
          gg <- gg + ggtitle(title)+ theme(plot.title = element_text(hjust = 0.5))
        }
      } else {
        colors.text <- ifelse(df$contribution.relative < 1-tol, color.use[1], ifelse(df$contribution.relative > 1+tol, color.use[2], "black"))
        gg <- ggplot(df, aes(x=name, y=contribution.scaled, fill = group)) + geom_bar(stat="identity",width = bar.w, position = position_dodge(0.8)) +
          theme_classic() + theme(axis.text=element_text(size=font.size), axis.title.y = element_text(size=font.size)) +
          xlab("") + ylab("Information flow") + coord_flip()#+ theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())
        gg <- gg + scale_fill_manual(name = "", values = color.use) + theme(axis.text.y = element_text(colour = colors.text))
        if (axis.gap) {
          gg <- gg + theme_bw() + theme(panel.grid = element_blank())
          gg.gap::gg.gap(gg,
                         ylim = ylim,
                         segments = segments,
                         tick_width = tick_width,
                         rel_heights = rel_heights)
        }
        if (!is.null(title)) {
          gg <- gg + ggtitle(title)+ theme(plot.title = element_text(hjust = 0.5))
        }
      }

    } else {
      df$contribution[df$group == object.names[comparison[1]]] <- -df$contribution[df$group == object.names[comparison[1]]]
      colors.text <- ifelse(df$contribution.relative < 1, color.use[1], ifelse(df$contribution.relative > 1, color.use[2], "black"))
      gg <- ggplot(df, aes(x=name, y=contribution, fill = group)) + geom_bar(stat="identity",width = bar.w) +
        theme_classic() + theme(axis.text=element_text(size=font.size),axis.title.y = element_text(size=font.size)) +
        xlab("") + ylab("Information flow") + coord_flip()#+
      gg <- gg + scale_fill_manual(name = "", values = color.use) + theme(axis.text.y = element_text(colour = colors.text))
      if (!is.null(title)) {
        gg <- gg + ggtitle(title)+ theme(plot.title = element_text(hjust = 0.5))
      }
    }
  }

  if (return.data) {
    df$contribution <- abs(df$contribution)
    df$contribution.scaled <- abs(df$contribution.scaled)
    return(list(signaling.contribution = df, gg.obj = gg))
  } else {
    return(gg)
  }
}

#' Comparing the number of inferred communication links between different datasets
#'
#' @param object A merged CellChat object
#' @param measure "count" or "weight". "count": comparing the number of interactions; "weight": comparing the total interaction weights (strength)
#' @param group a vector giving the groups of different datasets to define colors of the bar plot. Default: only one group and a single color
#' @param group.levels the factor level in the defined group
#' @param color.use defining the color for each group of datasets
#' @param color.alpha transparency
#' @param legend.title legend title
#' @param width bar width
#' @param title.name main title of the plot
#' @param xlabel label of x-axis
#' @param ylabel label of y-axis
#' @param remove.xtick whether remove xtick
#' @param size.text font size of the text
#' @param show.legend whether show the legend
#' @param x.lab.rot,angle.x,vjust.x,hjust.x adjusting parameters if rotating xtick.labels when x.lab.rot = TRUE
#' @import ggplot2
#' @return A ggplot object
#' @export
#'
compareInteractions <- function(object, measure = c("count", "weight"), group = NULL, group.levels = NULL, color.use = NULL, color.alpha = 1, legend.title = NULL, width=0.6, title.name = NULL,
                                xlabel = NULL, ylabel = NULL, remove.xtick = FALSE,
                                show.legend = TRUE, x.lab.rot = FALSE, angle.x = 45, vjust.x = NULL, hjust.x = 1, size.text = 10) {
  measure <- match.arg(measure)
  if (measure == "count") {
    df <- as.data.frame(sapply(object@net, function(x) sum(x$count)))
    if (is.null(ylabel)) {
      ylabel = "Number of inferred interactions"
    }
  } else if (measure == "weight") {
    df <- as.data.frame(sapply(object@net, function(x) sum(x$weight)))
    df[,1] <- round(df[,1],1)
    if (is.null(ylabel)) {
      ylabel = "Interaction strength"
    }
  }
  colnames(df) <- "count"

  df$dataset <- names(object@net)
  if (is.null(group)) {
    group <- 1
  }
  df$group <- group
  df$dataset <- factor(df$dataset, levels = names(object@net))
  if (is.null(group.levels)) {
    df$group <- factor(df$group)
  } else {
    df$group <- factor(df$group, levels = group.levels)
  }

  if (is.null(color.use)) {
    color.use <- ggPalette(length(unique(group)))
  }
  gg <- ggplot(df, aes(x=dataset, y=count, fill = group)) + geom_bar(stat="identity", width=width, position=position_dodge())
  #   theme_classic() #+ scale_x_discrete(limits = (levels(df$x)))
  gg <- gg + geom_text(aes(label=count), vjust=-0.3, size=3, position = position_dodge(0.9))
  gg <- gg + ylab(ylabel) + xlab(xlabel) + theme_classic() +
    labs(title = title.name) +  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
    theme(text = element_text(size = size.text), axis.text = element_text(colour="black"))
  gg <- gg + scale_fill_manual(values = alpha(color.use, alpha = color.alpha), drop = FALSE)
  #  gg <- gg + scale_color_manual(values = alpha(color.use, alpha = 1), drop = FALSE) + guides(colour = FALSE)
  if (remove.xtick) {
    gg <- gg + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
  }
  if (is.null(legend.title)) {
    gg <- gg + theme(legend.title = element_blank())
  } else {
    gg <- gg + guides(fill=guide_legend(legend.title))
  }
  if (!show.legend) {
    gg <- gg + theme(legend.position = "none")
  }
  if (x.lab.rot) {
    gg <- gg + theme(axis.text.x = element_text(angle = angle.x, hjust = hjust.x, vjust = vjust.x, size=size.text))
  }
  gg
  return(gg)
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
#' @importFrom igraph degree alpha.centrality
#' @return
alpha_centrality<-function(g){

  N<-length(igraph::V(g))

  r<-sort(igraph::alpha.centrality(g,exo=igraph::degree(g)/(N-1),alpha=1/N))/((N^2))

  return(c(r,max(c(0,1-sum(r)))))

}

#' Compute the structural distance between two signaling networks
#'
#' @param g a graph object of one signaling network
#' @param h a graph object of another signaling network
#' @param w1 parameter
#' @param w2 parameter
#' @param w3 parameter
#' @importFrom igraph graph_from_adjacency_matrix V graph.complementer
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


#' Identify all the significant interactions (L-R pairs) and related signaling genes for a given signaling pathway
#'
#' @param object CellChat object
#' @param signaling a char vector containing signaling pathway names for searching
#' @param geneLR.return whether return the related signaling genes of enriched L-R pairs
#' @param enriched.only whether only return the identified enriched signaling genes in the database. Default = TRUE, returning the significantly enriched signaling interactions
#' @param thresh threshold of the p-value for determining significant interaction
#' @param geneInfo a dataframe with gene official symbol (there should be one column named `Symbol`)
#' @param complex_input signaling complex information from CellChatDB
#' @importFrom dplyr select
#'
#' @return The returned value depends on the input argument:
#'
#' When `geneLR.return = FALSE`, it returns a data frame containing the significant interactions (L-R pairs)
#'
#' When `geneLR.return = TRUE`, it returns a list, the first element is a data frame containing the significant interactions (L-R pairs), and the second is a vector containing the related signaling genes of enriched L-R pairs, which can be used for examining the gene expression pattern using the function \code{\link{plotGeneExpression}}
#'
#' @export
#'
extractEnrichedLR <- function(object, signaling, geneLR.return = FALSE, enriched.only = TRUE, thresh = 0.05, geneInfo = NULL, complex_input = NULL) {
  DB <- object@DB
  if (is.null(geneInfo)) {
    geneInfo = DB$geneInfo
  } else {
    DB$geneInfo = geneInfo
  }
  if (is.null(complex_input)) {
    complex_input = DB$complex
  } else {
    DB$complex = complex_input
  }
  pairLR.all <- c()
  geneLR.all <- c()
  net0 <- slot(object, "net")
  for (ii in 1:length(signaling)) {
    signaling.i <- signaling[ii]
    if (!is.list(net0[[1]])) {
      net <- net0
      LR <- object@LR
      res <- extractEnrichedLR_internal(net, LR, DB, signaling = signaling.i, enriched.only = enriched.only, thresh = thresh)
    } else {
      geneLR.t <- c()
      pairLR.t <- c()
      for (i in 1:length(net0)) {
        net <- net0[[i]]
        LR <- object@LR[[i]]
        res.t <- extractEnrichedLR_internal(net, LR, DB, signaling = signaling.i, enriched.only = enriched.only, thresh = thresh)
        geneLR.t <- BiocGenerics::union(geneLR.t, as.character(res.t[[1]]))
        pairLR.t <- BiocGenerics::union(pairLR.t, as.character(res.t[[2]]))
      }
      res <- list(geneLR.t, pairLR.t)
    }
    geneLR.all <- c(geneLR.all, as.character(res[[1]]))
    pairLR.all <- c(pairLR.all, as.character(res[[2]]))
  }
  pairLR.all <- data.frame(interaction_name = pairLR.all, stringsAsFactors = FALSE)

  if (geneLR.return) {
    return(list(pairLR = pairLR.all, geneLR = geneLR.all))
  } else {
    return(pairLR.all)
  }
}

#' Identify all the significant interactions (L-R pairs) and related signaling genes for a given signaling pathway
#'
#' @param net,LR,DB object@net object@LR object@DB
#' @param signaling a char vector containing signaling pathway names for searching
#' @param enriched.only whether only return the identified enriched signaling genes in the database. Default = TRUE, returning the significantly enriched signaling interactions
#' @param thresh threshold of the p-value for determining significant interaction
#' @importFrom dplyr select
#'
#' @return a list: list(geneLR, pairLR.name.use)
extractEnrichedLR_internal <- function(net, LR, DB, signaling, enriched.only = TRUE, thresh = 0.05){
  pairLR <- searchPair(signaling = signaling, pairLR.use = LR$LRsig, key = "pathway_name", matching.exact = T, pair.only = T)
  pairLR.name.use = dplyr::select(DB$interaction[rownames(pairLR),],"interaction_name")
  if (enriched.only) {
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
      message(paste0('There is no significant communication of ', signaling))
    } else {
      pairLR <- pairLR[pairLR.name.use,]
    }
  }
  geneL <- unique(pairLR$ligand)
  geneR <- unique(pairLR$receptor)
  geneL <- extractGeneSubset(geneL, DB$complex, DB$geneInfo)
  geneR <- extractGeneSubset(geneR, DB$complex, DB$geneInfo)
  geneLR <- c(geneL, geneR)
  return(list(geneLR, pairLR.name.use))
}



#' Filter cell-cell communication if there are only few number of cells in certain cell groups
#'
#' @param object CellChat object
#' @param min.cells the minmum number of cells required in each cell group for cell-cell communication
#' @return CellChat object with an updated slot net
#' @export
#'
filterCommunication <- function(object, min.cells = 10) {
  net <- object@net
  cell.excludes <- which(as.numeric(table(object@idents)) <= min.cells)
  if (length(cell.excludes) > 0) {
    cat("The cell-cell communication related with the following cell groups are excluded due to the few number of cells: ", levels(object@idents)[cell.excludes],'\n')
    net$prob[cell.excludes,,] <- 0
    net$prob[,cell.excludes,] <- 0
    object@net <- net
  }
  return(object)
}


#' Compute the maximum value of certain measures in the inferred cell-cell communication networks
#'
#' To better control the node size and edge weights of the inferred networks across different datasets,
#' we compute the maximum number of cells per cell group and the maximum number of interactions (or interaction weights) across all datasets
#'
#' @param object.list List of CellChat objects
#' @param slot.name the slot name of object that is used to compute the maximum value.
#'
#' When slot.name = "idents", 'attribute' should be "idents", which will compute the maximum number of cells per cell group across all datasets
#'
#' When slot.name = "net", 'attribute' can be either "count" or "weight", which will compute he maximum number of interactions (or interaction weights) across all datasets
#'
#' When slot.name = "net" or "netP", 'attribute' can be a single pathway name or a ligand-receptor pair name
#'
#' @param attribute the attribute to compute the maximum values. `attribute` should have the same length as `slot.name`.
#'
#' `attribute` can only be "count", "weight","count.merged","weight.merged" or a single pathway name or a ligand-receptor pair name
#'
#' @return A numeric vector
#' @export
#'
getMaxWeight <- function(object.list, slot.name = c("idents", "net"), attribute = c("idents", "count")) {
  weight <- c()
  for (i in 1:length(slot.name)) {
    if (slot.name[i] == "idents") {
      weight.all <- sapply(object.list, function (x) {max(as.numeric(table(slot(x, slot.name[i]))))})
    } else if ((slot.name[i] == "net") & (attribute[i] %in% c("count", "weight","count.merged","weight.merged"))) {
      weight.all <- sapply(object.list, function (x) {slot(x, slot.name[i])[[attribute[i]]]})
    } else if (attribute[i] %in% c(object.list[[1]]@DB$interaction$pathway_name, object.list[[1]]@DB$interaction$interaction_name)) {
      weight.all <- sapply(object.list, function (x) {max(slot(x, slot.name[i])$prob[,,attribute[i]])})
    }
    weight[i] <- max(weight.all)
  }
  names(weight) <- attribute
  weight.max <- weight
  return(weight.max)
}


#' Compute the number of interactions/interaction strength between cell types based on their associated cell subpopulations
#'
#' @param object CellChat object
#' @param group.merged a factor defining the group for merging different clusters/subpopulations
#'
#' @return An updated slot `net` by adding three elements:
#'
#' `count.merged`: the number of interactions between cell types (i.e., merged cell groups)
#'
#' `weight.merged`: interaction strength between cell types (i.e., merged cell groups)
#'
#' `group.merged` the defined group for merging different clusters/subpopulations
#'
#' @export
#'
mergeInteractions <- function(object, group.merged) {
  if (!is.factor(group.merged)) {
    group.merged <- factor(group.merged)
  }
  count <- object@net$count
  count.merged <- matrix(0, nrow = nlevels(group.merged), ncol = nlevels(group.merged))
  rownames(count.merged) <- levels(group.merged); colnames(count.merged) <- levels(group.merged);
  weight <- object@net$weight
  weight.merged <- count.merged
  dimnames(weight.merged) <- dimnames(count.merged)
  for (i in levels(group.merged)) {
    for (j in levels(group.merged)) {
      count.merged[i, j] <- sum(count[group.merged == i, group.merged == j])
      weight.merged[i, j] <- sum(weight[group.merged == i, group.merged == j])
    }
  }
  object@net$count.merged <- count.merged
  object@net$weight.merged <- weight.merged
  object@net$group.merged <- group.merged
  return(object)
}


#' Subset the inferred cell-cell communications of interest
#'
#' NB: If all arguments are NULL, it returns a data frame consisting of all the inferred cell-cell communications
#'
#' @param object CellChat object
#' @param slot.name the slot name of object: slot.name = "net" when extracting the inferred communications at the level of ligands/receptors; slot.name = "netP" when extracting the inferred communications at the level of signaling pathways
#' @param sources.use a vector giving the index or the name of source cell groups
#' @param targets.use a vector giving the index or the name of target cell groups.
#' @param signaling a character vector giving the name of signaling pathways of interest
#' @param pairLR.use a data frame consisting of one column named either "interaction_name" or "pathway_name", defining the interactions of interest
#' @param thresh threshold of the p-value for determining significant interaction
#' @importFrom  dplyr select group_by summarize groups
#' @importFrom stringr str_split
#' @importFrom BiocGenerics as.data.frame
#' @importFrom reshape2 melt
#' @importFrom magrittr %>%
#'
#' @return If input object is created from a single dataset, a data frame of the inferred cell-cell communications of interest, consisting of source, target, interaction_name, pathway_name, prob and other information
#'
#' If input object is a merged object from multiple datasets, it will return a list and each element is a data frame for one dataset
#'
#' @export
#'
#' @examples
#'\dontrun{
#' # access all the inferred cell-cell communications
#' df.net <- subsetCommunication(cellchat)
#'
#' # access all the inferred cell-cell communications at the level of signaling pathways
#' df.net <- subsetCommunication(cellchat, slot.name = "netP")
#'
#' # Subset to certain cells with sources.use and targets.use
#' df.net <- subsetCommunication(cellchat, sources.use = c(1,2), targets.use = c(4,5))
#'
#' # Subset to certain signaling, e.g., WNT and TGFb
#' df.net <- subsetCommunication(cellchat, signaling = c("WNT", "TGFb"))
#'}
#'
subsetCommunication <- function(object, slot.name = "net",
                                sources.use = NULL, targets.use = NULL,
                                signaling = NULL,
                                pairLR.use = NULL,
                                thresh = 0.05) {
  if (!is.null(pairLR.use)) {
    if (!is.data.frame(pairLR.use)) {
      stop("pairLR.use should be a data frame with a signle column named either 'interaction_name' or 'pathway_name' ")
    } else if ("pathway_name" %in% colnames(pairLR.use)) {
      message("slot.name is set to be 'netP' when pairLR.use contains signaling pathways")
      slot.name = "netP"
    }
  }

  if (!is.null(pairLR.use) & !is.null(signaling)) {
    stop("Please do not assign values to 'signaling' when using 'pairLR.use'")
  }

  net0 <- slot(object, "net")
  if (!is.list(object@net[[1]])) {
    net <- net0
    LR <- object@LR
    cells.level <- levels(object@idents)
    df.net <- subsetCommunication_internal(net, LR, cells.level, slot.name = slot.name,
                                           sources.use = sources.use, targets.use = targets.use,
                                           signaling = signaling,
                                           pairLR.use = pairLR.use,
                                           thresh = thresh)
  } else {
    df.net <- vector("list", length(net0))
    names(df.net) <- names(net0)
    for (i in 1:length(net0)) {
      net <- net0[[i]]
      LR <- object@LR[[i]]
      cells.level <- levels(object@idents[[i]])
      df.net[[i]] <- subsetCommunication_internal(net, LR, cells.level, slot.name = slot.name,
                                                  sources.use = sources.use, targets.use = targets.use,
                                                  signaling = signaling,
                                                  pairLR.use = pairLR.use,
                                                  thresh = thresh)
    }
  }

  return(df.net)

}



#' Subset the inferred cell-cell communications of interest
#'
#' NB: If all arguments are NULL, it returns a data frame consisting of all the inferred cell-cell communications
#'
#' @param net,LR,cells.level object@net object@LR levels(object@idents)
#' @param slot.name the slot name of object: slot.name = "net" when extracting the inferred communications at the level of ligands/receptors; slot.name = "netP" when extracting the inferred communications at the level of signaling pathways
#' @param sources.use a vector giving the index or the name of source cell groups
#' @param targets.use a vector giving the index or the name of target cell groups.
#' @param signaling a character vector giving the name of signaling pathways of interest
#' @param pairLR.use a data frame consisting of one column named either "interaction_name" or "pathway_name", defining the interactions of interest
#' @param thresh threshold of the p-value for determining significant interaction
#' @importFrom  dplyr select group_by summarize groups
#' @importFrom stringr str_split
#' @importFrom BiocGenerics as.data.frame
#' @importFrom reshape2 melt
#' @importFrom magrittr %>%
#'
#' @return A data frame of the inferred cell-cell communications of interest, consisting of source, target, interaction_name, pathway_name, prob and other information

subsetCommunication_internal <- function(net, LR, cells.level, slot.name = "net",
                                         sources.use = NULL, targets.use = NULL,
                                         signaling = NULL,
                                         pairLR.use = NULL,
                                         thresh = 0.05) {
  prob <- net$prob
  pval <- net$pval
  prob[pval > thresh] <- 0
  net <- reshape2::melt(prob, value.name = "prob")
  colnames(net)[1:3] <- c("source","target","interaction_name")
  net.pval <- reshape2::melt(pval, value.name = "pval")
  net$pval <- net.pval$pval
  # remove the interactions with zero values
  net <- subset(net, prob > 0)

  pairLR <- dplyr::select(LR$LRsig, c("interaction_name_2", "pathway_name", "ligand",  "receptor" ,"annotation","evidence"))
  idx <- match(net$interaction_name, rownames(pairLR))
  net <- cbind(net, pairLR[idx,])

  if (!is.null(signaling)) {
    pairLR.use <- data.frame()
    for (i in 1:length(signaling)) {
      pairLR.use.i <- searchPair(signaling = signaling[i], pairLR.use = LR$LRsig, key = "pathway_name", matching.exact = T, pair.only = T)
      pairLR.use <- rbind(pairLR.use, pairLR.use.i)
    }
  }

  if (!is.null(pairLR.use)){
    net <- tryCatch({
      subset(net,interaction_name %in% pairLR.use$interaction_name)
    }, error = function(e) {
      subset(net, pathway_name %in% pairLR.use$pathway_name)
    })
  }

  if (slot.name == "netP") {
    net <- dplyr::select(net, c("source","target","pathway_name","prob", "pval","annotation"))
    net$source_target <- paste(net$source, net$target, sep = "_")
    # net$source_target_pathway <- paste(paste(net$source, net$target, sep = "_"), net$pathway_name, sep = "_")
    net <- net %>% group_by(source_target, pathway_name) %>% summarize(prob = sum(prob), .groups = 'drop')
    net.pval <- net %>% group_by(source_target, pathway_name) %>% summarize(pval = mean(pval), .groups = 'drop')
    a <- stringr::str_split(net$source_target, "_", simplify = T)
    net$source <- as.character(a[, 1])
    net$target <- as.character(a[, 2])
    net <- dplyr::select(net, -source_target)
    net$pval <- net.pval$pval
  }

  # keep the interactions associated with sources and targets of interest
  if (!is.null(sources.use)){
    if (is.numeric(sources.use)) {
      sources.use <- cells.level[sources.use]
    }
    net <- subset(net, source %in% sources.use)
  }
  if (!is.null(targets.use)){
    if (is.numeric(targets.use)) {
      targets.use <- cells.level[targets.use]
    }
    net <- subset(net, target %in% targets.use)
  }

  net <- BiocGenerics::as.data.frame(net, stringsAsFactors=FALSE)

  if (nrow(net) == 0) {
    warning("No significant signaling interactions are inferred!")
  } else {
    rownames(net) <- 1:nrow(net)
  }

  if (slot.name == "net") {
    net <- net[,c("source", "target", "ligand", "receptor",  "prob", "pval", "interaction_name", "interaction_name_2", "pathway_name","annotation","evidence")]
  } else if (slot.name == "netP") {
    net <- net[,c("source", "target", "pathway_name", "prob", "pval")]
  }

  return(net)

}




#' Heatmap showing the centrality scores/importance of cell groups as senders, receivers, mediators and influencers in a single intercellular communication network
#'
#' @param object CellChat object
#' @param signaling a character vector giving the name of signaling networks
#' @param slot.name the slot name of object that is used to compute centrality measures of signaling networks
#' @param measure centrality measures to show
#' @param measure.name the names of centrality measures to show
#' @param color.use the character vector defining the color of each cell group
#' @param color.heatmap a color name in brewer.pal
#' @param width width of heatmap
#' @param height height of heatmap
#' @param font.size fontsize in heatmap
#' @param font.size.title font size of the title
#' @param cluster.rows whether cluster rows
#' @param cluster.cols whether cluster columns
#' @importFrom methods slot
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#' @importFrom ComplexHeatmap Heatmap HeatmapAnnotation draw
#' @importFrom stats setNames
#'
#' @return
#' @export
#'
#' @examples
netAnalysis_signalingRole_network <- function(object, signaling, slot.name = "netP", measure = c("outdeg","indeg","flowbet","info"), measure.name = c("Sender","Receiver","Mediator","Influencer"),
                                    color.use = NULL, color.heatmap = "BuGn",
                                    width = 6.5, height = 1.4, font.size = 8, font.size.title = 10, cluster.rows = FALSE, cluster.cols = FALSE) {
  if (length(slot(object, slot.name)$centr) == 0) {
    stop("Please run `netAnalysis_computeCentrality` to compute the network centrality scores! ")
  }
  centr <- slot(object, slot.name)$centr[signaling]
  for(i in 1:length(centr)) {
    centr0 <- centr[[i]]
    mat <- matrix(unlist(centr0), ncol = length(centr0), byrow = FALSE)
    mat <- t(mat)
    rownames(mat) <- names(centr0); colnames(mat) <- names(centr0$outdeg)
    if (!is.null(measure)) {
      mat <- mat[measure,]
      if (!is.null(measure.name)) {
        rownames(mat) <- measure.name
      }
    }
    mat <- sweep(mat, 1L, apply(mat, 1, max), '/', check.margin = FALSE)

    if (is.null(color.use)) {
      color.use <- scPalette(length(colnames(mat)))
    }
    color.heatmap.use = grDevices::colorRampPalette((RColorBrewer::brewer.pal(n = 9, name = color.heatmap)))(100)

    df<- data.frame(group = colnames(mat)); rownames(df) <- colnames(mat)
    cell.cols.assigned <- setNames(color.use, unique(as.character(df$group)))
    col_annotation <- HeatmapAnnotation(df = df, col = list(group = cell.cols.assigned),which = "column",
                                        show_legend = FALSE, show_annotation_name = FALSE,
                                        simple_anno_size = grid::unit(0.2, "cm"))

    ht1 = Heatmap(mat, col = color.heatmap.use, na_col = "white", name = "Importance",
                  bottom_annotation = col_annotation,
                  cluster_rows = cluster.rows,cluster_columns = cluster.rows,
                  row_names_side = "left",row_names_rot = 0,row_names_gp = gpar(fontsize = font.size),column_names_gp = gpar(fontsize = font.size),
                  width = unit(width, "cm"), height = unit(height, "cm"),
                  column_title = paste0(names(centr[i]), " signaling pathway network"),column_title_gp = gpar(fontsize = font.size.title),column_names_rot = 45,
                  heatmap_legend_param = list(title = "Importance", title_gp = gpar(fontsize = 8, fontface = "plain"),title_position = "leftcenter-rot",
                                              border = NA, at = c(round(min(mat, na.rm = T), digits = 1), round(max(mat, na.rm = T), digits = 1)),
                                              legend_height = unit(20, "mm"),labels_gp = gpar(fontsize = 8),grid_width = unit(2, "mm"))
    )
    draw(ht1)
  }
}


#' 2D visualization of dominant senders (sources) and receivers (targets)
#'
#' @description
#' This scatter plot shows the dominant senders (sources) and receivers (targets) in a 2D space. Dot size is proportional to the number of inferred links (both outgoing and incoming) associated with each cell group.
#' Dot colors indicate different cell groups. Dot shapes indicate different categories of cell groups if `group`` is defined.
#'
#' @param object CellChat object
#' @param signaling a char vector containing signaling pathway names. signaling = NULL: Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
#' @param color.use defining the color for each cell group
#' @param slot.name the slot name of object that is used to compute centrality measures of signaling networks
#' @param group a vector to categorize the cell groups, e.g., categorize the cell groups into two major categories: immune cells and fibroblasts
#' @param weight.MinMax the Minmum/maximum weight, which is useful to control the dot size when comparing multiple datasets
#' @param point.shape point shape when group is not NULL
#' @param label.size font size of the text
#' @param dot.alpha transparency
#' @param dot.size a range defining the size of the symbol
#' @param xlabel label of x-axis
#' @param ylabel label of y-axis
#' @param title main title of the plot
#' @param font.size font size of the text
#' @param font.size.title font size of the title
#' @param do.label label the each point
#' @param show.legend whether show the legend
#' @param show.axes whether show the axes
#' @import ggplot2
#' @importFrom ggrepel geom_text_repel
#' @importFrom methods slot
#' @return ggplot object
#' @export
#'
netAnalysis_signalingRole_scatter <- function(object, signaling = NULL, color.use = NULL, slot.name = "netP", group = NULL, weight.MinMax = NULL, dot.size = c(2, 6), point.shape = c(21, 22, 24, 23, 25, 8, 3), label.size = 3, dot.alpha = 0.6,
                                        xlabel = "Outgoing interaction strength", ylabel = "Incoming interaction strength", title = NULL,
                                        font.size = 10, font.size.title = 10, do.label = T, show.legend = T, show.axes = T) {
  if (length(slot(object, slot.name)$centr) == 0) {
    stop("Please run `netAnalysis_computeCentrality` to compute the network centrality scores! ")
  }
  centr <- slot(object, slot.name)$centr
  outgoing <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
  incoming <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
  dimnames(outgoing) <- list(levels(object@idents), names(centr))
  dimnames(incoming) <- dimnames(outgoing)
  for (i in 1:length(centr)) {
    outgoing[,i] <- centr[[i]]$outdeg
    incoming[,i] <- centr[[i]]$indeg
  }
  if (is.null(signaling)) {
    message("Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways")
  } else {
    message("Signaling role analysis on the cell-cell communication network from user's input")
    signaling <- signaling[signaling %in% object@netP$pathways]
    if (length(signaling) == 0) {
      stop('There is no significant communication for the input signaling. All the significant signaling are shown in `object@netP$pathways`')
    }
    outgoing <- outgoing[ , signaling, drop = FALSE]
    incoming <- incoming[ , signaling, drop = FALSE]
  }
  outgoing.cells <- rowSums(outgoing)
  incoming.cells <- rowSums(incoming)

  num.link <- aggregateNet(object, signaling = signaling, return.object = FALSE, remove.isolate = FALSE)$count
  num.link <- rowSums(num.link) + colSums(num.link)-diag(num.link)
  df <- data.frame(x = outgoing.cells, y = incoming.cells, labels = names(incoming.cells),
                   Count = num.link)
  if (!is.null(group)) {
    df$Group <- group
  }
  if (is.null(color.use)) {
    color.use <- scPalette(nlevels(object@idents))
  }
  if (!is.null(group)) {
    gg <- ggplot(data = df, aes(x, y)) +
      geom_point(aes(size = Count, colour = labels, fill = labels, shape = Group))
  } else {
    gg <- ggplot(data = df, aes(x, y)) +
      geom_point(aes(size = Count, colour = labels, fill = labels))
  }

  gg <- gg + CellChat_theme_opts() +
    theme(text = element_text(size = font.size), legend.key.height = grid::unit(0.15, "in"))+
    # guides(colour = guide_legend(override.aes = list(size = 3)))+
    labs(title = title, x = xlabel, y = ylabel) + theme(plot.title = element_text(size= font.size.title, face="plain"))+
    # theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank()) +
    theme(axis.line.x = element_line(size = 0.25), axis.line.y = element_line(size = 0.25))
  gg <- gg + scale_fill_manual(values = ggplot2::alpha(color.use, alpha = dot.alpha), drop = FALSE) + guides(fill=FALSE)
  gg <- gg + scale_colour_manual(values = color.use, drop = FALSE) + guides(colour=FALSE)
  # gg <- gg + scale_colour_manual(values = ggplot2::alpha(color.use, alpha = dot.alpha), drop = FALSE) + guides(colour=FALSE)
  # gg <- gg + scale_shape_manual(values = point.shape[1:length(prob)])
  if (!is.null(group)) {
    gg <- gg + scale_shape_manual(values = point.shape[1:length(unique(df$Group))])
  }
  if (is.null(weight.MinMax)) {
    gg <- gg + scale_size_continuous(range = dot.size)
  } else {
    gg <- gg + scale_size_continuous(limits = weight.MinMax, range = dot.size)
  }
  if (do.label) {
    gg <- gg + ggrepel::geom_text_repel(mapping = aes(label = labels, colour = labels), size = label.size, show.legend = F,segment.size = 0.2, segment.alpha = 0.5)
  }

  if (!show.legend) {
    gg <- gg + theme(legend.position = "none")
  }

  if (!show.axes) {
    gg <- gg + theme_void()
  }

  gg

}


#' Heatmap showing the contribution of signals (signaling pathways or ligand-receptor pairs) to cell groups in terms of outgoing or incoming signaling
#'
#' In this heatmap, colobar represents the relative signaling strength of a signaling pathway across cell groups (NB: values are row-scaled).
#' The top colored bar plot shows the total signaling strength of a cell group by summarizing all signaling pathways displayed in the heatmap.
#' The right grey bar plot shows the total signaling strength of a signaling pathway by summarizing all cell groups displayed in the heatmap.
#'
#' @param object CellChat object
#' @param signaling a character vector giving the name of signaling networks
#' @param pattern "outgoing", "incoming" or "all". When pattern = "all", it aggregates the outgoing and incoming signaling strength together
#' @param slot.name the slot name of object that is used to compute centrality measures of signaling networks
#' @param color.use the character vector defining the color of each cell group
#' @param color.heatmap a color name in brewer.pal
#' @param title title name
#' @param width width of heatmap
#' @param height height of heatmap
#' @param font.size fontsize in heatmap
#' @param font.size.title font size of the title
#' @param cluster.rows whether cluster rows
#' @param cluster.cols whether cluster columns
#' @importFrom methods slot
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#' @importFrom ComplexHeatmap Heatmap HeatmapAnnotation anno_barplot rowAnnotation
#' @importFrom stats setNames
#'
#' @return
#' @export
#'
netAnalysis_signalingRole_heatmap <- function(object, signaling = NULL, pattern = c("outgoing", "incoming","all"), slot.name = "netP",
                                              color.use = NULL, color.heatmap = "BuGn",
                                              title = NULL, width = 10, height = 8, font.size = 8, font.size.title = 10, cluster.rows = FALSE, cluster.cols = FALSE){
  pattern <- match.arg(pattern)
  if (length(slot(object, slot.name)$centr) == 0) {
    stop("Please run `netAnalysis_computeCentrality` to compute the network centrality scores! ")
  }
  centr <- slot(object, slot.name)$centr
  outgoing <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
  incoming <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
  dimnames(outgoing) <- list(levels(object@idents), names(centr))
  dimnames(incoming) <- dimnames(outgoing)
  for (i in 1:length(centr)) {
    outgoing[,i] <- centr[[i]]$outdeg
    incoming[,i] <- centr[[i]]$indeg
  }
  if (pattern == "outgoing") {
    mat <- t(outgoing)
    legend.name <- "Outgoing"
  } else if (pattern == "incoming") {
    mat <- t(incoming)
    legend.name <- "Incoming"
  } else if (pattern == "all") {
    mat <- t(outgoing+ incoming)
    legend.name <- "Overall"
  }
  if (is.null(title)) {
    title <- paste0(legend.name, " signaling patterns")
  } else {
    title <- paste0(paste0(legend.name, " signaling patterns"), " - ",title)
  }

  if (!is.null(signaling)) {
    mat1 <- mat[rownames(mat) %in% signaling, , drop = FALSE]
    mat <- matrix(0, nrow = length(signaling), ncol = ncol(mat))
    idx <- match(rownames(mat1), signaling)
    mat[idx[!is.na(idx)], ] <- mat1
    dimnames(mat) <- list(signaling, colnames(mat1))
  }
  mat.ori <- mat
  mat <- sweep(mat, 1L, apply(mat, 1, max), '/', check.margin = FALSE)
  mat[mat == 0] <- NA


  if (is.null(color.use)) {
    color.use <- scPalette(length(colnames(mat)))
  }
  color.heatmap.use = grDevices::colorRampPalette((RColorBrewer::brewer.pal(n = 9, name = color.heatmap)))(100)

  df<- data.frame(group = colnames(mat)); rownames(df) <- colnames(mat)
  names(color.use) <- colnames(mat)
  col_annotation <- HeatmapAnnotation(df = df, col = list(group = color.use),which = "column",
                                      show_legend = FALSE, show_annotation_name = FALSE,
                                      simple_anno_size = grid::unit(0.2, "cm"))
  ha2 = HeatmapAnnotation(Strength = anno_barplot(colSums(mat.ori), border = FALSE,gp = gpar(fill = color.use, col=color.use)), show_annotation_name = FALSE)

  pSum <- rowSums(mat.ori)
  pSum.original <- pSum
  pSum <- -1/log(pSum)
  pSum[is.na(pSum)] <- 0
  idx1 <- which(is.infinite(pSum) | pSum < 0)
  values.assign <- seq(max(pSum)*1.1, max(pSum)*1.5, length.out = length(idx1))
  position <- sort(pSum.original[idx1], index.return = TRUE)$ix
  pSum[idx1] <- values.assign[match(1:length(idx1), position)]
  ha1 = rowAnnotation(Strength = anno_barplot(pSum, border = FALSE), show_annotation_name = FALSE)

  ht1 = Heatmap(mat, col = color.heatmap.use, na_col = "white", name = "Relative strength",
                bottom_annotation = col_annotation, top_annotation = ha2, right_annotation = ha1,
                cluster_rows = cluster.rows,cluster_columns = cluster.rows,
                row_names_side = "left",row_names_rot = 0,row_names_gp = gpar(fontsize = font.size),column_names_gp = gpar(fontsize = font.size),
                width = unit(width, "cm"), height = unit(height, "cm"),
                column_title = title,column_title_gp = gpar(fontsize = font.size.title),column_names_rot = 90,
                heatmap_legend_param = list(title = "Relative strength", title_gp = gpar(fontsize = 8, fontface = "plain"),title_position = "leftcenter-rot",
                                            border = NA, at = c(round(min(mat, na.rm = T), digits = 1), round(max(mat, na.rm = T), digits = 1)),
                                            legend_height = unit(20, "mm"),labels_gp = gpar(fontsize = 8),grid_width = unit(2, "mm"))
  )
  #  draw(ht1)
  return(ht1)
}





