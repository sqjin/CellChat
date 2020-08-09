#' ggplot theme in CellChat
#'
#' @return
#' @export
#'
#' @examples
#' @importFrom ggplot2 theme_classic element_rect theme element_blank element_line element_text
CellChat_theme_opts <- function() {
  theme(strip.background = element_rect(colour = "white", fill = "white")) +
    theme_classic() +
    theme(panel.border = element_blank()) +
    theme(axis.line.x = element_line(color = "black")) +
    theme(axis.line.y = element_line(color = "black")) +
    theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) +
    theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) +
    theme(panel.background = element_rect(fill = "white")) +
    theme(legend.key = element_blank()) + theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5))
}


#' Generate ggplot2 colors
#'
#' @param n number of colors to generate
#' @importFrom grDevices hcl
#' @export
#'
ggPalette <- function(n) {
  hues = seq(15, 375, length = n + 1)
  grDevices::hcl(h = hues, l = 65, c = 100)[1:n]
}

#' Generate colors from a customed color palette
#'
#' @param n number of colors
#'
#' @return A color palette for plotting
#' @importFrom grDevices colorRampPalette
#'
#' @export
#'
scPalette <- function(n) {
  colorSpace <- c('#E41A1C','#377EB8','#4DAF4A','#984EA3','#F29403','#F781BF','#BC9DCC','#A65628','#54B0E4','#222F75','#1B9E77','#B2DF8A',
                  '#E3BE00','#FB9A99','#E7298A','#910241','#00CDD1','#A6CEE3','#CE1261','#5E4FA2','#8CA77B','#00441B','#DEDC00','#B3DE69','#8DD3C7','#999999')
  if (n <= length(colorSpace)) {
    colors <- colorSpace[1:n]
  } else {
    colors <- grDevices::colorRampPalette(colorSpace)(n)
  }
  return(colors)
}

#' Visualize the inferred cell-cell communication network
#'
#' Automatically save plots in the current working directory
#'
#' @param object CellChat object
#' @param signaling a signaling pathway name
#' @param signaling.name alternative signaling pathway name to show on the plot
#' @param vertex.receiver a numeric vector giving the index of the cell groups as targets in the first hierarchy plot
#' @param color.use the character vector defining the color of each cell group
#' @param vertex.size The size of vertex
#' @param layout "hierarchy" or "circle"
#' @param height height of plot
#' @param thresh threshold of the p-value for determining significant interaction
#' @param pt.title font size of the text
#' @param title.space the space between the title and plot
#' @param vertex.label.cex The label size of vertex in the network
#' @importFrom svglite svglite
#' @importFrom grDevices dev.off pdf
#'
#' @return
#' @export
#'
#' @examples
#'
netVisual <- function(object, signaling, signaling.name = NULL, vertex.receiver = NULL, color.use = NULL, vertex.size = 20, layout = c("hierarchy","circle"), height = 5, thresh = 0.05, pt.title = 12, title.space = 6, vertex.label.cex = 0.8) {
  layout <- match.arg(layout)
  pairLR <- searchPair(signaling = signaling, pairLR.use = object@LR$LRsig, key = "pathway_name", matching.exact = T, pair.only = F)

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
  nRow <- length(pairLR.name.use)

  prob <- prob[,,pairLR.name.use]
  pval <- pval[,,pairLR.name.use]

  if (length(dim(prob)) == 2) {
    prob <- replicate(1, prob, simplify="array")
    pval <- replicate(1, pval, simplify="array")
  }
  prob <-(prob-min(prob))/(max(prob)-min(prob))

  if (layout == "hierarchy") {
    svglite(file = paste0(signaling.name, "_hierarchy_individual.svg"), width = 8, height = nRow*height)
    par(mfrow=c(nRow,2), mar = c(5, 4, 4, 2) +0.1)
    for (i in 1:length(pairLR.name.use)) {
      #signalName_i <- paste0(pairLR$ligand[i], "-",pairLR$receptor[i], sep = "")
      signalName_i <- pairLR$interaction_name_2[i]
      prob.i <- prob[,,i]
      netVisual_hierarchy1(prob.i, vertex.receiver = vertex.receiver, color.use = color.use, vertex.size = vertex.size, signaling.name = signalName_i, vertex.label.cex = vertex.label.cex)
      netVisual_hierarchy2(prob.i, vertex.receiver = setdiff(1:nrow(prob.i),vertex.receiver), color.use = color.use, vertex.size = vertex.size, signaling.name = signalName_i, vertex.label.cex = vertex.label.cex)
    }
    dev.off()

    prob.sum <- apply(prob, c(1,2), sum)
    prob.sum <-(prob.sum-min(prob.sum))/(max(prob.sum)-min(prob.sum))
    svglite(file = paste0(signaling.name, "_hierarchy_aggregate.svg"), width = 7, height = 1*height)
    par(mfrow=c(1,2), ps = pt.title)
    netVisual_hierarchy1(prob.sum, vertex.receiver = vertex.receiver, color.use = color.use, vertex.size = vertex.size, signaling.name = NULL, vertex.label.cex = vertex.label.cex)
    netVisual_hierarchy2(prob.sum, vertex.receiver = setdiff(1:nrow(prob.sum),vertex.receiver), color.use = color.use, vertex.size = vertex.size, signaling.name = NULL, vertex.label.cex = vertex.label.cex)
    graphics::mtext(paste0(signaling.name, " signaling pathway network"), side = 3, outer = TRUE, cex = 1, line = -title.space)
    dev.off()
  } else if (layout == "circle") {
    svglite(file = paste0(signaling.name,"_", layout, "_individual.svg"), width = height, height = nRow*height)
    par(mfrow=c(nRow,1))
    for (i in 1:length(pairLR.name.use)) {
      #signalName_i <- paste0(pairLR$ligand[i], "-",pairLR$receptor[i], sep = "")
      signalName_i <- pairLR$interaction_name_2[i]
      prob.i <- prob[,,i]
      netVisual_circle(prob.i, top = 1, color.use = color.use, vertex.size = vertex.size, signaling.name = signalName_i, vertex.label.cex = vertex.label.cex)
    }
    dev.off()

    prob.sum <- apply(prob, c(1,2), sum)
    prob.sum <-(prob.sum-min(prob.sum))/(max(prob.sum)-min(prob.sum))
    svglite(file = paste0(signaling.name,"_", layout,  "_aggregate.svg"), width = height, height = 1*height)
    netVisual_circle(prob.sum, top = 1, color.use = color.use, vertex.size = vertex.size, signaling.name = paste0(signaling.name, " signaling pathway network"), vertex.label.cex = vertex.label.cex)
    dev.off()
  }


}


#' Visualize the inferred signaling network of signaling pathways by aggregating all L-R pairs
#'
#' @param object CellChat object
#' @param signaling a signaling pathway name
#' @param signaling.name alternative signaling pathway name to show on the plot
#' @param vertex.receiver a numeric vector giving the index of the cell groups as targets in the first hierarchy plot
#' @param color.use the character vector defining the color of each cell group
#' @param vertex.size The size of vertex
#' @param layout "hierarchy" or "circle"
#' @param thresh threshold of the p-value for determining significant interaction
#' @param pt.title font size of the text
#' @param title.space the space between the title and plot
#' @param vertex.label.cex The label size of vertex in the network
#' @importFrom grDevices dev.off pdf
#'
#' @return
#' @export
#'
#' @examples
#'
netVisual_aggregate <- function(object, signaling, signaling.name = NULL, vertex.receiver = NULL, color.use = NULL, vertex.size = 20, layout = c("hierarchy","circle"), thresh = 0.05,
                                pt.title = 12, title.space = 6, vertex.label.cex = 0.8) {
  layout <- match.arg(layout)
  pairLR <- searchPair(signaling = signaling, pairLR.use = object@LR$LRsig, key = "pathway_name", matching.exact = T, pair.only = T)

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
  nRow <- length(pairLR.name.use)

  prob <- prob[,,pairLR.name.use]
  pval <- pval[,,pairLR.name.use]

  if (length(dim(prob)) == 2) {
    prob <- replicate(1, prob, simplify="array")
    pval <- replicate(1, pval, simplify="array")
  }
  prob <-(prob-min(prob))/(max(prob)-min(prob))

  if (layout == "hierarchy") {
    prob.sum <- apply(prob, c(1,2), sum)
    prob.sum <-(prob.sum-min(prob.sum))/(max(prob.sum)-min(prob.sum))
    par(mfrow=c(1,2), ps = pt.title)
    netVisual_hierarchy1(prob.sum, vertex.receiver = vertex.receiver, color.use = color.use, vertex.size = vertex.size, signaling.name = NULL, vertex.label.cex = vertex.label.cex)
    netVisual_hierarchy2(prob.sum, vertex.receiver = setdiff(1:nrow(prob.sum),vertex.receiver), color.use = color.use, vertex.size = vertex.size, signaling.name = NULL, vertex.label.cex = vertex.label.cex)
    graphics::mtext(paste0(signaling.name, " signaling pathway network"), side = 3, outer = TRUE, cex = 1, line = -title.space)
  } else if (layout == "circle") {
    prob.sum <- apply(prob, c(1,2), sum)
    prob.sum <-(prob.sum-min(prob.sum))/(max(prob.sum)-min(prob.sum))
    netVisual_circle(prob.sum, top = 1, color.use = color.use, vertex.size = vertex.size, signaling.name = paste0(signaling.name, " signaling pathway network"), vertex.label.cex = vertex.label.cex)
  }


}



#' Visualize the inferred signaling network of individual L-R pairs
#'
#' @param object CellChat object
#' @param signaling a signaling pathway name
#' @param signaling.name alternative signaling pathway name to show on the plot
#' @param vertex.receiver a numeric vector giving the index of the cell groups as targets in the first hierarchy plot
#' @param color.use the character vector defining the color of each cell group
#' @param vertex.size The size of vertex
#' @param layout "hierarchy" or "circle"
#' @param height height of plot
#' @param thresh threshold of the p-value for determining significant interaction
#' @importFrom grDevices dev.off pdf
#'
#' @return
#' @export
#'
#' @examples
#'
netVisual_individual <- function(object, signaling, signaling.name = NULL, vertex.receiver = NULL, color.use = NULL, vertex.size = 20, layout = c("hierarchy","circle"), height = 5, thresh = 0.05) {
  layout <- match.arg(layout)
  pairLR <- searchPair(signaling = signaling, pairLR.use = object@LR$LRsig, key = "pathway_name", matching.exact = T, pair.only = T)

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

  nRow <- length(pairLR.name.use)

  prob <- prob[,,pairLR.name.use]
  pval <- pval[,,pairLR.name.use]

  if (length(dim(prob)) == 2) {
    prob <- replicate(1, prob, simplify="array")
    pval <- replicate(1, pval, simplify="array")
  }

  prob <-(prob-min(prob))/(max(prob)-min(prob))

  if (layout == "hierarchy") {
    par(mfrow=c(nRow,2), mar = c(5, 4, 4, 2) +0.1)
    for (i in 1:length(pairLR.name.use)) {
      signalName_i <- paste0(pairLR$ligand[i], "-",pairLR$receptor[i], sep = "")
      prob.i <- prob[,,i]
      netVisual_hierarchy1(prob.i, vertex.receiver = vertex.receiver, color.use = color.use, vertex.size = vertex.size, signaling.name = signalName_i)
      netVisual_hierarchy2(prob.i, vertex.receiver = setdiff(1:nrow(prob.i),vertex.receiver), color.use = color.use, vertex.size = vertex.size, signaling.name = signalName_i)
    }

  } else if (layout == "circle") {
    par(mfrow=c(nRow,1))
    for (i in 1:length(pairLR.name.use)) {
      signalName_i <- paste0(pairLR$ligand[i], "-",pairLR$receptor[i], sep = "")
      prob.i <- prob[,,i]
      netVisual_circle(prob.i, top = 1, color.use = color.use, vertex.size = vertex.size, signaling.name = signalName_i)
    }

  }
}




#' hierarchy plot of cell-cell communications sending to cell groups in vertex.receiver
#'
#' This function loads the significant interactions as a weighted matrix, and colors
#' represent different types of cells as a structure. The width of edges represent
#' the strength of the communication.
#'
#' @param net a weighted matrix defining the signaling network
#' @param vertex.receiver  a numeric vector giving the index of the cell groups as targets in the first hierarchy plot
#' @param color.use the character vector defining the color of each cell group
#' @param signaling.name alternative signaling pathway name to show on the plot
#' @param weight.scale whether rescale the edge weights
#' @param label.dist the distance between labels and dot position
#' @param space.v the space between different columns in the plot
#' @param space.h the space between different rows in the plot
#' @param label Whether or not shows the label of edges (number of connections
#' between different cell types)
#' @param edge.curved Specifies whether to draw curved edges, or not.
#' This can be a logical or a numeric vector or scalar.
#' First the vector is replicated to have the same length as the number of
#' edges in the graph. Then it is interpreted for each edge separately.
#' A numeric value specifies the curvature of the edge; zero curvature means
#' straight edges, negative values means the edge bends clockwise, positive
#' values the opposite. TRUE means curvature 0.5, FALSE means curvature zero
#' @param shape The shape of the vertex, currently “circle”, “square”,
#' “csquare”, “rectangle”, “crectangle”, “vrectangle”, “pie” (see
#' vertex.shape.pie), ‘sphere’, and “none” are supported, and only by the
#' plot.igraph command. “none” does not draw the vertices at all, although
#' vertex label are plotted (if given). See shapes for details about vertex
#' shapes and vertex.shape.pie for using pie charts as vertices.
#' @param vertex.size The size of vertex
#' @param margin The amount of empty space below, over, at the left and right
#'  of the plot, it is a numeric vector of length four. Usually values between
#'  0 and 0.5 are meaningful, but negative values are also possible, that will
#'  make the plot zoom in to a part of the graph. If it is shorter than four
#'  then it is recycled.
#' @param vertex.label.cex The label size of vertex
#' @param vertex.label.color The color of label for vertex
#' @param arrow.width The width of arrows
#' @param arrow.size the size of arrow
#' @param alpha.edge the transprency of edge
#' @param edge.label.color The color for single arrow
#' @param edge.label.cex The size of label for arrows
#' @param edge.max.width The maximum arrow size
#' @import igraph
#' @importFrom grDevices adjustcolor
#' @importFrom shape Arrows
#' @return A network graph of the significant interactions
#' @export
netVisual_hierarchy1 <-function(net, vertex.receiver, color.use = NULL, signaling.name = NULL, weight.scale = FALSE, label.dist = 2.8, space.v = 1.5, space.h = 1.6, shape= NULL, label=FALSE,edge.curved=0, vertex.size=20,margin=0.2,
                                vertex.label.cex=0.6,vertex.label.color= "black",arrow.width=1,arrow.size = 0.2,edge.label.color='black',edge.label.cex=0.5,edge.max.width=8,alpha.edge = 0.6){
  options(warn = -1)
  if (is.null(color.use)) {
    color.use <- scPalette(nrow(net))
  }
  vertex.size <- vertex.size/max(vertex.size)*15+7
  m <- length(vertex.receiver)
  net2 <- net
  net2 <- net2[,vertex.receiver]
  # Expand out to symmetric (M+N)x(M+N) matrix
  m1 <- nrow(net2); n1 <- ncol(net2)
  net3 <- rbind(cbind(matrix(0, m1, m1), net2), matrix(0, n1, m1+n1))

  row.names(net3) <- c(row.names(net)[vertex.receiver], row.names(net)[setdiff(1:m1,vertex.receiver)], rep("",m))
  colnames(net3) <- row.names(net3)
  color.use3 <- c(color.use[vertex.receiver], color.use[setdiff(1:m1,vertex.receiver)], rep("#FFFFFF",length(vertex.receiver)))
  color.use3.frame <- c(color.use[vertex.receiver], color.use[setdiff(1:m1,vertex.receiver)], color.use[vertex.receiver])

  if (length(vertex.size) != 1) {
    vertex.size = c(vertex.size[vertex.receiver], vertex.size[setdiff(1:m1,vertex.receiver)],vertex.size[vertex.receiver])
  }
  if (is.null(shape)) {
    shape <- c(rep("circle",m), rep("circle", m1-m), rep("circle",m))
  }

  g <- graph_from_adjacency_matrix(net3, mode = "directed", weighted = T)
  edge.start <- ends(g, es=E(g), names=FALSE)
  coords <- matrix(NA, nrow(net3), 2)
  coords[1:m,1] <- 0; coords[(m+1):m1,1] <- space.h; coords[(m1+1):nrow(net3),1] <- space.h/2;
  coords[1:m,2] <- seq(space.v, 0, by = -space.v/(m-1)); coords[(m+1):m1,2] <- seq(space.v, 0, by = -space.v/(m1-m-1));coords[(m1+1):nrow(net3),2] <- seq(space.v, 0, by = -space.v/(n1-1));
  coords_scale<-coords

  V(g)$size<-vertex.size
  V(g)$color<-color.use3[V(g)]
  V(g)$frame.color <- color.use3.frame[V(g)]
  V(g)$label.color <- vertex.label.color
  V(g)$label.cex<-vertex.label.cex
  if(label){
    E(g)$label<-E(g)$weight
  }
  if (weight.scale == TRUE) {
    E(g)$width<-0.3+edge.max.width/(max(E(g)$weight)-min(E(g)$weight))*(E(g)$weight-min(E(g)$weight))
  }else{
    E(g)$width<-0.3+edge.max.width*E(g)$weight
  }

  E(g)$arrow.width<-arrow.width
  E(g)$arrow.size<-arrow.size
  E(g)$label.color<-edge.label.color
  E(g)$label.cex<-edge.label.cex
  E(g)$color<-adjustcolor(V(g)$color[edge.start[,1]],alpha.edge)

  label.dist <- c(rep(space.h*label.dist,m), rep(space.h*label.dist, m1-m),rep(0, nrow(net3)-m1))
  label.locs <- c(rep(-pi, m), rep(0, m1-m),rep(-pi, nrow(net3)-m1))
  # text.pos <- cbind(c(-space.h/1.5, space.h/10, space.h/1.2), space.v-space.v/10)
  text.pos <- cbind(c(-space.h/1.5, space.h/22, space.h/1.5), space.v-space.v/7)
  add.vertex.shape("fcircle", clip=igraph.shape.noclip,plot=mycircle, parameters=list(vertex.frame.color=1, vertex.frame.width=1))
  plot(g,edge.curved=edge.curved,layout=coords_scale,margin=margin,rescale=T,vertex.shape="fcircle", vertex.frame.width = c(rep(1,m1), rep(2,nrow(net3)-m1)),
       vertex.label.degree=label.locs, vertex.label.dist=label.dist, vertex.label.family="Arial")
  text(text.pos, c("Source","Target","Source"), cex = 0.8, col = c("#c51b7d","#c51b7d","#2f6661"))
  arrow.pos1 <- c(-space.h/1.5, space.v-space.v/4, space.h/100000, space.v-space.v/4)
  arrow.pos2 <- c(space.h/1.5, space.v-space.v/4, space.h/20, space.v-space.v/4)
  shape::Arrows(arrow.pos1[1], arrow.pos1[2], arrow.pos1[3], arrow.pos1[4], col = "#c51b7d",arr.lwd = 0.0001,arr.length = 0.2, lwd = 0.8,arr.type="triangle")
  shape::Arrows(arrow.pos2[1], arrow.pos2[2], arrow.pos2[3], arrow.pos2[4], col = "#2f6661",arr.lwd = 0.0001,arr.length = 0.2, lwd = 0.8,arr.type="triangle")
  if (!is.null(signaling.name)) {
    title.pos = c(space.h/8, space.v)
    text(title.pos[1],title.pos[2],paste0(signaling.name, " signaling network"), cex = 1)
  }
}


#' Hierarchy plot of cell-cell communication sending to cell groups not in vertex.receiver
#'
#' This function loads the significant interactions as a weighted matrix, and colors
#' represent different types of cells as a structure. The width of edges represent the strength of the communication.
#'
#' @param net a weighted matrix defining the signaling network
#' @param vertex.receiver  a numeric vector giving the index of the cell groups as targets in the first hierarchy plot
#' @param color.use the character vector defining the color of each cell group
#' @param signaling.name alternative signaling pathway name to show on the plot
#' @param weight.scale whether rescale the edge weights
#' @param label.dist the distance between labels and dot position
#' @param space.v the space between different columns in the plot
#' @param space.h the space between different rows in the plot
#' @param label Whether or not shows the label of edges (number of connections between different cell types)
#' @param edge.curved Specifies whether to draw curved edges, or not.
#' This can be a logical or a numeric vector or scalar.
#' First the vector is replicated to have the same length as the number of
#' edges in the graph. Then it is interpreted for each edge separately.
#' A numeric value specifies the curvature of the edge; zero curvature means
#' straight edges, negative values means the edge bends clockwise, positive
#' values the opposite. TRUE means curvature 0.5, FALSE means curvature zero
#' @param shape The shape of the vertex, currently “circle”, “square”,
#' “csquare”, “rectangle”, “crectangle”, “vrectangle”, “pie” (see
#' vertex.shape.pie), ‘sphere’, and “none” are supported, and only by the
#' plot.igraph command. “none” does not draw the vertices at all, although
#' vertex label are plotted (if given). See shapes for details about vertex
#' shapes and vertex.shape.pie for using pie charts as vertices.
#' @param vertex.size The size of vertex
#' @param margin The amount of empty space below, over, at the left and right
#'  of the plot, it is a numeric vector of length four. Usually values between
#'  0 and 0.5 are meaningful, but negative values are also possible, that will
#'  make the plot zoom in to a part of the graph. If it is shorter than four
#'  then it is recycled.
#' @param vertex.label.cex The label size of vertex
#' @param vertex.label.color The color of label for vertex
#' @param arrow.width The width of arrows
#' @param arrow.size the size of arrow
#' @param alpha.edge the transprency of edge
#' @param edge.label.color The color for single arrow
#' @param edge.label.cex The size of label for arrows
#' @param edge.max.width The maximum arrow size
#' @import igraph
#' @importFrom shape Arrows
#' @return A network graph of the significant interactions
#' @export
netVisual_hierarchy2 <-function(net, vertex.receiver, color.use = NULL, signaling.name = NULL, weight.scale = FALSE, label.dist = 2.8, space.v = 1.5, space.h = 1.6, shape= NULL, label=FALSE,edge.curved=0, vertex.size=20,margin=0.2,
                                vertex.label.cex=0.6,vertex.label.color= "black",arrow.width=1,arrow.size = 0.2,edge.label.color='black',edge.label.cex=0.5,edge.max.width=8,alpha.edge = 0.6){
  options(warn = -1)
  if (is.null(color.use)) {
    color.use <- scPalette(nrow(net))
  }
  vertex.size <- vertex.size/max(vertex.size)*15+6
  m <- length(vertex.receiver)
  m0 <- nrow(net)-length(vertex.receiver)
  net2 <- net
  net2 <- net2[,vertex.receiver]
  # Expand out to symmetric (M+N)x(M+N) matrix
  m1 <- nrow(net2); n1 <- ncol(net2)
  net3 <- rbind(cbind(matrix(0, m1, m1), net2), matrix(0, n1, m1+n1))
  row.names(net3) <- c(row.names(net)[setdiff(1:m1,vertex.receiver)],row.names(net)[vertex.receiver],  rep("",m))
  colnames(net3) <- row.names(net3)
  color.use3 <- c(color.use[setdiff(1:m1,vertex.receiver)],color.use[vertex.receiver],  rep("#FFFFFF",length(vertex.receiver)))
  color.use3.frame <- c(color.use[setdiff(1:m1,vertex.receiver)], color.use[vertex.receiver], color.use[vertex.receiver])


  if (length(vertex.size) != 1) {
    vertex.size = c(vertex.size[setdiff(1:m1,vertex.receiver)], vertex.size[vertex.receiver], vertex.size[vertex.receiver])
  }
  if (is.null(shape)) {
    shape <- rep("circle",nrow(net3))
  }

  g <- graph_from_adjacency_matrix(net3, mode = "directed", weighted = T)
  edge.start <- ends(g, es=E(g), names=FALSE)
  coords <- matrix(NA, nrow(net3), 2)
  coords[1:m0,1] <- 0; coords[(m0+1):m1,1] <- space.h; coords[(m1+1):nrow(net3),1] <- space.h/2;
  coords[1:m0,2] <- seq(space.v, 0, by = -space.v/(m0-1)); coords[(m0+1):m1,2] <- seq(space.v, 0, by = -space.v/(m1-m0-1));coords[(m1+1):nrow(net3),2] <- seq(space.v, 0, by = -space.v/(n1-1));
  coords_scale<-coords

  V(g)$size<-vertex.size
  V(g)$color<-color.use3[V(g)]
  V(g)$frame.color <- color.use3.frame[V(g)]
  V(g)$label.color <- vertex.label.color
  V(g)$label.cex<-vertex.label.cex
  if(label){
    E(g)$label<-E(g)$weight
  }
  if (weight.scale == TRUE) {
    E(g)$width<-0.3+edge.max.width/(max(E(g)$weight)-min(E(g)$weight))*(E(g)$weight-min(E(g)$weight))
  }else{
    E(g)$width<-0.3+edge.max.width*E(g)$weight
  }
  E(g)$arrow.width<-arrow.width
  E(g)$arrow.size<-arrow.size
  E(g)$label.color<-edge.label.color
  E(g)$label.cex<-edge.label.cex
  E(g)$color<-adjustcolor(V(g)$color[edge.start[,1]],alpha.edge)

  label.dist <- c(rep(space.h*label.dist,m), rep(space.h*label.dist, m1-m),rep(0, nrow(net3)-m1))
  label.locs <- c(rep(-pi, m0), rep(0, m1-m0),rep(-pi, nrow(net3)-m1))
  #text.pos <- cbind(c(-space.h/1.5, space.h/10, space.h/1.2), space.v-space.v/10)
  text.pos <- cbind(c(-space.h/1.5, space.h/22, space.h/1.5), space.v-space.v/7)
  add.vertex.shape("fcircle", clip=igraph.shape.noclip,plot=mycircle, parameters=list(vertex.frame.color=1, vertex.frame.width=1))
  plot(g,edge.curved=edge.curved,layout=coords_scale,margin=margin,rescale=T,vertex.shape="fcircle", vertex.frame.width = c(rep(1,m1), rep(2,nrow(net3)-m1)),
       vertex.label.degree=label.locs, vertex.label.dist=label.dist, vertex.label.family="Arial")
  text(text.pos, c("Source","Target","Source"), cex = 0.8, col = c("#c51b7d","#2f6661","#2f6661"))

  arrow.pos1 <- c(-space.h/1.5, space.v-space.v/4, space.h/100000, space.v-space.v/4)
  arrow.pos2 <- c(space.h/1.5, space.v-space.v/4, space.h/20, space.v-space.v/4)
  shape::Arrows(arrow.pos1[1], arrow.pos1[2], arrow.pos1[3], arrow.pos1[4], col = "#c51b7d",arr.lwd = 0.0001,arr.length = 0.2, lwd = 0.8,arr.type="triangle")
  shape::Arrows(arrow.pos2[1], arrow.pos2[2], arrow.pos2[3], arrow.pos2[4], col = "#2f6661",arr.lwd = 0.0001,arr.length = 0.2, lwd = 0.8,arr.type="triangle")

  if (!is.null(signaling.name)) {
    title.pos = c(space.h/8, space.v)
    text(title.pos[1],title.pos[2],paste0(signaling.name, " signaling network"), cex = 1)
  }
}



#' Circle plot of cell-cell communication network
#'
#' The width of edges represent the strength of the communication.
#'
#' @param net A weighted matrix representing the connections
#' @param color.use Colors represent different cell groups
#' @param signaling.name the name of the title
#' @param top the fraction of interactions to show
#' @param weight.scale whether scale the weight
#' @param label.edge Whether or not shows the label of edges
#' @param edge.curved Specifies whether to draw curved edges, or not.
#' This can be a logical or a numeric vector or scalar.
#' First the vector is replicated to have the same length as the number of
#' edges in the graph. Then it is interpreted for each edge separately.
#' A numeric value specifies the curvature of the edge; zero curvature means
#' straight edges, negative values means the edge bends clockwise, positive
#' values the opposite. TRUE means curvature 0.5, FALSE means curvature zero
#' @param shape The shape of the vertex, currently “circle”, “square”,
#' “csquare”, “rectangle”, “crectangle”, “vrectangle”, “pie” (see
#' vertex.shape.pie), ‘sphere’, and “none” are supported, and only by the
#' plot.igraph command. “none” does not draw the vertices at all, although
#' vertex label are plotted (if given). See shapes for details about vertex
#' shapes and vertex.shape.pie for using pie charts as vertices.
#' @param layout The layout specification. It must be a call to a layout
#' specification function.
#' @param vertex.size The size of vertex
#' @param margin The amount of empty space below, over, at the left and right
#'  of the plot, it is a numeric vector of length four. Usually values between
#'  0 and 0.5 are meaningful, but negative values are also possible, that will
#'  make the plot zoom in to a part of the graph. If it is shorter than four
#'  then it is recycled.
#' @param vertex.label.cex The label size of vertex
#' @param vertex.label.color The color of label for vertex
#' @param arrow.width The width of arrows
#' @param arrow.size the size of arrow
#' @param alpha.edge the transprency of edge
#' @param edge.label.color The color for single arrow
#' @param edge.label.cex The size of label for arrows
#' @param edge.max.width The maximum arrow size
#' @import igraph
#' @return A network graph of the significant interactions
#' @export
netVisual_circle <-function(net, color.use = NULL,signaling.name = NULL, top = 1, weight.scale = FALSE,label.edge = FALSE,edge.curved=0.2,shape='circle',layout=in_circle(),vertex.size=20,margin=0.2,
                            vertex.label.cex=0.8,vertex.label.color= "black",arrow.width=1,arrow.size = 0.2,edge.label.color='black',edge.label.cex=0.5,edge.max.width=8,alpha.edge = 0.6){
  options(warn = -1)
  thresh <- stats::quantile(net, probs = 1-top)
  net[net < thresh] <- 0
  g <- graph_from_adjacency_matrix(net, mode = "directed", weighted = T)
  edge.start <- igraph::ends(g, es=E(g), names=FALSE)
  coords<-layout_(g,layout)
  if(nrow(coords)!=1){
    coords_scale=scale(coords)
  }else{
    coords_scale<-coords
  }
  if (is.null(color.use)) {
    color.use = scPalette(length(V(g)))
  }

  vertex.size <- vertex.size/max(vertex.size)*15+5

  loop.angle<-ifelse(coords_scale[V(g),1]>0,-atan(coords_scale[V(g),2]/coords_scale[V(g),1]),pi-atan(coords_scale[V(g),2]/coords_scale[V(g),1]))
  V(g)$size<-vertex.size
  V(g)$color<-color.use[V(g)]
  V(g)$frame.color <- color.use[V(g)]
  V(g)$label.color <- vertex.label.color
  V(g)$label.cex<-vertex.label.cex
  if(label.edge){
    E(g)$label<-E(g)$weight
  }
  if (weight.scale == TRUE) {
    E(g)$width<-0.3+edge.max.width/(max(E(g)$weight)-min(E(g)$weight))*(E(g)$weight-min(E(g)$weight))
  }else{
    E(g)$width<-0.3+edge.max.width*E(g)$weight
  }

  E(g)$arrow.width<-arrow.width
  E(g)$arrow.size<-arrow.size
  E(g)$label.color<-edge.label.color
  E(g)$label.cex<-edge.label.cex
  E(g)$color<- grDevices::adjustcolor(V(g)$color[edge.start[,1]],alpha.edge)

  if(sum(edge.start[,2]==edge.start[,1])!=0){
    E(g)$loop.angle[which(edge.start[,2]==edge.start[,1])]<-loop.angle[edge.start[which(edge.start[,2]==edge.start[,1]),1]]
  }
  radian.rescale <- function(x, start=0, direction=1) {
    c.rotate <- function(x) (x + start) %% (2 * pi) * direction
    c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
  }
  label.locs <- radian.rescale(x=1:length(V(g)), direction=-1, start=0)
  label.dist <- vertex.size/max(vertex.size)+2
  plot(g,edge.curved=edge.curved,vertex.shape=shape,layout=coords_scale,margin=margin, vertex.label.dist=label.dist,
       vertex.label.degree=label.locs, vertex.label.family="Arial") # "sans"
  if (!is.null(signaling.name)) {
    text(0,1.5,signaling.name, cex = 1.1)
  }

}


#' generate circle symbol
#'
#' @param coords coordinates of points
#' @param v vetex
#' @param params parameters
#' @importFrom graphics symbols
#' @return
mycircle <- function(coords, v=NULL, params) {
  vertex.color <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.size  <- 1/200 * params("vertex", "size")
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size <- vertex.size[v]
  }
  vertex.frame.color <- params("vertex", "frame.color")
  if (length(vertex.frame.color) != 1 && !is.null(v)) {
    vertex.frame.color <- vertex.frame.color[v]
  }
  vertex.frame.width <- params("vertex", "frame.width")
  if (length(vertex.frame.width) != 1 && !is.null(v)) {
    vertex.frame.width <- vertex.frame.width[v]
  }

  mapply(coords[,1], coords[,2], vertex.color, vertex.frame.color,
         vertex.size, vertex.frame.width,
         FUN=function(x, y, bg, fg, size, lwd) {
           symbols(x=x, y=y, bg=bg, fg=fg, lwd=lwd,
                   circles=size, add=TRUE, inches=FALSE)
         })
}


#' Heatmap showing the centrality scores/importance of cell groups as senders, receivers, mediators and influencers in the intercellular communication network
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
netVisual_signalingRole <- function(object, signaling, slot.name = "netP", measure = c("outdeg","indeg","flowbet","info"), measure.name = c("Sender","Receiver","Mediator","Influencer"),
                             color.use = NULL, color.heatmap = "BuGn",
                             width = 6.5, height = 1.4, font.size = 8, font.size.title = 10, cluster.rows = FALSE, cluster.cols = FALSE) {
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


#' Show all the significant interactions (L-R pairs) from some cell groups to other cell groups
#'
#' @param object CellChat object
#' @param from a vector giving the index or the name of source cell groups
#' @param to a corresponding vector giving the index or the name of target cell groups. Note: The length of 'from' and 'to' must be the same, giving the corresponding pair of cell groups for communication.
#' @param bidirection whether show the bidirectional communication, i.e., both 'from'->'to' and 'to'->'from'.
#' @param pairLR.use0 ligand-receptor pairs to use; default is all the significant interactions
#' @param color.heatmap color map
#' @param thresh threshold of the p-value for determining significant interaction
#'
#' @return
#' @export
#'
netVisual_bubble <- function(object, from, to, bidirection = FALSE, pairLR.use0 = NULL,  color.heatmap = viridis::viridis(50), thresh = 0.05){
  pairwiseLR <- object@net$pairwiseRank
  if (is.null(pairwiseLR)) {
    stop("The interactions between pairwise cell groups have not been extracted!
         Please first run `object <- rankNetPairwise(object)`")
  }
  group.names.all <- names(pairwiseLR)
  if (!is.numeric(from)) {
    from <- match(from, group.names.all)
    if (sum(is.na(from)) > 0) {
      message("Some input cell group names in 'from' do not exist!")
      from <- from[!is.na(from)]
    }
  }
  if (!is.numeric(to)) {
    to <- match(to, group.names.all)
    if (sum(is.na(to)) > 0) {
      message("Some input cell group names in 'to' do not exist!")
      to <- to[!is.na(to)]
    }
  }
  if (length(from) != length(to)) {
    stop("The length of 'from' and 'to' must be the same!")
  }
  if (bidirection) {
    from <- c(from, to)
    to <- c(to, from)
  }
  if (is.null(pairLR.use0)) {
    k <- 0
    pairLR.use0 <- list()
    for (i in 1:length(from)){
        pairwiseLR_ij <- pairwiseLR[[from[i]]][[to[i]]]
        idx <- pairwiseLR_ij$pval < thresh
        if (length(idx) > 0) {
          k <- k +1
          pairLR.use0[[k]] <- pairwiseLR_ij[idx,]
        }
    }
    pairLR.use0 <- do.call(rbind, pairLR.use0)
  }

  k <- 0
  pval <- matrix(nrow = length(rownames(pairLR.use0)), ncol = length(from))
  prob <- pval
  group.names <- c()
  for (i in 1:length(from)) {
      k <- k+1
      pairwiseLR_ij <- pairwiseLR[[from[i]]][[to[i]]]
      pairwiseLR_ij <- pairwiseLR_ij[rownames(pairLR.use0),]
      pval_ij <- pairwiseLR_ij$pval
      prob_ij <- pairwiseLR_ij$prob
      pval_ij[pval_ij > 0.05] = 1
      pval_ij[pval_ij > 0.01 & pval_ij <= 0.05] = 2
      pval_ij[pval_ij <= 0.01] = 3
      pval[,k] <- pval_ij
      prob[,k] <- prob_ij
      group.names <- c(group.names, paste(group.names.all[from[i]], group.names.all[to[i]], sep = " - "))
  }
  prob[which(prob == 0)] <- NA
  # remove rows that are entirely NA
  pval <- pval[rowSums(is.na(prob)) != ncol(prob), ,drop = FALSE]
  pairLR.use0 <- pairLR.use0[rowSums(is.na(prob)) != ncol(prob), ,drop = FALSE]
  prob <- prob[rowSums(is.na(prob)) != ncol(prob), ,drop = FALSE]

  prob <- -1/log(prob)
  prob[is.infinite(prob) & prob < 0] <- max(prob+1)

  rownames(prob) <- as.character(pairwiseLR_ij[rownames(pairLR.use0),]$interaction_name_2)
  colnames(prob) <- group.names
  dimnames(pval) <- dimnames(prob)

  df.prob = reshape2::melt(prob, value.name = "probs")
  df.pval = reshape2::melt(pval, value.name = "pvals")
  df <- data.frame(cbind(df.prob, pvals = df.pval$pvals))
  df$Var1 <- as.character(df$Var1)
  df <- with(df, df[order(Var1),])

  g <- ggplot(df, aes(x = Var2, y = Var1, color = probs, size = pvals)) +
    geom_point(pch = 16) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust=1),
          axis.title.x = element_blank(),
          axis.title.y = element_blank()) +
    scale_x_discrete(position = "bottom") +
    scale_radius(range = c(1,3), breaks = c(1,2,3),labels = c("p > 0.05", "0.01 < p < 0.05","p < 0.01"), name = "p-value")

  if(length(color.heatmap) == 1){
    g <- g + scale_colour_gradientn(colors = colorRampPalette(c("white", color.heatmap))(100), na.value = "white",limits=c(0, quantile(df$probs, 1,na.rm= T)),
                                    breaks = c(0, quantile(df$probs, 1,na.rm= T)), labels = c("0","max")) +
      guides(fill = guide_colourbar(barwidth = 0.1, title = "Commun. Prob."))
  } else {
    g <- g + scale_colour_gradientn(colors = colorRampPalette(color.heatmap)(99), na.value = "white",limits=c(0, quantile(df$probs, 1,na.rm= T)),
                                    breaks = c(0, quantile(df$probs, 1,na.rm= T)), labels = c("0","max")) +
      guides(color = guide_colourbar(barwidth = 0.5, title = "Commun. Prob."))
  }
  g <- g + theme(legend.title = element_text(size = 8), legend.text = element_text(size = 6))
  g
  return(g)

}


#' River plot showing the associations of latent patterns with cell groups and ligand-receptor pairs or signaling pathways
#'
#' @param object CellChat object
#' @param slot.name the slot name of object that is used to compute centrality measures of signaling networks
#' @param pattern "outgoing" or "incoming"
#' @param cutoff the threshold for filtering out weak links
#' @param color.use the character vector defining the color of each cell group
#' @param color.use.pattern the character vector defining the color of each pattern
#' @param color.use.signaling the character vector defining the color of each signaling
#' @param do.order whether reorder the cell groups or signaling according to their similarity
#' @param main.title the title of plot
#' @param font.size font size of the text
#' @param font.size.title font size of the title
#' @importFrom methods slot
#' @importFrom stats cutree dist hclust
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#' @import ggalluvial
# #' @importFrom ggalluvial geom_stratum geom_flow to_lodes_form
#' @importFrom ggplot2 geom_text scale_x_discrete scale_fill_manual theme ggtitle
#' @importFrom cowplot plot_grid ggdraw draw_label
#' @return
#' @export
#'
#' @examples
netAnalysis_river <- function(object, slot.name = "netP", pattern = c("outgoing","incoming"), cutoff = 0.5,
                              color.use = NULL, color.use.pattern = NULL, color.use.signaling = "grey50",
                              do.order = FALSE, main.title = NULL,
                              font.size = 2.5, font.size.title = 12){
  requireNamespace("ggalluvial")
  res.pattern <- methods::slot(object, slot.name)$pattern[[pattern]]
  data1 = res.pattern$pattern$cell
  data2 = res.pattern$pattern$signaling
  if (is.null(color.use.pattern)) {
    nPatterns <- length(unique(data1$Pattern))
    if (pattern == "outgoing") {
      color.use.pattern = ggPalette(nPatterns*2)[seq(1,nPatterns*2, by = 2)]
    } else if (pattern == "incoming") {
      color.use.pattern = ggPalette(nPatterns*2)[seq(2,nPatterns*2, by = 2)]
    }
  }
  if (is.null(main.title)) {
    if (pattern == "outgoing") {
      main.title = "Outgoing communication patterns of secreting cells"
    } else if (pattern == "incoming") {
      main.title = "Incoming communication patterns of target cells"
    }
  }

  if (is.null(data2)) {
    data1$Contribution[data1$Contribution < cutoff] <- 0
    plot.data <- data1
    nPatterns<-length(unique(plot.data$Pattern))
    nCellGroup<-length(unique(plot.data$CellGroup))
    if (is.null(color.use)) {
      color.use <- scPalette(nCellGroup)
    }
    if (is.null(color.use.pattern)){
      color.use.pattern <- ggPalette(nPatterns)
    }

    plot.data.long <- to_lodes_form(plot.data, axes = 1:2, id = "connection")
    if (do.order) {
      mat = tapply(plot.data[["Contribution"]], list(plot.data[["CellGroup"]], plot.data[["Pattern"]]), sum)
      d <- dist(as.matrix(mat))
      hc <- hclust(d, "ave")
      k <- length(unique(grep("Pattern", plot.data.long$stratum[plot.data.long$Contribution != 0], value = T)))
      cluster <- hc %>% cutree(k)
      order.name <- order(cluster)
      plot.data.long$stratum <- factor(plot.data.long$stratum, levels = c(names(cluster)[order.name], colnames(mat)))
      color.use <- color.use[order.name]
    }
    color.use.all <- c(color.use, color.use.pattern)
    gg <- ggplot(plot.data.long,aes(x = factor(x, levels = c("CellGroup", "Pattern")),y=Contribution,
                                    stratum = stratum, alluvium = connection,
                                    fill = stratum, label = stratum)) +
      geom_flow(width = 1/3,aes.flow = "backward") +
      geom_stratum(width=1/3,size=0.1,color="black", alpha = 0.8, linetype = 1) +
      geom_text(stat = "stratum", size = font.size) +
      scale_x_discrete(limits = c(),  labels=c("Cell groups", "Patterns")) +
      scale_fill_manual(values = alpha(color.use.all, alpha = 0.8), drop = FALSE) +
      theme_bw()+
      theme(legend.position = "none",
            axis.title = element_blank(),
            axis.text.y= element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor  = element_blank(),
            panel.border = element_blank(),
            axis.ticks = element_blank(),axis.text=element_text(size=10))+
      ggtitle(main.title)

  } else {
    data1$Contribution[data1$Contribution < cutoff] <- 0
    plot.data <- data1
    nPatterns<-length(unique(plot.data$Pattern))
    nCellGroup<-length(unique(plot.data$CellGroup))
    if (is.null(color.use)) {
      color.use <- scPalette(nCellGroup)
      # color.use <- grDevices::colorRampPalette((RColorBrewer::brewer.pal(n = 9, name = "Set1")))(nCellGroup)
    }
    if (is.null(color.use.pattern)){
      color.use.pattern <- ggPalette(nPatterns)
    }

    ## connect cell groups with patterns
    plot.data.long <- to_lodes_form(plot.data, axes = 1:2, id = "connection")
    if (do.order) {
      mat = tapply(plot.data[["Contribution"]], list(plot.data[["CellGroup"]], plot.data[["Pattern"]]), sum)
      d <- dist(as.matrix(mat))
      hc <- hclust(d, "ave")
      k <- length(unique(grep("Pattern", plot.data.long$stratum[plot.data.long$Contribution != 0], value = T)))
      cluster <- hc %>% cutree(k)
      order.name <- order(cluster)
      plot.data.long$stratum <- factor(plot.data.long$stratum, levels = c(names(cluster)[order.name], colnames(mat)))
      color.use <- color.use[order.name]
    }
    color.use.all <- c(color.use, color.use.pattern)
    StatStratum <- ggalluvial::StatStratum
    gg1 <- ggplot(plot.data.long,aes(x = factor(x, levels = c("CellGroup", "Pattern")),y=Contribution,
                                     stratum = stratum, alluvium = connection,
                                     fill = stratum, label = stratum)) +
      geom_flow(width = 1/3,aes.flow = "backward") +
      geom_stratum(width=1/3,size=0.1,color="black", alpha = 0.8, linetype = 1) +
      geom_text(stat = "stratum", size = font.size) +
      scale_x_discrete(limits = c(),  labels=c("Cell groups", "Patterns")) +
      scale_fill_manual(values = alpha(color.use.all, alpha = 0.8), drop = FALSE) +
      theme_bw()+
      theme(legend.position = "none",
            axis.title = element_blank(),
            axis.text.y= element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor  = element_blank(),
            panel.border = element_blank(),
            axis.ticks = element_blank(),axis.text=element_text(size=10)) +
      theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

    ## connect patterns with signaling
    data2$Contribution[data2$Contribution < cutoff] <- 0
    plot.data <- data2
    nPatterns<-length(unique(plot.data$Pattern))
    nSignaling<-length(unique(plot.data$Signaling))
    if (length(color.use.signaling) == 1) {
      color.use.all <- c(color.use.pattern, rep(color.use.signaling, nSignaling))
    } else {
      color.use.all <- c(color.use.pattern, color.use.signaling)
    }


    plot.data.long <- ggalluvial::to_lodes_form(plot.data, axes = 1:2, id = "connection")
    if (do.order) {
      mat = tapply(plot.data[["Contribution"]], list(plot.data[["Signaling"]], plot.data[["Pattern"]]), sum)
      d <- dist(as.matrix(mat))
      hc <- hclust(d, "ave")
      k <- length(unique(grep("Pattern", plot.data.long$stratum[plot.data.long$Contribution != 0], value = T)))
      cluster <- hc %>% cutree(k)
      order.name <- order(cluster)
      plot.data.long$stratum <- factor(plot.data.long$stratum, levels = c(colnames(mat),names(cluster)[order.name]))
    }

    gg2 <- ggplot(plot.data.long,aes(x = factor(x, levels = c("Pattern", "Signaling")),y= Contribution,
                                     stratum = stratum, alluvium = connection,
                                     fill = stratum, label = stratum)) +
      geom_flow(width = 1/3,aes.flow = "forward") +
      geom_stratum(width=1/3,size=0.1,color="black", alpha = 0.8, linetype = 1) +
      geom_text(stat = "stratum", size = font.size) + # 2.5
      scale_x_discrete(limits = c(),  labels=c("Patterns", "Signaling")) +
      scale_fill_manual(values = alpha(color.use.all, alpha = 0.8), drop = FALSE) +
      theme_bw()+
      theme(legend.position = "none",
            axis.title = element_blank(),
            axis.text.y= element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor  = element_blank(),
            panel.border = element_blank(),
            axis.ticks = element_blank(),axis.text=element_text(size= 10))+
      theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

    ## connect cell groups with signaling
    # data1 = data1[data1$Contribution > 0,]
    # data2 = data2[data2$Contribution > 0,]
    data3 = merge(data1, data2, by.x="Pattern", by.y="Pattern")
    data3$Contribution <- data3$Contribution.x * data3$Contribution.y
    data3 <- data3[,colnames(data3) %in% c("CellGroup","Signaling","Contribution")]

    # plot.data <- data3
    # nSignaling<-length(unique(plot.data$Signaling))
    # nCellGroup<-length(unique(plot.data$CellGroup))
    #
    # if (length(color.use.signaling) == 1) {
    #   color.use.signaling <- rep(color.use.signaling, nSignaling)
    # }
    #
    #
    # ## connect cell groups with patterns
    # plot.data.long <- to_lodes_form(plot.data, axes = 1:2, id = "connection")
    # if (do.order) {
    #   mat = tapply(plot.data[["Contribution"]], list(plot.data[["CellGroup"]], plot.data[["Signaling"]]), sum)
    #   d <- dist(as.matrix(mat))
    #   hc <- hclust(d, "ave")
    #   k <- length(unique(grep("Signaling", plot.data.long$stratum[plot.data.long$Contribution != 0], value = T)))
    #   cluster <- hc %>% cutree(k)
    #   order.name <- order(cluster)
    #   plot.data.long$stratum <- factor(plot.data.long$stratum, levels = c(names(cluster)[order.name], colnames(mat)))
    #   color.use <- color.use[order.name]
    # }
    # color.use.all <- c(color.use, color.use.signaling)

    # gg3 <- ggplot(plot.data.long, aes(x = factor(x, levels = c("CellGroup", "Signaling")),y=Contribution,
    #                                  stratum = stratum, alluvium = connection,
    #                                  fill = stratum, label = stratum)) +
    #   geom_flow(width = 1/3,aes.flow = "forward") +
    #   geom_stratum(width=1/3,size=0.1,color="black", alpha = 0.8, linetype = 1) +
    #   geom_text(stat = "stratum", size = 2.5) +
    #   scale_x_discrete(limits = c(),  labels=c("Cell groups", "Signaling")) +
    #   scale_fill_manual(values = alpha(color.use.all, alpha = 0.8), drop = FALSE) +
    #   theme_bw()+
    #   theme(legend.position = "none",
    #         axis.title = element_blank(),
    #         axis.text.y= element_blank(),
    #         panel.grid.major = element_blank(),
    #         panel.grid.minor  = element_blank(),
    #         panel.border = element_blank(),
    #         axis.ticks = element_blank(),axis.text=element_text(size=10)) +
    #   theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))


    gg <- cowplot::plot_grid(gg1, gg2,align = "h", nrow = 1)
    title <- cowplot::ggdraw() + cowplot::draw_label(main.title,size = font.size.title)
    gg <- cowplot::plot_grid(title, gg, ncol=1, rel_heights=c(0.1, 1))
  }
  return(gg)
}


#' Dot plots showing the associations of latent patterns with cell groups and ligand-receptor pairs or signaling pathways
#'
#' @param object CellChat object
#' @param slot.name the slot name of object that is used to compute centrality measures of signaling networks
#' @param pattern "outgoing" or "incoming"
#' @param cutoff the threshold for filtering out weak links
#' @param color.use the character vector defining the color of each cell group
#' @param pathway.show the character vector defining the signaling to show
#' @param group.show the character vector defining the cell group to show
#' @param shape the shape of the symbol: 21 for circle and 22 for square
#' @param dot.size a range defining the size of the symbol
#' @param dot.alpha transparency
#' @param main.title the title of plot
#' @param font.size font size of the text
#' @param font.size.title font size of the title
#' @importFrom methods slot
#' @import ggplot2
#' @importFrom dplyr group_by top_n
#' @return
#' @export
#'
#' @examples
netAnalysis_dot <- function(object, slot.name = "netP", pattern = c("outgoing","incoming"), cutoff = NULL, color.use = NULL,
                            pathway.show = NULL, group.show = NULL,
                            shape = 21, dot.size = c(1, 3), dot.alpha = 1, main.title = NULL,
                            font.size = 10, font.size.title = 12){
  pattern <- match.arg(pattern)
  patternSignaling <- methods::slot(object, slot.name)$pattern[[pattern]]
  data1 = patternSignaling$pattern$cell
  data2 = patternSignaling$pattern$signaling
  data = patternSignaling$data
  if (is.null(main.title)) {
    if (pattern == "outgoing") {
      main.title = "Outgoing communication patterns of secreting cells"
    } else if (pattern == "incoming") {
      main.title = "Incoming communication patterns of target cells"
    }
  }
  if (is.null(color.use)) {
    color.use <- scPalette(nlevels(data1$CellGroup))
  }
  if (is.null(cutoff)) {
    cutoff <- 1/length(unique(data1$Pattern))
  }
  options(warn = -1)
  data1$Contribution[data1$Contribution < cutoff] <- 0
  data2$Contribution[data2$Contribution < cutoff] <- 0
  data3 = merge(data1, data2, by.x="Pattern", by.y="Pattern")
  data3$Contribution <- data3$Contribution.x * data3$Contribution.y
  data3 <- data3[,colnames(data3) %in% c("CellGroup","Signaling","Contribution")]
  if (!is.null(pathway.show)) {
    data3 <- data3[data3$Signaling %in% pathway.show, ]
    pathway.add <- pathway.show[which(pathway.show %in% data3$Signaling == 0)]
    if (length(pathway.add) > 1) {
      data.add <- expand.grid(CellGroup = levels(data1$CellGroup), Signaling = pathway.add)
      data.add$Contribution <- 0
      data3 <- rbind(data3, data.add)
    }
    data3$Signaling <- factor(data3$Signaling, levels = pathway.show)
  }
  if (!is.null(group.show)) {
    data3$CellGroup <- as.character(data3$CellGroup)
    data3 <- data3[data3$CellGroup %in% group.show, ]
    data3$CellGroup <- factor(data3$CellGroup, levels = group.show)
  }

  data <- as.data.frame(as.table(data));
  data <- data[data[,3] != 0, ]
  data12 <- paste0(data[,1],data[,2])
  data312 <- paste0(data3[,1],data3[,2])
  idx1 <- which(match(data312, data12, nomatch = 0) ==0)
  data3$Contribution[idx1] <- 0
  data3$id <- data312
  data3 <- data3 %>% group_by(id) %>% top_n(1, Contribution)

  data3$Contribution[which(data3$Contribution == 0)] <- NA

  df <- data3
  gg <- ggplot(data = df, aes(x = Signaling, y = CellGroup)) +
    geom_point(aes(size =  Contribution, fill = CellGroup, colour = CellGroup), shape = shape) +
    scale_size_continuous(range = dot.size) +
    theme_linedraw() +
    scale_x_discrete(position = "bottom") +
    ggtitle(main.title) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(text = element_text(size = font.size),plot.title = element_text(size=font.size.title, face="plain"),
          axis.text.x = element_text(angle = 45, hjust=1),
          axis.text.y = element_text(angle = 0, hjust=1),
          axis.title.x = element_blank(),
          axis.title.y = element_blank()) +
    theme(axis.line.x = element_line(size = 0.25), axis.line.y = element_line(size = 0.25)) +
    theme(panel.grid.major = element_line(colour="grey90", size = (0.1)))
  gg <- gg + scale_y_discrete(limits = rev(levels(data3$CellGroup)))
  gg <- gg + scale_fill_manual(values = ggplot2::alpha(color.use, alpha = dot.alpha), drop = FALSE, na.value = "white")
  gg <- gg + scale_colour_manual(values = color.use, drop = FALSE, na.value = "white")
  gg <- gg + guides(colour=FALSE) + guides(fill=FALSE)
  gg <- gg + theme(legend.title = element_text(size = 10), legend.text = element_text(size = 8))
  gg
  return(gg)
}


#' 2D visualization of the learned manifold of signaling networks
#'
#' @param object CellChat object
#' @param slot.name the slot name of object that is used to compute centrality measures of signaling networks
#' @param type "functional","structural"
#' @param pathway.remove a character vector defining the signaling to remove
#' @param pathway.remove.show whether show the removed signaling names
#' @param color.use defining the color for each cell group
#' @param dot.size a range defining the size of the symbol
#' @param dot.alpha transparency
#' @param xlabel label of x-axis
#' @param ylabel label of y-axis
#' @param title main title of the plot
#' @param font.size font size of the text
#' @param font.size.title font size of the title
#' @param label.size font size of the text
#' @param do.label label the each point
#' @param show.legend whether show the legend
#' @param show.axes whether show the axes
#' @import ggplot2
#' @importFrom ggrepel geom_text_repel
#' @importFrom methods slot
#' @return
#' @export
#'
#' @examples
netVisual_embedding <- function(object, slot.name = "netP", type = c("functional","structural"), color.use = NULL, pathway.remove = NULL, pathway.remove.show = TRUE, dot.size = c(2, 6), label.size = 2, dot.alpha = 0.5,
                                xlabel = "Dim 1", ylabel = "Dim 2", title = NULL,
                                font.size = 10, font.size.title = 12, do.label = T, show.legend = T, show.axes = T) {
  type <- match.arg(type)
  Y <- methods::slot(object, slot.name)$similarity[[type]]$dr
  Groups <- methods::slot(object, slot.name)$similarity[[type]]$group
  prob <- methods::slot(object, slot.name)$prob
  if (is.null(pathway.remove)) {
    similarity <- methods::slot(object, slot.name)$similarity[[type]]$matrix
    pathway.remove <- rownames(similarity)[which(colSums(similarity) == 1)]
  }

  if (length(pathway.remove) > 0) {
    pathway.remove.idx <- which(dimnames(prob)[[3]] %in% pathway.remove)
    prob <- prob[ , , -pathway.remove.idx]
  }

  prob_sum <- apply(prob, 3, sum)
  df <- data.frame(x = Y[,1], y = Y[, 2], Commun.Prob. = prob_sum/max(prob_sum), labels = as.character(unlist(dimnames(prob)[3])), Groups = as.factor(Groups))
  if (is.null(color.use)) {
    color.use <- ggPalette(length(unique(Groups)))
  }
  gg <- ggplot(data = df, aes(x, y)) +
    geom_point(aes(size = Commun.Prob.,fill = Groups, colour = Groups), shape = 21) +
    CellChat_theme_opts() +
    theme(text = element_text(size = font.size), legend.key.height = grid::unit(0.15, "in"))+
    guides(colour = guide_legend(override.aes = list(size = 3)))+
    labs(title = title, x = xlabel, y = ylabel) + theme(plot.title = element_text(size= font.size.title, face="plain"))+
    scale_size_continuous(limits = c(0,1), range = dot.size, breaks = c(0.1,0.5,0.9)) +
    theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank()) +
    theme(axis.line.x = element_line(size = 0.25), axis.line.y = element_line(size = 0.25))
  gg <- gg + scale_fill_manual(values = ggplot2::alpha(color.use, alpha = dot.alpha), drop = FALSE)
  gg <- gg + scale_colour_manual(values = color.use, drop = FALSE)
  if (do.label) {
    gg <- gg + ggrepel::geom_text_repel(mapping = aes(label = labels, colour = Groups), size = label.size, show.legend = F,segment.size = 0.2, segment.alpha = 0.5)
  }

  if (length(pathway.remove) > 0 & pathway.remove.show) {
    gg <- gg + annotate(geom = 'text', label =  paste("Isolate pathways: ", paste(pathway.remove, collapse = ', ')), x = -Inf, y = Inf, hjust = 0, vjust = 1, size = label.size,fontface="italic")
  }
  if (!show.legend) {
    gg <- gg + theme(legend.position = "none")
  }

  if (!show.axes) {
    gg <- gg + theme_void()
  }
  gg
}


#' Zoom into the 2D visualization of the learned manifold learning of the signaling networks
#'
#' @param object CellChat object
#' @param slot.name the slot name of object that is used to compute centrality measures of signaling networks
#' @param type "functional","structural"
#' @param pathway.remove a character vector defining the signaling to remove
#' @param color.use defining the color for each cell group
#' @param nCol the number of columns of the plot
#' @param dot.size a range defining the size of the symbol
#' @param dot.alpha transparency
#' @param xlabel label of x-axis
#' @param ylabel label of y-axis
#' @param label.size font size of the text
#' @param do.label label the each point
#' @param show.legend whether show the legend
#' @param show.axes whether show the axes
#' @import ggplot2
#' @importFrom ggrepel geom_text_repel
#' @importFrom cowplot plot_grid
#' @importFrom methods slot
#' @return
#' @export
#'
#' @examples
netVisual_embeddingZoomIn <- function(object, slot.name = "netP", type = c("functional","structural"), color.use = NULL, pathway.remove = NULL,  nCol = 1, dot.size = c(2, 6), label.size = 2.8, dot.alpha = 0.5,
                                      xlabel = NULL, ylabel = NULL, do.label = T, show.legend = F, show.axes = T) {
  Y <- methods::slot(object, slot.name)$similarity[[type]]$dr
  clusters <- methods::slot(object, slot.name)$similarity[[type]]$group
  prob <- methods::slot(object, slot.name)$prob
  if (is.null(pathway.remove)) {
    similarity <- methods::slot(object, slot.name)$similarity[[type]]$matrix
    pathway.remove <- rownames(similarity)[which(colSums(similarity) == 1)]
  }

  if (length(pathway.remove) > 0) {
    pathway.remove.idx <- which(dimnames(prob)[[3]] %in% pathway.remove)
    prob <- prob[ , , -pathway.remove.idx]
  }

  prob_sum <- apply(prob, 3, sum)
  df <- data.frame(x = Y[,1], y = Y[, 2], Commun.Prob. = prob_sum/max(prob_sum), labels = as.character(unlist(dimnames(prob)[3])), clusters = as.factor(clusters))

  if (is.null(color.use)) {
    color.use <- ggPalette(length(unique(clusters)))
  }

  # zoom into each cluster and do labels
  ggAll <- vector("list", length(unique(clusters)))
  for (i in 1:length(unique(clusters))) {
    clusterID = i
    title <- paste0("Group ",  clusterID)
    df2 <- df[df$clusters %in% clusterID,]
    gg <- ggplot(data = df2, aes(x, y)) +
      geom_point(aes(size = Commun.Prob.), shape = 21, colour = alpha(color.use[clusterID], alpha = 1), fill = alpha(color.use[clusterID], alpha = dot.alpha)) +
      CellChat_theme_opts() +
      theme(text = element_text(size = 10), legend.key.height = grid::unit(0.15, "in"))+
      labs(title = title, x = xlabel, y = ylabel) + theme(plot.title = element_text(size=12))+
      scale_size_continuous(limits = c(0,1), range = dot.size, breaks = c(0.1,0.5,0.9)) +
      theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank()) +
      theme(axis.line.x = element_line(size = 0.25), axis.line.y = element_line(size = 0.25))
    if (do.label) {
      gg <- gg + ggrepel::geom_text_repel(mapping = aes(label = labels), colour = color.use[clusterID], size = label.size, segment.size = 0.2, segment.alpha = 0.5)
    }

    if (!show.legend) {
      gg <- gg + theme(legend.position = "none")
    }

    if (!show.axes) {
      gg <- gg + theme_void()
    }
    ggAll[[i]] <- gg
  }
  gg.combined <- cowplot::plot_grid(plotlist = ggAll, ncol = nCol)

  gg.combined

}



#' 2D visualization of the joint manifold learning of signaling networks from two datasets
#'
#' @param object CellChat object
#' @param slot.name the slot name of object that is used to compute centrality measures of signaling networks
#' @param type "functional","structural"
#' @param pathway.remove a character vector defining the signaling to remove
#' @param pathway.remove.show whether show the removed signaling names
#' @param color.use defining the color for each cell group
#' @param point.shape a numeric vector giving the point shapes. By default point.shape <- c(21, 0, 24, 23, 25, 10, 12), see available shapes at http://www.sthda.com/english/wiki/r-plot-pch-symbols-the-different-point-shapes-available-in-r
#' @param dot.size a range defining the size of the symbol
#' @param dot.alpha transparency
#' @param xlabel label of x-axis
#' @param ylabel label of y-axis
#' @param title main title of the plot
#' @param label.size font size of the text
#' @param do.label label the each point
#' @param show.legend whether show the legend
#' @param show.axes whether show the axes
#' @import ggplot2
#' @importFrom ggrepel geom_text_repel
#' @importFrom methods slot
#' @return
#' @export
#'
#' @examples
netVisual_embeddingPairwise <- function(object, slot.name = "netP", type = c("functional","structural"), color.use = NULL, point.shape = NULL, pathway.remove = NULL, pathway.remove.show = TRUE, dot.size = c(2, 6), label.size = 2.5, dot.alpha = 0.5,
                                        xlabel = "Dim 1", ylabel = "Dim 2", title = NULL,do.label = T, show.legend = T, show.axes = T) {
  type <- match.arg(type)
  Y <- methods::slot(object, slot.name)$similarity[[type]]$dr
  clusters <- methods::slot(object, slot.name)$similarity[[type]]$group
  object.names <- setdiff(names(methods::slot(object, slot.name)), "similarity")
  prob <- list()
  for (i in 1:length(setdiff(names(methods::slot(object, slot.name)), "similarity"))) {
    object.net <- methods::slot(object, slot.name)[[i]]
    prob[[i]] = object.net$prob
  }

  if (is.null(point.shape)) {
    point.shape <- c(21, 0, 24, 23, 25, 10, 12)
  }

  if (is.null(pathway.remove)) {
    similarity <- methods::slot(object, slot.name)$similarity[[type]]$matrix
    pathway.remove <- rownames(similarity)[which(colSums(similarity) == 1)]
  }

  if (length(pathway.remove) > 0) {
    for (i in 1:length(prob)) {
      probi <- prob[[i]]
      pathway.remove.idx <- which(dimnames(probi)[[3]] %in% pathway.remove)
      if (length(pathway.remove.idx) > 0) {
        probi <- probi[ , , -pathway.remove.idx]
      }
      prob[[i]] <- probi
    }
  }
  prob_sum.each <- list()
  signalingAll <- c()
  for (i in 1:length(prob)) {
    probi <- prob[[i]]
    prob_sum.each[[i]] <- apply(probi, 3, sum)
    signalingAll <- c(signalingAll, paste0(names(prob_sum.each[[i]]),"-",object.names[i]))
  }
  prob_sum <- unlist(prob_sum.each)
  names(prob_sum) <- signalingAll

  group <- sub(".*-", "", names(prob_sum))
  labels = sub("-.*", "", names(prob_sum))

  df <- data.frame(x = Y[,1], y = Y[, 2], Commun.Prob. = prob_sum/max(prob_sum),
                   labels = as.character(labels), clusters = as.factor(clusters), group = factor(group, levels = unique(group)))
  # color dots (light inside color and dark border) based on clustering and no labels
  if (is.null(color.use)) {
    color.use <- ggPalette(length(unique(clusters)))
  }
  gg <- ggplot(data = df, aes(x, y)) +
    geom_point(aes(size = Commun.Prob.,fill = clusters, colour = clusters, shape = group)) +
    CellChat_theme_opts() +
    theme(text = element_text(size = 10), legend.key.height = grid::unit(0.15, "in"))+
    guides(colour = guide_legend(override.aes = list(size = 3)))+
    labs(title = title, x = xlabel, y = ylabel) +
    scale_size_continuous(limits = c(0,1), range = dot.size, breaks = c(0.1,0.5,0.9)) +
    theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank()) +
    theme(axis.line.x = element_line(size = 0.25), axis.line.y = element_line(size = 0.25))
  gg <- gg + scale_fill_manual(values = ggplot2::alpha(color.use, alpha = dot.alpha), drop = FALSE) #+ scale_alpha(group, range = c(0.1, 1))
  gg <- gg + scale_colour_manual(values = color.use, drop = FALSE)
  gg <- gg + scale_shape_manual(values = point.shape[1:length(prob)])
  if (do.label) {
      gg <- gg + ggrepel::geom_text_repel(mapping = aes(label = labels, colour = clusters, alpha=group), size = label.size, show.legend = F,segment.size = 0.2, segment.alpha = 0.5) + scale_alpha_discrete(range = c(1, 0.6))
  }

  if (length(pathway.remove) > 0 & pathway.remove.show) {
    gg <- gg + annotate(geom = 'text', label =  paste("Isolate pathways: ", paste(pathway.remove, collapse = ', ')), x = -Inf, y = Inf, hjust = 0, vjust = 1, size = label.size,fontface="italic")
  }

  if (!show.legend) {
    gg <- gg + theme(legend.position = "none")
  }

  if (!show.axes) {
    gg <- gg + theme_void()
  }
  gg
}



#' Zoom into the 2D visualization of the joint manifold learning of signaling networks from two datasets
#'
#' @param object CellChat object
#' @param slot.name the slot name of object that is used to compute centrality measures of signaling networks
#' @param type "functional","structural"
#' @param pathway.remove a character vector defining the signaling to remove
#' @param color.use defining the color for each cell group
#' @param nCol number of columns in the plot
#' @param point.shape a numeric vector giving the point shapes. By default point.shape <- c(21, 0, 24, 23, 25, 10, 12), see available shapes at http://www.sthda.com/english/wiki/r-plot-pch-symbols-the-different-point-shapes-available-in-r
#' @param dot.size a range defining the size of the symbol
#' @param dot.alpha transparency
#' @param xlabel label of x-axis
#' @param ylabel label of y-axis
#' @param label.size font size of the text
#' @param do.label label the each point
#' @param show.legend whether show the legend
#' @param show.axes whether show the axes
#' @import ggplot2
#' @importFrom ggrepel geom_text_repel
#' @importFrom methods slot
#' @return
#' @export
#'
#' @examples
netVisual_embeddingPairwiseZoomIn <- function(object, slot.name = "netP", type = c("functional","structural"), color.use = NULL, nCol = 1, point.shape = NULL, pathway.remove = NULL, dot.size = c(2, 6), label.size = 2.8, dot.alpha = 0.5,
                                              xlabel = NULL, ylabel = NULL, do.label = T, show.legend = F, show.axes = T) {

  type <- match.arg(type)
  Y <- methods::slot(object, slot.name)$similarity[[type]]$dr
  clusters <- methods::slot(object, slot.name)$similarity[[type]]$group

  object.names <- setdiff(names(methods::slot(object, slot.name)), "similarity")
  prob <- list()
  for (i in 1:length(setdiff(names(methods::slot(object, slot.name)), "similarity"))) {
    object.net <- methods::slot(object, slot.name)[[i]]
    prob[[i]] = object.net$prob
  }

  if (is.null(point.shape)) {
    point.shape <- c(21, 0, 24, 23, 25, 10, 12)
  }

  if (is.null(pathway.remove)) {
    similarity <- methods::slot(object, slot.name)$similarity[[type]]$matrix
    pathway.remove <- rownames(similarity)[which(colSums(similarity) == 1)]
  }

  if (length(pathway.remove) > 0) {
    for (i in 1:length(prob)) {
      probi <- prob[[i]]
      pathway.remove.idx <- which(dimnames(probi)[[3]] %in% pathway.remove)
      if (length(pathway.remove.idx) > 0) {
        probi <- probi[ , , -pathway.remove.idx]
      }
      prob[[i]] <- probi
    }
  }

  prob_sum.each <- list()
  signalingAll <- c()
  for (i in 1:length(prob)) {
    probi <- prob[[i]]
    prob_sum.each[[i]] <- apply(probi, 3, sum)
    signalingAll <- c(signalingAll, paste0(names(prob_sum.each[[i]]),"-",object.names[i]))
  }
  prob_sum <- unlist(prob_sum.each)
  names(prob_sum) <- signalingAll

  group <- sub(".*-", "", names(prob_sum))
  labels = sub("-.*", "", names(prob_sum))

  df <- data.frame(x = Y[,1], y = Y[, 2], Commun.Prob. = prob_sum/max(prob_sum),
                   labels = as.character(labels), clusters = as.factor(clusters), group = factor(group, levels = unique(group)))
  if (is.null(color.use)) {
    color.use <- ggPalette(length(unique(clusters)))
  }

  # zoom into each cluster and do labels
  ggAll <- vector("list", length(unique(clusters)))
  for (i in 1:length(unique(clusters))) {
    clusterID = i
    title <- paste0("Cluster ",  clusterID)
    df2 <- df[df$clusters %in% clusterID,]
    gg <- ggplot(data = df2, aes(x, y)) +
      geom_point(aes(size = Commun.Prob., shape = group),fill = alpha(color.use[clusterID], alpha = dot.alpha), colour = alpha(color.use[clusterID], alpha = 1)) +
      CellChat_theme_opts() +
      theme(text = element_text(size = 10), legend.key.height = grid::unit(0.15, "in"))+
      guides(colour = guide_legend(override.aes = list(size = 3)))+
      labs(title = title, x = xlabel, y = ylabel) +
      scale_size_continuous(limits = c(0,1), range = dot.size, breaks = c(0.1,0.5,0.9)) +
      theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank()) +
      theme(axis.line.x = element_line(size = 0.25), axis.line.y = element_line(size = 0.25))
    idx <- match(unique(df2$group), levels(df$group), nomatch = 0)
    gg <- gg + scale_shape_manual(values= point.shape[idx])
    if (do.label) {
        gg <- gg + ggrepel::geom_text_repel(mapping = aes(label = labels), colour = color.use[clusterID], size = label.size, show.legend = F,segment.size = 0.2, segment.alpha = 0.5) + scale_alpha_discrete(range = c(1, 0.6))
    }

    if (!show.legend) {
      gg <- gg + theme(legend.position = "none")
    }

    if (!show.axes) {
      gg <- gg + theme_void()
    }
    ggAll[[i]] <- gg
  }
  gg.combined <- cowplot::plot_grid(plotlist = ggAll, ncol = nCol)

  gg.combined

}



#' Show the description of CellChatDB databse
#'
#' @param CellChatDB CellChatDB databse
#' @param nrow the number of rows in the plot
#' @importFrom dplyr group_by summarise
#'
#' @return
#' @export
#'
showDatabaseCategory <- function(CellChatDB, nrow = 1) {
  interaction_input <- CellChatDB$interaction
  geneIfo <- CellChatDB$geneInfo
  df <- interaction_input %>% group_by(annotation) %>% summarise(value=n())
  df$group <- factor(df$annotation, levels = c("Secreted Signaling","ECM-Receptor","Cell-Cell Contact"))
  gg1 <- pieChart(df)
  binary <- (interaction_input$ligand %in% geneIfo$Symbol) & (interaction_input$receptor %in% geneIfo$Symbol)
  df <- data.frame(group = rep("Heterodimers", dim(interaction_input)[1]),stringsAsFactors = FALSE)
  df$group[binary] <- rep("Others",sum(binary),1)
  df <- df %>% group_by(group) %>% summarise(value=n())
  df$group <- factor(df$group, levels = c("Heterodimers","Others"))
  gg2 <- pieChart(df)

  kegg <- grepl("KEGG", interaction_input$evidence)
  df <- data.frame(group = rep("Literature", dim(interaction_input)[1]),stringsAsFactors = FALSE)
  df$group[kegg] <- rep("KEGG",sum(kegg),1)
  df <- df %>% group_by(group) %>% summarise(value=n())
  df$group <- factor(df$group, levels = c("KEGG","Literature"))
  gg3 <- pieChart(df)

  gg <- cowplot::plot_grid(gg1, gg2, gg3, nrow = nrow, align = "h", rel_widths = c(1, 1,1))
  return(gg)
}


#' Plot pie chart
#'
#' @param df a dataframe
#' @param label.size a character
#' @param color.use the name of the variable in CellChatDB interaction_input
#' @param title the title of plot
#' @import ggplot2
#' @importFrom scales percent
#' @importFrom dplyr arrange desc mutate
#' @importFrom ggrepel geom_text_repel
#' @return
#' @export
#'
pieChart <- function(df, label.size = 2.5, color.use = NULL, title = "") {
  df %>% arrange(dplyr::desc(value)) %>%
    mutate(prop = scales::percent(value/sum(value))) -> df

  gg <- ggplot(df, aes(x="", y=value, fill=forcats::fct_inorder(group))) +
    geom_bar(stat="identity", width=1) +
    coord_polar("y", start=0)+theme_void() +
    ggrepel::geom_text_repel(aes(label = prop), size= label.size, show.legend = F, nudge_x = 0)
  gg <- gg + theme(legend.position="bottom", legend.direction = "vertical")

  if(!is.null(color.use)) {
    gg <- gg + scale_color_manual(color.use)
  }

  if (!is.null(title)) {
    gg <- gg + guides(fill = guide_legend(title = title))
  }
  gg
}
