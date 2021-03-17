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
#' Automatically save plots in the current working directory.
#'
#' @param object CellChat object
#' @param signaling a signaling pathway name
#' @param signaling.name alternative signaling pathway name to show on the plot
#' @param color.use the character vector defining the color of each cell group
#' @param vertex.receiver a numeric vector giving the index of the cell groups as targets in the first hierarchy plot
#' @param top the fraction of interactions to show (0 < top <= 1)
#' @param sources.use a vector giving the index or the name of source cell groups
#' @param targets.use a vector giving the index or the name of target cell groups.
#' @param remove.isolate whether remove the isolate nodes in the communication network
#' @param weight.scale whether scale the edge weight
#' @param vertex.weight The weight of vertex: either a scale value or a vector
#' @param vertex.weight.max the maximum weight of vertex; defualt = max(vertex.weight)
#' @param vertex.size.max the maximum vertex size for visualization
#' @param edge.weight.max.individual the maximum weight of edge when plotting the individual L-R netwrok; defualt = max(net)
#' @param edge.weight.max.aggregate the maximum weight of edge when plotting the aggregated signaling pathway network
#' @param edge.width.max The maximum edge width for visualization
#' @param layout "hierarchy", "circle" or "chord"
#' @param height height of plot
#' @param thresh threshold of the p-value for determining significant interaction
#' @param pt.title font size of the text
#' @param title.space the space between the title and plot
#' @param vertex.label.cex The label size of vertex in the network
#' @param out.format the format of output figures: svg, png and pdf
#' @param from,to,bidirection Deprecated. Use `sources.use`,`targets.use`
#' @param vertex.size Deprecated. Use `vertex.weight`
#'
#' Parameters below are set for "chord" diagram. Please also check the function `netVisual_chord_cell` for more parameters.
#' @param group A named group labels for making multiple-group Chord diagrams. The sector names should be used as the names in the vector.
#' The order of group controls the sector orders and if group is set as a factor, the order of levels controls the order of groups.
#' @param cell.order a char vector defining the cell type orders (sector orders)
#' @param small.gap Small gap between sectors.
#' @param big.gap Gap between the different sets of sectors, which are defined in the `group` parameter
#' @param scale scale each sector to same width; default = FALSE; however, it is set to be TRUE when remove.isolate = TRUE
#' @param reduce if the ratio of the width of certain grid compared to the whole circle is less than this value, the grid is removed on the plot. Set it to value less than zero if you want to keep all tiny grid.
#' @param show.legend whether show the figure legend
#' @param legend.pos.x,legend.pos.y adjust the legend position
#' @param nCol number of columns when displaying the network mediated by ligand-receptor using "circle" or "chord"
#'
#' @param ... other parameters (e.g.,vertex.label.cex, vertex.label.color, alpha.edge, label.edge, edge.label.color, edge.label.cex, edge.curved)
#'  passing to `netVisual_hierarchy1`,`netVisual_hierarchy2`,`netVisual_circle`. NB: some parameters might be not supported
#' @importFrom svglite svglite
#' @importFrom grDevices dev.off pdf
#'
#' @return
#' @export
#'
#' @examples
#'
netVisual <- function(object, signaling, signaling.name = NULL, color.use = NULL, vertex.receiver = NULL, sources.use = NULL, targets.use = NULL, top = 1, remove.isolate = FALSE,
                      vertex.weight = NULL, vertex.weight.max = NULL, vertex.size.max = 15,
                      weight.scale = TRUE, edge.weight.max.individual = NULL, edge.weight.max.aggregate = NULL, edge.width.max=8,
                      layout = c("hierarchy","circle","chord"), height = 5, thresh = 0.05, pt.title = 12, title.space = 6, vertex.label.cex = 0.8,from = NULL, to = NULL, bidirection = NULL,vertex.size = NULL,
                      out.format = c("svg","png"),
                      group = NULL,cell.order = NULL,small.gap = 1, big.gap = 10, scale = FALSE, reduce = -1, show.legend = FALSE, legend.pos.x = 20,legend.pos.y = 20, nCol = NULL,
                      ...) {
  layout <- match.arg(layout)
  if (!is.null(vertex.size)) {
    warning("'vertex.size' is deprecated. Use `vertex.weight`")
  }
  if (is.null(vertex.weight)) {
    vertex.weight <- as.numeric(table(object@idents))
  }
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

  if (is.null(nCol)) {
    nCol <- min(length(pairLR.name.use), 2)
  }

  if (length(dim(prob)) == 2) {
    prob <- replicate(1, prob, simplify="array")
    pval <- replicate(1, pval, simplify="array")
  }
#  prob <-(prob-min(prob))/(max(prob)-min(prob))
  if (is.null(edge.weight.max.individual)) {
    edge.weight.max.individual = max(prob)
  }
  prob.sum <- apply(prob, c(1,2), sum)
  #  prob.sum <-(prob.sum-min(prob.sum))/(max(prob.sum)-min(prob.sum))
  if (is.null(edge.weight.max.aggregate)) {
    edge.weight.max.aggregate = max(prob.sum)
  }

  if (layout == "hierarchy") {
    if (is.element("svg", out.format)) {
      svglite::svglite(file = paste0(signaling.name, "_hierarchy_individual.svg"), width = 8, height = nRow*height)
      par(mfrow=c(nRow,2), mar = c(5, 4, 4, 2) +0.1)
      for (i in 1:length(pairLR.name.use)) {
        #signalName_i <- paste0(pairLR$ligand[i], "-",pairLR$receptor[i], sep = "")
        signalName_i <- pairLR$interaction_name_2[i]
        prob.i <- prob[,,i]
        netVisual_hierarchy1(prob.i, vertex.receiver = vertex.receiver, sources.use = sources.use, targets.use = targets.use, remove.isolate = remove.isolate, top = top, color.use = color.use, vertex.weight = vertex.weight, vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max, weight.scale = weight.scale, edge.weight.max = edge.weight.max.individual, edge.width.max=edge.width.max, title.name = signalName_i, vertex.label.cex = vertex.label.cex,...)
        netVisual_hierarchy2(prob.i, vertex.receiver = setdiff(1:nrow(prob.i),vertex.receiver), sources.use = sources.use, targets.use = targets.use, remove.isolate = remove.isolate, top = top, color.use = color.use, vertex.weight = vertex.weight, vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max, weight.scale = weight.scale, edge.weight.max = edge.weight.max.individual, edge.width.max=edge.width.max, title.name = signalName_i, vertex.label.cex = vertex.label.cex,...)
      }
      dev.off()
    }
    if (is.element("png", out.format)) {
      grDevices::png(paste0(signaling.name, "_hierarchy_individual.png"), width = 8, height = nRow*height, units = "in",res = 300)
      par(mfrow=c(nRow,2), mar = c(5, 4, 4, 2) +0.1)
      for (i in 1:length(pairLR.name.use)) {
        signalName_i <- pairLR$interaction_name_2[i]
        prob.i <- prob[,,i]
        netVisual_hierarchy1(prob.i, vertex.receiver = vertex.receiver, sources.use = sources.use, targets.use = targets.use, remove.isolate = remove.isolate, top = top, color.use = color.use, vertex.weight = vertex.weight, vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max, weight.scale = weight.scale, edge.weight.max = edge.weight.max.individual, edge.width.max=edge.width.max, title.name = signalName_i, vertex.label.cex = vertex.label.cex,...)
        netVisual_hierarchy2(prob.i, vertex.receiver = setdiff(1:nrow(prob.i),vertex.receiver), sources.use = sources.use, targets.use = targets.use, remove.isolate = remove.isolate, top = top, color.use = color.use, vertex.weight = vertex.weight, vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max, weight.scale = weight.scale, edge.weight.max =edge.weight.max.individual, edge.width.max=edge.width.max, title.name = signalName_i, vertex.label.cex = vertex.label.cex,...)
      }
      dev.off()
    }
    if (is.element("pdf", out.format)) {
     # grDevices::pdf(paste0(signaling.name, "_hierarchy_individual.pdf"), width = 8, height = nRow*height)
      grDevices::cairo_pdf(paste0(signaling.name, "_hierarchy_individual.pdf"), width = 8, height = nRow*height)
      par(mfrow=c(nRow,2), mar = c(5, 4, 4, 2) +0.1)
      for (i in 1:length(pairLR.name.use)) {
        signalName_i <- pairLR$interaction_name_2[i]
        prob.i <- prob[,,i]
        netVisual_hierarchy1(prob.i, vertex.receiver = vertex.receiver, sources.use = sources.use, targets.use = targets.use, remove.isolate = remove.isolate, top = top, color.use = color.use, vertex.weight = vertex.weight, vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max, weight.scale = weight.scale, edge.weight.max = edge.weight.max.individual, edge.width.max=edge.width.max, title.name = signalName_i, vertex.label.cex = vertex.label.cex,...)
        netVisual_hierarchy2(prob.i, vertex.receiver = setdiff(1:nrow(prob.i),vertex.receiver), sources.use = sources.use, targets.use = targets.use, remove.isolate = remove.isolate, top = top, color.use = color.use, vertex.weight = vertex.weight, vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max, weight.scale = weight.scale, edge.weight.max =edge.weight.max.individual, edge.width.max=edge.width.max, title.name = signalName_i, vertex.label.cex = vertex.label.cex,...)
      }
      dev.off()
    }


    if (is.element("svg", out.format)) {
      svglite::svglite(file = paste0(signaling.name, "_hierarchy_aggregate.svg"), width = 7, height = 1*height)
      par(mfrow=c(1,2), ps = pt.title)
      netVisual_hierarchy1(prob.sum, vertex.receiver = vertex.receiver, sources.use = sources.use, targets.use = targets.use, remove.isolate = remove.isolate, top = top, color.use = color.use, vertex.weight = vertex.weight, vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max, weight.scale = weight.scale, edge.weight.max = edge.weight.max.aggregate, edge.width.max=edge.width.max,title.name = NULL, vertex.label.cex = vertex.label.cex,...)
      netVisual_hierarchy2(prob.sum, vertex.receiver = setdiff(1:nrow(prob.sum),vertex.receiver), sources.use = sources.use, targets.use = targets.use, remove.isolate = remove.isolate, top = top, color.use = color.use, vertex.weight = vertex.weight, vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max, weight.scale = weight.scale, edge.weight.max = edge.weight.max.aggregate, edge.width.max=edge.width.max,title.name = NULL, vertex.label.cex = vertex.label.cex,...)
      graphics::mtext(paste0(signaling.name, " signaling pathway network"), side = 3, outer = TRUE, cex = 1, line = -title.space)
      dev.off()
    }
    if (is.element("png", out.format)) {
      grDevices::png(paste0(signaling.name, "_hierarchy_aggregate.png"), width = 7, height = 1*height, units = "in",res = 300)
      par(mfrow=c(1,2), ps = pt.title)
      netVisual_hierarchy1(prob.sum, vertex.receiver = vertex.receiver, sources.use = sources.use, targets.use = targets.use, remove.isolate = remove.isolate, top = top, color.use = color.use, vertex.weight = vertex.weight, vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max, weight.scale = weight.scale, edge.weight.max = edge.weight.max.aggregate, edge.width.max=edge.width.max, title.name = NULL, vertex.label.cex = vertex.label.cex,...)
      netVisual_hierarchy2(prob.sum, vertex.receiver = setdiff(1:nrow(prob.sum),vertex.receiver), sources.use = sources.use, targets.use = targets.use, remove.isolate = remove.isolate, top = top, color.use = color.use, vertex.weight = vertex.weight, vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max, weight.scale = weight.scale, edge.weight.max = edge.weight.max.aggregate, edge.width.max=edge.width.max,title.name = NULL, vertex.label.cex = vertex.label.cex,...)
      graphics::mtext(paste0(signaling.name, " signaling pathway network"), side = 3, outer = TRUE, cex = 1, line = -title.space)
      dev.off()
    }
    if (is.element("pdf", out.format)) {
      # grDevices::pdf(paste0(signaling.name, "_hierarchy_aggregate.pdf"), width = 7, height = 1*height)
      grDevices::cairo_pdf(paste0(signaling.name, "_hierarchy_aggregate.pdf"), width = 7, height = 1*height)
      par(mfrow=c(1,2), ps = pt.title)
      netVisual_hierarchy1(prob.sum, vertex.receiver = vertex.receiver, sources.use = sources.use, targets.use = targets.use, remove.isolate = remove.isolate, top = top, color.use = color.use, vertex.weight = vertex.weight, vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max, weight.scale = weight.scale, edge.weight.max = edge.weight.max.aggregate, edge.width.max=edge.width.max, title.name = NULL, vertex.label.cex = vertex.label.cex,...)
      netVisual_hierarchy2(prob.sum, vertex.receiver = setdiff(1:nrow(prob.sum),vertex.receiver), sources.use = sources.use, targets.use = targets.use, remove.isolate = remove.isolate, top = top, color.use = color.use, vertex.weight = vertex.weight, vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max, weight.scale = weight.scale, edge.weight.max = edge.weight.max.aggregate, edge.width.max=edge.width.max, title.name = NULL, vertex.label.cex = vertex.label.cex,...)
      graphics::mtext(paste0(signaling.name, " signaling pathway network"), side = 3, outer = TRUE, cex = 1, line = -title.space)
      dev.off()
    }

  } else if (layout == "circle") {
    if (is.element("svg", out.format)) {
      svglite::svglite(file = paste0(signaling.name,"_", layout, "_individual.svg"), width = height, height = nRow*height)
     # par(mfrow=c(nRow,1))
      par(mfrow = c(ceiling(length(pairLR.name.use)/nCol), nCol), xpd=TRUE)
      for (i in 1:length(pairLR.name.use)) {
        #signalName_i <- paste0(pairLR$ligand[i], "-",pairLR$receptor[i], sep = "")
        signalName_i <- pairLR$interaction_name_2[i]
        prob.i <- prob[,,i]
        netVisual_circle(prob.i, sources.use = sources.use, targets.use = targets.use, remove.isolate = remove.isolate, top = top, color.use = color.use, vertex.weight = vertex.weight, vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max, weight.scale = weight.scale, edge.weight.max = edge.weight.max.individual, edge.width.max=edge.width.max, title.name = signalName_i, vertex.label.cex = vertex.label.cex,...)
      }
      dev.off()
    }
    if (is.element("png", out.format)) {
      grDevices::png(paste0(signaling.name,"_", layout, "_individual.png"), width = height, height = nRow*height, units = "in",res = 300)
     # par(mfrow=c(nRow,1))
      par(mfrow = c(ceiling(length(pairLR.name.use)/nCol), nCol), xpd=TRUE)
      for (i in 1:length(pairLR.name.use)) {
        #signalName_i <- paste0(pairLR$ligand[i], "-",pairLR$receptor[i], sep = "")
        signalName_i <- pairLR$interaction_name_2[i]
        prob.i <- prob[,,i]
        netVisual_circle(prob.i, sources.use = sources.use, targets.use = targets.use, remove.isolate = remove.isolate, top = top, color.use = color.use, vertex.weight = vertex.weight, vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max, weight.scale = weight.scale, edge.weight.max = edge.weight.max.individual, edge.width.max=edge.width.max, title.name = signalName_i, vertex.label.cex = vertex.label.cex,...)
      }
      dev.off()
    }
    if (is.element("pdf", out.format)) {
      # grDevices::pdf(paste0(signaling.name,"_", layout, "_individual.pdf"), width = height, height = nRow*height)
      grDevices::cairo_pdf(paste0(signaling.name,"_", layout, "_individual.pdf"), width = height, height = nRow*height)
     # par(mfrow=c(nRow,1))
      par(mfrow = c(ceiling(length(pairLR.name.use)/nCol), nCol), xpd=TRUE)
      for (i in 1:length(pairLR.name.use)) {
        #signalName_i <- paste0(pairLR$ligand[i], "-",pairLR$receptor[i], sep = "")
        signalName_i <- pairLR$interaction_name_2[i]
        prob.i <- prob[,,i]
        netVisual_circle(prob.i, sources.use = sources.use, targets.use = targets.use, remove.isolate = remove.isolate, top = top, color.use = color.use, vertex.weight = vertex.weight, vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max, weight.scale = weight.scale, edge.weight.max = edge.weight.max.individual, edge.width.max=edge.width.max,title.name = signalName_i, vertex.label.cex = vertex.label.cex,...)
      }
      dev.off()
    }

  #  prob.sum <- apply(prob, c(1,2), sum)
  #  prob.sum <-(prob.sum-min(prob.sum))/(max(prob.sum)-min(prob.sum))
    if (is.element("svg", out.format)) {
      svglite(file = paste0(signaling.name,"_", layout,  "_aggregate.svg"), width = height, height = 1*height)
      netVisual_circle(prob.sum, sources.use = sources.use, targets.use = targets.use, remove.isolate = remove.isolate, top = top, color.use = color.use, vertex.weight = vertex.weight, vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max, weight.scale = weight.scale, edge.weight.max = edge.weight.max.aggregate, edge.width.max=edge.width.max,title.name = paste0(signaling.name, " signaling pathway network"), vertex.label.cex = vertex.label.cex,...)
      dev.off()
    }
    if (is.element("png", out.format)) {
      grDevices::png(paste0(signaling.name,"_", layout,  "_aggregate.png"), width = height, height = 1*height, units = "in",res = 300)
      netVisual_circle(prob.sum, sources.use = sources.use, targets.use = targets.use, remove.isolate = remove.isolate, top = top, color.use = color.use, vertex.weight = vertex.weight, vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max, weight.scale = weight.scale, edge.weight.max = edge.weight.max.aggregate, edge.width.max=edge.width.max,title.name = paste0(signaling.name, " signaling pathway network"), vertex.label.cex = vertex.label.cex,...)
      dev.off()
    }
    if (is.element("pdf", out.format)) {
     # grDevices::pdf(paste0(signaling.name,"_", layout,  "_aggregate.pdf"), width = height, height = 1*height)
      grDevices::cairo_pdf(paste0(signaling.name,"_", layout,  "_aggregate.pdf"), width = height, height = 1*height)
      netVisual_circle(prob.sum, sources.use = sources.use, targets.use = targets.use, remove.isolate = remove.isolate, top = top, color.use = color.use, vertex.weight = vertex.weight, vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max, weight.scale = weight.scale, edge.weight.max = edge.weight.max.aggregate, edge.width.max=edge.width.max, title.name = paste0(signaling.name, " signaling pathway network"), vertex.label.cex = vertex.label.cex,...)
      dev.off()
    }
  } else if (layout == "chord") {
    if (is.element("svg", out.format)) {

      svglite::svglite(file = paste0(signaling.name,"_", layout, "_individual.svg"), width = height, height = nRow*height)
      par(mfrow = c(ceiling(length(pairLR.name.use)/nCol), nCol), xpd=TRUE)
    #  gg <- vector("list", length(pairLR.name.use))
      for (i in 1:length(pairLR.name.use)) {
        title.name <- pairLR$interaction_name_2[i]
        net <- prob[,,i]
        netVisual_chord_cell_internal(net, color.use = color.use, sources.use = sources.use, targets.use = targets.use, remove.isolate = remove.isolate,
                                                 group = group, cell.order = cell.order,
                                                 lab.cex = vertex.label.cex,small.gap = small.gap, big.gap = big.gap,
                                                 scale = scale, reduce = reduce,
                                                 title.name = title.name, show.legend = show.legend, legend.pos.x = legend.pos.x,legend.pos.y=legend.pos.y)
      }
      dev.off()
    }
    if (is.element("png", out.format)) {
      grDevices::png(paste0(signaling.name,"_", layout, "_individual.png"), width = height, height = nRow*height, units = "in",res = 300)
      par(mfrow = c(ceiling(length(pairLR.name.use)/nCol), nCol), xpd=TRUE)
      #  gg <- vector("list", length(pairLR.name.use))
      for (i in 1:length(pairLR.name.use)) {
        title.name <- pairLR$interaction_name_2[i]
        net <- prob[,,i]
        netVisual_chord_cell_internal(net, color.use = color.use, sources.use = sources.use, targets.use = targets.use, remove.isolate = remove.isolate,
                                      group = group, cell.order = cell.order,
                                      lab.cex = vertex.label.cex,small.gap = small.gap, big.gap = big.gap,
                                      scale = scale, reduce = reduce,
                                      title.name = title.name, show.legend = show.legend, legend.pos.x = legend.pos.x,legend.pos.y=legend.pos.y)
      }
      dev.off()
    }
    if (is.element("pdf", out.format)) {
      # grDevices::pdf(paste0(signaling.name,"_", layout, "_individual.pdf"), width = height, height = nRow*height)
      grDevices::cairo_pdf(paste0(signaling.name,"_", layout, "_individual.pdf"), width = height, height = nRow*height)
      par(mfrow = c(ceiling(length(pairLR.name.use)/nCol), nCol), xpd=TRUE)
      #  gg <- vector("list", length(pairLR.name.use))
      for (i in 1:length(pairLR.name.use)) {
        title.name <- pairLR$interaction_name_2[i]
        net <- prob[,,i]
        netVisual_chord_cell_internal(net, color.use = color.use, sources.use = sources.use, targets.use = targets.use, remove.isolate = remove.isolate,
                                      group = group, cell.order = cell.order,
                                      lab.cex = vertex.label.cex,small.gap = small.gap, big.gap = big.gap,
                                      scale = scale, reduce = reduce,
                                      title.name = title.name, show.legend = show.legend, legend.pos.x = legend.pos.x,legend.pos.y=legend.pos.y)
      }
      dev.off()
    }

  #  prob.sum <- apply(prob, c(1,2), sum)
    if (is.element("svg", out.format)) {
      svglite(file = paste0(signaling.name,"_", layout,  "_aggregate.svg"), width = height, height = 1*height)
      netVisual_chord_cell_internal(prob.sum, color.use = color.use, sources.use = sources.use, targets.use = targets.use, remove.isolate = remove.isolate,
                                    group = group, cell.order = cell.order,
                                    lab.cex = vertex.label.cex,small.gap = small.gap, big.gap = big.gap,
                                    scale = scale, reduce = reduce,
                                    title.name = paste0(signaling.name, " signaling pathway network"), show.legend = show.legend, legend.pos.x = legend.pos.x,legend.pos.y=legend.pos.y)
      dev.off()
    }
    if (is.element("png", out.format)) {
      grDevices::png(paste0(signaling.name,"_", layout,  "_aggregate.png"), width = height, height = 1*height, units = "in",res = 300)
      netVisual_chord_cell_internal(prob.sum, color.use = color.use, sources.use = sources.use, targets.use = targets.use, remove.isolate = remove.isolate,
                                    group = group, cell.order = cell.order,
                                    lab.cex = vertex.label.cex,small.gap = small.gap, big.gap = big.gap,
                                    scale = scale, reduce = reduce,
                                    title.name = paste0(signaling.name, " signaling pathway network"), show.legend = show.legend, legend.pos.x = legend.pos.x,legend.pos.y=legend.pos.y)
      dev.off()
    }
    if (is.element("pdf", out.format)) {
      # grDevices::pdf(paste0(signaling.name,"_", layout,  "_aggregate.pdf"), width = height, height = 1*height)
      grDevices::cairo_pdf(paste0(signaling.name,"_", layout,  "_aggregate.pdf"), width = height, height = 1*height)
      netVisual_chord_cell_internal(prob.sum, color.use = color.use, sources.use = sources.use, targets.use = targets.use, remove.isolate = remove.isolate,
                                    group = group, cell.order = cell.order,
                                    lab.cex = vertex.label.cex,small.gap = small.gap, big.gap = big.gap,
                                    scale = scale, reduce = reduce,
                                    title.name = paste0(signaling.name, " signaling pathway network"), show.legend = show.legend, legend.pos.x = legend.pos.x,legend.pos.y=legend.pos.y)
      dev.off()
    }
  }

}


#' Visualize the inferred signaling network of signaling pathways by aggregating all L-R pairs
#'
#' @param object CellChat object
#' @param signaling a signaling pathway name
#' @param signaling.name alternative signaling pathway name to show on the plot
#' @param color.use the character vector defining the color of each cell group
#' @param vertex.receiver a numeric vector giving the index of the cell groups as targets in the first hierarchy plot
#' @param sources.use a vector giving the index or the name of source cell groups
#' @param targets.use a vector giving the index or the name of target cell groups.
#' @param remove.isolate whether remove the isolate nodes in the communication network
#' @param top the fraction of interactions to show
#' @param weight.scale whether scale the edge weight
#' @param vertex.weight The weight of vertex: either a scale value or a vector
#' @param vertex.weight.max the maximum weight of vertex; defualt = max(vertex.weight)
#' @param vertex.size.max the maximum vertex size for visualization
#' @param edge.weight.max the maximum weight of edge; defualt = max(net)
#' @param edge.width.max The maximum edge width for visualization
#' @param layout "hierarchy", "circle" or "chord"
#' @param thresh threshold of the p-value for determining significant interaction
#' @param pt.title font size of the text
#' @param title.space the space between the title and plot
#' @param vertex.label.cex The label size of vertex in the network
#' @param from,to,bidirection Deprecated. Use `sources.use`,`targets.use`
#' @param vertex.size Deprecated. Use `vertex.weight`
#'
#' Parameters below are set for "chord" diagram. Please also check the function `netVisual_chord_cell` for more parameters.
#' @param group A named group labels for making multiple-group Chord diagrams. The sector names should be used as the names in the vector.
#' The order of group controls the sector orders and if group is set as a factor, the order of levels controls the order of groups.
#' @param cell.order a char vector defining the cell type orders (sector orders)
#' @param small.gap Small gap between sectors.
#' @param big.gap Gap between the different sets of sectors, which are defined in the `group` parameter
#' @param scale scale each sector to same width; default = FALSE; however, it is set to be TRUE when remove.isolate = TRUE
#' @param reduce if the ratio of the width of certain grid compared to the whole circle is less than this value, the grid is removed on the plot. Set it to value less than zero if you want to keep all tiny grid.
#' @param show.legend whether show the figure legend
#' @param legend.pos.x,legend.pos.y adjust the legend position
#'
#' @param ... other parameters (e.g.,vertex.label.cex, vertex.label.color, alpha.edge, label.edge, edge.label.color, edge.label.cex, edge.curved)
#'  passing to `netVisual_hierarchy1`,`netVisual_hierarchy2`,`netVisual_circle`. NB: some parameters might be not supported
#' @importFrom grDevices recordPlot
#'
#' @return  an object of class "recordedplot"
#' @export
#'
#'
netVisual_aggregate <- function(object, signaling, signaling.name = NULL, color.use = NULL, vertex.receiver = NULL, sources.use = NULL, targets.use = NULL, top = 1, remove.isolate = FALSE,
                                vertex.weight = NULL, vertex.weight.max = NULL, vertex.size.max = 15,
                                weight.scale = TRUE, edge.weight.max = NULL, edge.width.max=8,
                                layout = c("hierarchy","circle","chord"), thresh = 0.05, from = NULL, to = NULL, bidirection = NULL, vertex.size = NULL,
                                pt.title = 12, title.space = 6, vertex.label.cex = 0.8,
                                group = NULL,cell.order = NULL,small.gap = 1, big.gap = 10, scale = FALSE, reduce = -1, show.legend = FALSE, legend.pos.x = 20,legend.pos.y = 20,
                                ...) {
  layout <- match.arg(layout)
  if (!is.null(vertex.size)) {
    warning("'vertex.size' is deprecated. Use `vertex.weight`")
  }
  if (is.null(vertex.weight)) {
    vertex.weight <- as.numeric(table(object@idents))
  }
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
 # prob <-(prob-min(prob))/(max(prob)-min(prob))

  if (layout == "hierarchy") {
    prob.sum <- apply(prob, c(1,2), sum)
   # prob.sum <-(prob.sum-min(prob.sum))/(max(prob.sum)-min(prob.sum))
    if (is.null(edge.weight.max)) {
      edge.weight.max = max(prob.sum)
    }
    par(mfrow=c(1,2), ps = pt.title)
    netVisual_hierarchy1(prob.sum, vertex.receiver = vertex.receiver, sources.use = sources.use, targets.use = targets.use, remove.isolate = remove.isolate, top = top, color.use = color.use, vertex.weight = vertex.weight, vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max, weight.scale = weight.scale, edge.weight.max = edge.weight.max, edge.width.max=edge.width.max, title.name = NULL, vertex.label.cex = vertex.label.cex,...)
    netVisual_hierarchy2(prob.sum, vertex.receiver = setdiff(1:nrow(prob.sum),vertex.receiver), sources.use = sources.use, targets.use = targets.use, remove.isolate = remove.isolate, top = top, color.use = color.use, vertex.weight = vertex.weight, vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max, weight.scale = weight.scale, edge.weight.max = edge.weight.max, edge.width.max=edge.width.max, title.name = NULL, vertex.label.cex = vertex.label.cex,...)
    graphics::mtext(paste0(signaling.name, " signaling pathway network"), side = 3, outer = TRUE, cex = 1, line = -title.space)
    # https://www.andrewheiss.com/blog/2016/12/08/save-base-graphics-as-pseudo-objects-in-r/
    # grid.echo()
    # gg <-  grid.grab()
    gg <- recordPlot()
  } else if (layout == "circle") {
    prob.sum <- apply(prob, c(1,2), sum)
   # prob.sum <-(prob.sum-min(prob.sum))/(max(prob.sum)-min(prob.sum))
    gg <- netVisual_circle(prob.sum, sources.use = sources.use, targets.use = targets.use, remove.isolate = remove.isolate, top = top, color.use = color.use, vertex.weight = vertex.weight, vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max, weight.scale = weight.scale, edge.weight.max = edge.weight.max, edge.width.max=edge.width.max,title.name = paste0(signaling.name, " signaling pathway network"), vertex.label.cex = vertex.label.cex,...)
  } else if (layout == "chord") {
    prob.sum <- apply(prob, c(1,2), sum)
    gg <- netVisual_chord_cell_internal(prob.sum, color.use = color.use, sources.use = sources.use, targets.use = targets.use, remove.isolate = remove.isolate,
                                        group = group, cell.order = cell.order,
                                        lab.cex = vertex.label.cex,small.gap = small.gap, big.gap = big.gap,
                                        scale = scale, reduce = reduce,
                                        title.name = paste0(signaling.name, " signaling pathway network"), show.legend = show.legend, legend.pos.x = legend.pos.x, legend.pos.y= legend.pos.y)
  }

return(gg)

}



#' Visualize the inferred signaling network of individual L-R pairs
#'
#' @param object CellChat object
#' @param signaling a signaling pathway name
#' @param signaling.name alternative signaling pathway name to show on the plot
#' @param pairLR.use a char vector or a data frame consisting of one column named "interaction_name", defining the L-R pairs of interest
#' @param color.use the character vector defining the color of each cell group
#' @param vertex.receiver a numeric vector giving the index of the cell groups as targets in the first hierarchy plot
#' @param sources.use a vector giving the index or the name of source cell groups
#' @param targets.use a vector giving the index or the name of target cell groups.
#' @param remove.isolate whether remove the isolate nodes in the communication network
#' @param top the fraction of interactions to show
#' @param weight.scale whether scale the edge weight
#' @param vertex.weight The weight of vertex: either a scale value or a vector
#' @param vertex.weight.max the maximum weight of vertex; defualt = max(vertex.weight)
#' @param vertex.size.max the maximum vertex size for visualization
#' @param vertex.label.cex The label size of vertex in the network
#' @param edge.weight.max the maximum weight of edge; defualt = max(net)
#' @param edge.width.max The maximum edge width for visualization
#' @param layout "hierarchy", "circle" or "chord"
#' @param height height of plot
#' @param thresh threshold of the p-value for determining significant interaction
# #' @param from,to,bidirection Deprecated. Use `sources.use`,`targets.use`
# #' @param vertex.size Deprecated. Use `vertex.weight`
#'
#' Parameters below are set for "chord" diagram. Please also check the function `netVisual_chord_cell` for more parameters.
#' @param group A named group labels for making multiple-group Chord diagrams. The sector names should be used as the names in the vector.
#' The order of group controls the sector orders and if group is set as a factor, the order of levels controls the order of groups.
#' @param cell.order a char vector defining the cell type orders (sector orders)
#' @param small.gap Small gap between sectors.
#' @param big.gap Gap between the different sets of sectors, which are defined in the `group` parameter
#' @param scale scale each sector to same width; default = FALSE; however, it is set to be TRUE when remove.isolate = TRUE
#' @param reduce if the ratio of the width of certain grid compared to the whole circle is less than this value, the grid is removed on the plot. Set it to value less than zero if you want to keep all tiny grid.
#' @param show.legend whether show the figure legend
#' @param legend.pos.x,legend.pos.y adjust the legend position
#' @param nCol number of columns when displaying the figures using "circle" or "chord"
#'
#' @param ... other parameters (e.g.,vertex.label.cex, vertex.label.color, alpha.edge, label.edge, edge.label.color, edge.label.cex, edge.curved)
#'  passing to `netVisual_hierarchy1`,`netVisual_hierarchy2`,`netVisual_circle`. NB: some parameters might be not supported
#' @importFrom grDevices dev.off pdf
#'
#' @return  an object of class "recordedplot"
#' @export
#'
#'
netVisual_individual <- function(object, signaling, signaling.name = NULL, pairLR.use = NULL, color.use = NULL, vertex.receiver = NULL, sources.use = NULL, targets.use = NULL, top = 1, remove.isolate = FALSE,
                                 vertex.weight = NULL, vertex.weight.max = NULL, vertex.size.max = 15, vertex.label.cex = 0.8,
                                 weight.scale = TRUE, edge.weight.max = NULL, edge.width.max=8,
                                 layout = c("hierarchy","circle","chord"), height = 5, thresh = 0.05, #from = NULL, to = NULL, bidirection = NULL,vertex.size = NULL,
                                 group = NULL,cell.order = NULL,small.gap = 1, big.gap = 10, scale = FALSE, reduce = -1, show.legend = FALSE, legend.pos.x = 20, legend.pos.y = 20, nCol = NULL,
                                 ...) {
  layout <- match.arg(layout)
  # if (!is.null(vertex.size)) {
  #   warning("'vertex.size' is deprecated. Use `vertex.weight`")
  # }
  if (is.null(vertex.weight)) {
    vertex.weight <- as.numeric(table(object@idents))
  }
  pairLR <- searchPair(signaling = signaling, pairLR.use = object@LR$LRsig, key = "pathway_name", matching.exact = T, pair.only = F)

  if (is.null(signaling.name)) {
    signaling.name <- signaling
  }
  net <- object@net

  pairLR.use.name <- dimnames(net$prob)[[3]]
  pairLR.name <- intersect(rownames(pairLR), pairLR.use.name)
  if (!is.null(pairLR.use)) {
    if (is.data.frame(pairLR.use)) {
      pairLR.name <- intersect(pairLR.name, as.character(pairLR.use$interaction_name))
    } else {
      pairLR.name <- intersect(pairLR.name, as.character(pairLR.use))
    }

    if (length(pairLR.name) == 0) {
      stop("There is no significant communication for the input L-R pairs!")
    }
  }

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

  if (is.null(nCol)) {
    nCol <- min(length(pairLR.name.use), 2)
  }

  if (length(dim(prob)) == 2) {
    prob <- replicate(1, prob, simplify="array")
    pval <- replicate(1, pval, simplify="array")
  }

 # prob <-(prob-min(prob))/(max(prob)-min(prob))
  if (is.null(edge.weight.max)) {
    edge.weight.max = max(prob)
  }

  if (layout == "hierarchy") {
    par(mfrow=c(nRow,2), mar = c(5, 4, 4, 2) +0.1)
    for (i in 1:length(pairLR.name.use)) {
      signalName_i <- pairLR$interaction_name_2[i]
      prob.i <- prob[,,i]
      netVisual_hierarchy1(prob.i, vertex.receiver = vertex.receiver, sources.use = sources.use, targets.use = targets.use, remove.isolate = remove.isolate, top = top, color.use = color.use, vertex.weight = vertex.weight, vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max, weight.scale = weight.scale, edge.weight.max = edge.weight.max, edge.width.max=edge.width.max, title.name = signalName_i,...)
      netVisual_hierarchy2(prob.i, vertex.receiver = setdiff(1:nrow(prob.i),vertex.receiver), sources.use = sources.use, targets.use = targets.use, remove.isolate = remove.isolate, top = top, color.use = color.use, vertex.weight = vertex.weight, vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max, weight.scale = weight.scale, edge.weight.max = edge.weight.max, edge.width.max=edge.width.max, title.name = signalName_i,...)
    }
    # grid.echo()
    # gg <-  grid.grab()
    gg <- recordPlot()

  } else if (layout == "circle") {
   # par(mfrow=c(nRow,1))
    par(mfrow = c(ceiling(length(pairLR.name.use)/nCol), nCol), xpd=TRUE)
    gg <- vector("list", length(pairLR.name.use))
    for (i in 1:length(pairLR.name.use)) {
      signalName_i <- pairLR$interaction_name_2[i]
      prob.i <- prob[,,i]
      gg[[i]] <- netVisual_circle(prob.i, sources.use = sources.use, targets.use = targets.use, remove.isolate = remove.isolate, top = top, color.use = color.use, vertex.weight = vertex.weight, vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max, weight.scale = weight.scale, edge.weight.max = edge.weight.max, edge.width.max=edge.width.max, title.name = signalName_i,...)
    }
  } else if (layout == "chord") {
    par(mfrow = c(ceiling(length(pairLR.name.use)/nCol), nCol), xpd=TRUE)
    gg <- vector("list", length(pairLR.name.use))
    for (i in 1:length(pairLR.name.use)) {
      title.name <- pairLR$interaction_name_2[i]
      net <- prob[,,i]
      gg[[i]] <- netVisual_chord_cell_internal(net, color.use = color.use, sources.use = sources.use, targets.use = targets.use, remove.isolate = remove.isolate,
                                               group = group, cell.order = cell.order,
                                               lab.cex = vertex.label.cex,small.gap = small.gap, big.gap = big.gap,
                                               scale = scale, reduce = reduce,
                                               title.name = title.name, show.legend = show.legend, legend.pos.x = legend.pos.x, legend.pos.y = legend.pos.y)
    }
  }
  return(gg)
}



#' Hierarchy plot of cell-cell communications sending to cell groups in vertex.receiver
#'
#' The width of edges represent the strength of the communication.
#'
#' @param net a weighted matrix defining the signaling network
#' @param vertex.receiver  a numeric vector giving the index of the cell groups as targets in the first hierarchy plot
#' @param color.use the character vector defining the color of each cell group
#' @param title.name alternative signaling pathway name to show on the plot
#' @param sources.use a vector giving the index or the name of source cell groups
#' @param targets.use a vector giving the index or the name of target cell groups.
#' @param remove.isolate whether remove the isolate nodes in the communication network
#' @param top the fraction of interactions to show
#' @param weight.scale whether rescale the edge weights
#' @param vertex.weight The weight of vertex: either a scale value or a vector
#' @param vertex.weight.max the maximum weight of vertex; defualt = max(vertex.weight)
#' @param vertex.size.max the maximum vertex size for visualization
#' @param edge.weight.max the maximum weight of edge; defualt = max(net)
#' @param edge.width.max The maximum edge width for visualization
#' @param label.dist the distance between labels and dot position
#' @param space.v the space between different columns in the plot
#' @param space.h the space between different rows in the plot
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
#' @param label.edge whether label edge
#' @param edge.label.color The color for single arrow
#' @param edge.label.cex The size of label for arrows
#' @param vertex.size Deprecated. Use `vertex.weight`
#' @importFrom igraph graph_from_adjacency_matrix ends E V layout_
#' @importFrom grDevices adjustcolor recordPlot
#' @importFrom shape Arrows
#' @return  an object of class "recordedplot"
#' @export
netVisual_hierarchy1 <- function(net, vertex.receiver, color.use = NULL, title.name = NULL,  sources.use = NULL, targets.use = NULL, remove.isolate = FALSE, top = 1,
                                 weight.scale = FALSE, vertex.weight=20, vertex.weight.max = NULL, vertex.size.max = 15,
                                 edge.weight.max = NULL, edge.width.max=8, alpha.edge = 0.6,
                                 label.dist = 2.8, space.v = 1.5, space.h = 1.6, shape= NULL, label.edge=FALSE,edge.curved=0, margin=0.2,
                                vertex.label.cex=0.6,vertex.label.color= "black",arrow.width=1,arrow.size = 0.2,edge.label.color='black',edge.label.cex=0.5, vertex.size = NULL){
  if (!is.null(vertex.size)) {
    warning("'vertex.size' is deprecated. Use `vertex.weight`")
  }
  options(warn = -1)
  thresh <- stats::quantile(net, probs = 1-top)
  net[net < thresh] <- 0

  if ((!is.null(sources.use)) | (!is.null(targets.use))) {
    df.net <- reshape2::melt(net, value.name = "value")
    colnames(df.net)[1:2] <- c("source","target")
    # keep the interactions associated with sources and targets of interest
    if (!is.null(sources.use)){
      if (is.numeric(sources.use)) {
        sources.use <- levels(object@idents)[sources.use]
      }
      df.net <- subset(df.net, source %in% sources.use)
    }
    if (!is.null(targets.use)){
      if (is.numeric(targets.use)) {
        targets.use <- levels(object@idents)[targets.use]
      }
      df.net <- subset(df.net, target %in% targets.use)
    }
    cells.level <- levels(object@idents)
    df.net$source <- factor(df.net$source, levels = cells.level)
    df.net$target <- factor(df.net$target, levels = cells.level)
    df.net$value[is.na(df.net$value)] <- 0
    net <- tapply(df.net[["value"]], list(df.net[["source"]], df.net[["target"]]), sum)
  }
  net[is.na(net)] <- 0

  if (remove.isolate) {
    idx1 <- which(Matrix::rowSums(net) == 0)
    idx2 <- which(Matrix::colSums(net) == 0)
    idx <- intersect(idx1, idx2)
    net <- net[-idx, ]
    net <- net[, -idx]
  }

  if (is.null(color.use)) {
    color.use <- scPalette(nrow(net))
  }

  if (is.null(vertex.weight.max)) {
    vertex.weight.max <- max(vertex.weight)
  }
  vertex.weight <- vertex.weight/vertex.weight.max*vertex.size.max+6

  m <- length(vertex.receiver)
  net2 <- net
  reorder.row <- c(vertex.receiver, setdiff(1:nrow(net),vertex.receiver))
  net2 <- net2[reorder.row,vertex.receiver]
  # Expand out to symmetric (M+N)x(M+N) matrix
  m1 <- nrow(net2); n1 <- ncol(net2)
  net3 <- rbind(cbind(matrix(0, m1, m1), net2), matrix(0, n1, m1+n1))

  row.names(net3) <- c(row.names(net)[vertex.receiver], row.names(net)[setdiff(1:m1,vertex.receiver)], rep("",m))
  colnames(net3) <- row.names(net3)
  color.use3 <- c(color.use[vertex.receiver], color.use[setdiff(1:m1,vertex.receiver)], rep("#FFFFFF",length(vertex.receiver)))
  color.use3.frame <- c(color.use[vertex.receiver], color.use[setdiff(1:m1,vertex.receiver)], color.use[vertex.receiver])

  if (length(vertex.weight) != 1) {
    vertex.weight = c(vertex.weight[vertex.receiver], vertex.weight[setdiff(1:m1,vertex.receiver)],vertex.weight[vertex.receiver])
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

  igraph::V(g)$size<-vertex.weight
  igraph::V(g)$color<-color.use3[igraph::V(g)]
  igraph::V(g)$frame.color <- color.use3.frame[igraph::V(g)]
  igraph::V(g)$label.color <- vertex.label.color
  igraph::V(g)$label.cex<-vertex.label.cex
  if(label.edge){
    E(g)$label<-E(g)$weight
    igraph::E(g)$label <- round(igraph::E(g)$label, digits = 1)
  }
  if (is.null(edge.weight.max)) {
    edge.weight.max <- max(igraph::E(g)$weight)
  }
  if (weight.scale == TRUE) {
    # E(g)$width<-0.3+edge.max.width/(max(E(g)$weight)-min(E(g)$weight))*(E(g)$weight-min(E(g)$weight))
    E(g)$width<- 0.3+E(g)$weight/edge.weight.max*edge.width.max
  }else{
    E(g)$width<-0.3+edge.width.max*E(g)$weight
  }

  E(g)$arrow.width<-arrow.width
  E(g)$arrow.size<-arrow.size
  E(g)$label.color<-edge.label.color
  E(g)$label.cex<-edge.label.cex
  E(g)$color<-adjustcolor(igraph::V(g)$color[edge.start[,1]],alpha.edge)

  label.dist <- c(rep(space.h*label.dist,m), rep(space.h*label.dist, m1-m),rep(0, nrow(net3)-m1))
  label.locs <- c(rep(-pi, m), rep(0, m1-m),rep(-pi, nrow(net3)-m1))
  # text.pos <- cbind(c(-space.h/1.5, space.h/10, space.h/1.2), space.v-space.v/10)
  text.pos <- cbind(c(-space.h/1.5, space.h/22, space.h/1.5), space.v-space.v/7)
  igraph::add.vertex.shape("fcircle", clip=igraph::igraph.shape.noclip,plot=mycircle, parameters=list(vertex.frame.color=1, vertex.frame.width=1))
  plot(g,edge.curved=edge.curved,layout=coords_scale,margin=margin,rescale=T,vertex.shape="fcircle", vertex.frame.width = c(rep(1,m1), rep(2,nrow(net3)-m1)),
       vertex.label.degree=label.locs, vertex.label.dist=label.dist, vertex.label.family="Helvetica")
  text(text.pos, c("Source","Target","Source"), cex = 0.8, col = c("#c51b7d","#c51b7d","#2f6661"))
  arrow.pos1 <- c(-space.h/1.5, space.v-space.v/4, space.h/100000, space.v-space.v/4)
  arrow.pos2 <- c(space.h/1.5, space.v-space.v/4, space.h/20, space.v-space.v/4)
  shape::Arrows(arrow.pos1[1], arrow.pos1[2], arrow.pos1[3], arrow.pos1[4], col = "#c51b7d",arr.lwd = 0.0001,arr.length = 0.2, lwd = 0.8,arr.type="triangle")
  shape::Arrows(arrow.pos2[1], arrow.pos2[2], arrow.pos2[3], arrow.pos2[4], col = "#2f6661",arr.lwd = 0.0001,arr.length = 0.2, lwd = 0.8,arr.type="triangle")
  if (!is.null(title.name)) {
    title.pos = c(space.h/8, space.v)
    text(title.pos[1],title.pos[2],paste0(title.name, " signaling network"), cex = 1)
  }
  # https://www.andrewheiss.com/blog/2016/12/08/save-base-graphics-as-pseudo-objects-in-r/
  # grid.echo()
  # gg <-  grid.grab()
  gg <- recordPlot()
  return(gg)
}


#' Hierarchy plot of cell-cell communication sending to cell groups not in vertex.receiver
#'
#' This function loads the significant interactions as a weighted matrix, and colors
#' represent different types of cells as a structure. The width of edges represent the strength of the communication.
#'
#' @param net a weighted matrix defining the signaling network
#' @param vertex.receiver  a numeric vector giving the index of the cell groups as targets in the first hierarchy plot
#' @param color.use the character vector defining the color of each cell group
#' @param title.name alternative signaling pathway name to show on the plot
#' @param sources.use a vector giving the index or the name of source cell groups
#' @param targets.use a vector giving the index or the name of target cell groups.
#' @param remove.isolate whether remove the isolate nodes in the communication network
#' @param top the fraction of interactions to show
#' @param weight.scale whether rescale the edge weights
#' @param vertex.weight The weight of vertex: either a scale value or a vector
#' @param vertex.weight.max the maximum weight of vertex; defualt = max(vertex.weight)
#' @param vertex.size.max the maximum vertex size for visualization
#' @param edge.weight.max the maximum weight of edge; defualt = max(net)
#' @param edge.width.max The maximum edge width for visualization
#' @param label.dist the distance between labels and dot position
#' @param space.v the space between different columns in the plot
#' @param space.h the space between different rows in the plot
#' @param label.edge Whether or not shows the label of edges (number of connections between different cell types)
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
#' @param vertex.size Deprecated. Use `vertex.weight`
#' @importFrom igraph graph_from_adjacency_matrix ends E V layout_
#' @importFrom grDevices adjustcolor recordPlot
#' @importFrom shape Arrows
#' @return  an object of class "recordedplot"
#' @export
netVisual_hierarchy2 <-function(net, vertex.receiver, color.use = NULL, title.name = NULL, sources.use = NULL, targets.use = NULL, remove.isolate = FALSE, top = 1,
                                weight.scale = FALSE, vertex.weight=20, vertex.weight.max = NULL, vertex.size.max = 15,
                                edge.weight.max = NULL, edge.width.max=8,alpha.edge = 0.6,
                                label.dist = 2.8, space.v = 1.5, space.h = 1.6, shape= NULL, label.edge=FALSE,edge.curved=0, margin=0.2,
                                vertex.label.cex=0.6,vertex.label.color= "black",arrow.width=1,arrow.size = 0.2,edge.label.color='black',edge.label.cex=0.5, vertex.size = NULL){
  if (!is.null(vertex.size)) {
    warning("'vertex.size' is deprecated. Use `vertex.weight`")
  }
  options(warn = -1)
  thresh <- stats::quantile(net, probs = 1-top)
  net[net < thresh] <- 0

  if ((!is.null(sources.use)) | (!is.null(targets.use))) {
    df.net <- reshape2::melt(net, value.name = "value")
    colnames(df.net)[1:2] <- c("source","target")
    # keep the interactions associated with sources and targets of interest
    if (!is.null(sources.use)){
      if (is.numeric(sources.use)) {
        sources.use <- levels(object@idents)[sources.use]
      }
      df.net <- subset(df.net, source %in% sources.use)
    }
    if (!is.null(targets.use)){
      if (is.numeric(targets.use)) {
        targets.use <- levels(object@idents)[targets.use]
      }
      df.net <- subset(df.net, target %in% targets.use)
    }
    cells.level <- levels(object@idents)
    df.net$source <- factor(df.net$source, levels = cells.level)
    df.net$target <- factor(df.net$target, levels = cells.level)
    df.net$value[is.na(df.net$value)] <- 0
    net <- tapply(df.net[["value"]], list(df.net[["source"]], df.net[["target"]]), sum)
  }
  net[is.na(net)] <- 0

  if (remove.isolate) {
    idx1 <- which(Matrix::rowSums(net) == 0)
    idx2 <- which(Matrix::colSums(net) == 0)
    idx <- intersect(idx1, idx2)
    net <- net[-idx, ]
    net <- net[, -idx]
  }


  if (is.null(color.use)) {
    color.use <- scPalette(nrow(net))
  }

  if (is.null(vertex.weight.max)) {
    vertex.weight.max <- max(vertex.weight)
  }
  vertex.weight <- vertex.weight/vertex.weight.max*vertex.size.max+6

  m <- length(vertex.receiver)
  m0 <- nrow(net)-length(vertex.receiver)
  net2 <- net
  reorder.row <- c(setdiff(1:nrow(net),vertex.receiver), vertex.receiver)
  net2 <- net2[reorder.row,vertex.receiver]
  # Expand out to symmetric (M+N)x(M+N) matrix
  m1 <- nrow(net2); n1 <- ncol(net2)
  net3 <- rbind(cbind(matrix(0, m1, m1), net2), matrix(0, n1, m1+n1))
  row.names(net3) <- c(row.names(net)[setdiff(1:m1,vertex.receiver)],row.names(net)[vertex.receiver],  rep("",m))
  colnames(net3) <- row.names(net3)
  color.use3 <- c(color.use[setdiff(1:m1,vertex.receiver)],color.use[vertex.receiver],  rep("#FFFFFF",length(vertex.receiver)))
  color.use3.frame <- c(color.use[setdiff(1:m1,vertex.receiver)], color.use[vertex.receiver], color.use[vertex.receiver])


  if (length(vertex.weight) != 1) {
    vertex.weight = c(vertex.weight[setdiff(1:m1,vertex.receiver)], vertex.weight[vertex.receiver], vertex.weight[vertex.receiver])
  }
  if (is.null(shape)) {
    shape <- rep("circle",nrow(net3))
  }

  g <- graph_from_adjacency_matrix(net3, mode = "directed", weighted = T)
  edge.start <- ends(g, es=igraph::E(g), names=FALSE)
  coords <- matrix(NA, nrow(net3), 2)
  coords[1:m0,1] <- 0; coords[(m0+1):m1,1] <- space.h; coords[(m1+1):nrow(net3),1] <- space.h/2;
  coords[1:m0,2] <- seq(space.v, 0, by = -space.v/(m0-1)); coords[(m0+1):m1,2] <- seq(space.v, 0, by = -space.v/(m1-m0-1));coords[(m1+1):nrow(net3),2] <- seq(space.v, 0, by = -space.v/(n1-1));
  coords_scale<-coords

  igraph::V(g)$size<-vertex.weight
  igraph::V(g)$color<-color.use3[igraph::V(g)]
  igraph::V(g)$frame.color <- color.use3.frame[igraph::V(g)]
  igraph::V(g)$label.color <- vertex.label.color
  igraph::V(g)$label.cex<-vertex.label.cex
  if(label.edge){
    igraph::E(g)$label<-igraph::E(g)$weight
    igraph::E(g)$label <- round(igraph::E(g)$label, digits = 1)
  }
  if (is.null(edge.weight.max)) {
    edge.weight.max <- max(igraph::E(g)$weight)
  }
  if (weight.scale == TRUE) {
   # E(g)$width<-0.3+edge.max.width/(max(E(g)$weight)-min(E(g)$weight))*(E(g)$weight-min(E(g)$weight))
    igraph::E(g)$width<- 0.3+igraph::E(g)$weight/edge.weight.max*edge.width.max
  }else{
    igraph::E(g)$width<-0.3+edge.width.max*igraph::E(g)$weight
  }
  igraph::E(g)$arrow.width<-arrow.width
  igraph::E(g)$arrow.size<-arrow.size
  igraph::E(g)$label.color<-edge.label.color
  igraph::E(g)$label.cex<-edge.label.cex
  igraph::E(g)$color<-adjustcolor(igraph::V(g)$color[edge.start[,1]],alpha.edge)

  label.dist <- c(rep(space.h*label.dist,m), rep(space.h*label.dist, m1-m),rep(0, nrow(net3)-m1))
  label.locs <- c(rep(-pi, m0), rep(0, m1-m0),rep(-pi, nrow(net3)-m1))
  #text.pos <- cbind(c(-space.h/1.5, space.h/10, space.h/1.2), space.v-space.v/10)
  text.pos <- cbind(c(-space.h/1.5, space.h/22, space.h/1.5), space.v-space.v/7)
  igraph::add.vertex.shape("fcircle", clip=igraph::igraph.shape.noclip,plot=mycircle, parameters=list(vertex.frame.color=1, vertex.frame.width=1))
  plot(g,edge.curved=edge.curved,layout=coords_scale,margin=margin,rescale=T,vertex.shape="fcircle", vertex.frame.width = c(rep(1,m1), rep(2,nrow(net3)-m1)),
       vertex.label.degree=label.locs, vertex.label.dist=label.dist, vertex.label.family="Helvetica")
  text(text.pos, c("Source","Target","Source"), cex = 0.8, col = c("#c51b7d","#2f6661","#2f6661"))

  arrow.pos1 <- c(-space.h/1.5, space.v-space.v/4, space.h/100000, space.v-space.v/4)
  arrow.pos2 <- c(space.h/1.5, space.v-space.v/4, space.h/20, space.v-space.v/4)
  shape::Arrows(arrow.pos1[1], arrow.pos1[2], arrow.pos1[3], arrow.pos1[4], col = "#c51b7d",arr.lwd = 0.0001,arr.length = 0.2, lwd = 0.8,arr.type="triangle")
  shape::Arrows(arrow.pos2[1], arrow.pos2[2], arrow.pos2[3], arrow.pos2[4], col = "#2f6661",arr.lwd = 0.0001,arr.length = 0.2, lwd = 0.8,arr.type="triangle")

  if (!is.null(title.name)) {
    title.pos = c(space.h/8, space.v)
    text(title.pos[1],title.pos[2],paste0(title.name, " signaling network"), cex = 1)
  }
  # https://www.andrewheiss.com/blog/2016/12/08/save-base-graphics-as-pseudo-objects-in-r/
  # grid.echo()
  # gg <-  grid.grab()
  gg <- recordPlot()
  return(gg)
}


#' Circle plot of cell-cell communication network
#'
#' The width of edges represent the strength of the communication.
#'
#' @param net A weighted matrix representing the connections
#' @param color.use Colors represent different cell groups
#' @param title.name the name of the title
#' @param sources.use a vector giving the index or the name of source cell groups
#' @param targets.use a vector giving the index or the name of target cell groups.
#' @param remove.isolate whether remove the isolate nodes in the communication network
#' @param top the fraction of interactions to show
#' @param weight.scale whether scale the weight
#' @param vertex.weight The weight of vertex: either a scale value or a vector
#' @param vertex.weight.max the maximum weight of vertex; defualt = max(vertex.weight)
#' @param vertex.size.max the maximum vertex size for visualization
#' @param vertex.label.cex The label size of vertex
#' @param vertex.label.color The color of label for vertex
#' @param edge.weight.max the maximum weight of edge; defualt = max(net)
#' @param edge.width.max The maximum edge width for visualization
#' @param label.edge Whether or not shows the label of edges
#' @param alpha.edge the transprency of edge
#' @param edge.label.color The color for single arrow
#' @param edge.label.cex The size of label for arrows
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
#' @param margin The amount of empty space below, over, at the left and right
#'  of the plot, it is a numeric vector of length four. Usually values between
#'  0 and 0.5 are meaningful, but negative values are also possible, that will
#'  make the plot zoom in to a part of the graph. If it is shorter than four
#'  then it is recycled.
#' @param arrow.width The width of arrows
#' @param arrow.size the size of arrow
# #' @param from,to,bidirection Deprecated. Use `sources.use`,`targets.use`
#' @param vertex.size Deprecated. Use `vertex.weight`
#' @importFrom igraph graph_from_adjacency_matrix ends E V layout_ in_circle
#' @importFrom grDevices recordPlot
#' @return  an object of class "recordedplot"
#' @export
netVisual_circle <-function(net, color.use = NULL,title.name = NULL, sources.use = NULL, targets.use = NULL, remove.isolate = FALSE, top = 1,
                            weight.scale = FALSE, vertex.weight = 20, vertex.weight.max = NULL, vertex.size.max = 15, vertex.label.cex=1,vertex.label.color= "black",
                            edge.weight.max = NULL, edge.width.max=8, alpha.edge = 0.6, label.edge = FALSE,edge.label.color='black',edge.label.cex=0.8,
                            edge.curved=0.2,shape='circle',layout=in_circle(), margin=0.2, vertex.size = NULL,
                            arrow.width=1,arrow.size = 0.2){
  if (!is.null(vertex.size)) {
    warning("'vertex.size' is deprecated. Use `vertex.weight`")
  }
  options(warn = -1)
  thresh <- stats::quantile(net, probs = 1-top)
  net[net < thresh] <- 0

  if ((!is.null(sources.use)) | (!is.null(targets.use))) {
    if (is.null(rownames(net))) {
      stop("The input weighted matrix should have rownames!")
    }
    cells.level <- rownames(net)
    df.net <- reshape2::melt(net, value.name = "value")
    colnames(df.net)[1:2] <- c("source","target")
    # keep the interactions associated with sources and targets of interest
    if (!is.null(sources.use)){
      if (is.numeric(sources.use)) {
        sources.use <- cells.level[sources.use]
      }
      df.net <- subset(df.net, source %in% sources.use)
    }
    if (!is.null(targets.use)){
      if (is.numeric(targets.use)) {
        targets.use <- cells.level[targets.use]
      }
      df.net <- subset(df.net, target %in% targets.use)
    }
    df.net$source <- factor(df.net$source, levels = cells.level)
    df.net$target <- factor(df.net$target, levels = cells.level)
    df.net$value[is.na(df.net$value)] <- 0
    net <- tapply(df.net[["value"]], list(df.net[["source"]], df.net[["target"]]), sum)
  }
  net[is.na(net)] <- 0


  if (remove.isolate) {
    idx1 <- which(Matrix::rowSums(net) == 0)
    idx2 <- which(Matrix::colSums(net) == 0)
    idx <- intersect(idx1, idx2)
    net <- net[-idx, ]
    net <- net[, -idx]
  }

  g <- graph_from_adjacency_matrix(net, mode = "directed", weighted = T)
  edge.start <- igraph::ends(g, es=igraph::E(g), names=FALSE)
  coords<-layout_(g,layout)
  if(nrow(coords)!=1){
    coords_scale=scale(coords)
  }else{
    coords_scale<-coords
  }
  if (is.null(color.use)) {
    color.use = scPalette(length(igraph::V(g)))
  }
  if (is.null(vertex.weight.max)) {
    vertex.weight.max <- max(vertex.weight)
  }
  vertex.weight <- vertex.weight/vertex.weight.max*vertex.size.max+5

  loop.angle<-ifelse(coords_scale[igraph::V(g),1]>0,-atan(coords_scale[igraph::V(g),2]/coords_scale[igraph::V(g),1]),pi-atan(coords_scale[igraph::V(g),2]/coords_scale[igraph::V(g),1]))
  igraph::V(g)$size<-vertex.weight
  igraph::V(g)$color<-color.use[igraph::V(g)]
  igraph::V(g)$frame.color <- color.use[igraph::V(g)]
  igraph::V(g)$label.color <- vertex.label.color
  igraph::V(g)$label.cex<-vertex.label.cex
  if(label.edge){
    igraph::E(g)$label<-igraph::E(g)$weight
    igraph::E(g)$label <- round(igraph::E(g)$label, digits = 1)
  }
  if (is.null(edge.weight.max)) {
    edge.weight.max <- max(igraph::E(g)$weight)
  }
  if (weight.scale == TRUE) {
    #E(g)$width<-0.3+edge.width.max/(max(E(g)$weight)-min(E(g)$weight))*(E(g)$weight-min(E(g)$weight))
    igraph::E(g)$width<- 0.3+igraph::E(g)$weight/edge.weight.max*edge.width.max
  }else{
    igraph::E(g)$width<-0.3+edge.width.max*igraph::E(g)$weight
  }

  igraph::E(g)$arrow.width<-arrow.width
  igraph::E(g)$arrow.size<-arrow.size
  igraph::E(g)$label.color<-edge.label.color
  igraph::E(g)$label.cex<-edge.label.cex
  igraph::E(g)$color<- grDevices::adjustcolor(igraph::V(g)$color[edge.start[,1]],alpha.edge)

  if(sum(edge.start[,2]==edge.start[,1])!=0){
    igraph::E(g)$loop.angle[which(edge.start[,2]==edge.start[,1])]<-loop.angle[edge.start[which(edge.start[,2]==edge.start[,1]),1]]
  }
  radian.rescale <- function(x, start=0, direction=1) {
    c.rotate <- function(x) (x + start) %% (2 * pi) * direction
    c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
  }
  label.locs <- radian.rescale(x=1:length(igraph::V(g)), direction=-1, start=0)
  label.dist <- vertex.weight/max(vertex.weight)+2
  plot(g,edge.curved=edge.curved,vertex.shape=shape,layout=coords_scale,margin=margin, vertex.label.dist=label.dist,
       vertex.label.degree=label.locs, vertex.label.family="Helvetica", edge.label.family="Helvetica") # "sans"
  if (!is.null(title.name)) {
    text(0,1.5,title.name, cex = 1.1)
  }
  # https://www.andrewheiss.com/blog/2016/12/08/save-base-graphics-as-pseudo-objects-in-r/
  # grid.echo()
  # gg <-  grid.grab()
  gg <- recordPlot()
  return(gg)
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


#' Circle plot showing differential cell-cell communication network between two datasets
#'
#' The width of edges represent the relative number of interactions or interaction strength.
#' Red (or blue) colored edges represent increased (or decreased) signaling in the second dataset compared to the first one.
#'
#' @param object A merged CellChat objects
#' @param comparison a numerical vector giving the datasets for comparison in object.list; e.g., comparison = c(1,2)
#' @param measure "count" or "weight". "count": comparing the number of interactions; "weight": comparing the total interaction weights (strength)
#' @param color.use Colors represent different cell groups
#' @param title.name the name of the title
#' @param sources.use a vector giving the index or the name of source cell groups
#' @param targets.use a vector giving the index or the name of target cell groups.
#' @param remove.isolate whether remove the isolate nodes in the communication network
#' @param top the fraction of interactions to show
#' @param weight.scale whether scale the weight
#' @param vertex.weight The weight of vertex: either a scale value or a vector
#' @param vertex.weight.max the maximum weight of vertex; defualt = max(vertex.weight)
#' @param vertex.size.max the maximum vertex size for visualization
#' @param vertex.label.cex The label size of vertex
#' @param vertex.label.color The color of label for vertex
#' @param edge.weight.max the maximum weight of edge; defualt = max(net)
#' @param edge.width.max The maximum edge width for visualization
#' @param label.edge Whether or not shows the label of edges
#' @param alpha.edge the transprency of edge
#' @param edge.label.color The color for single arrow
#' @param edge.label.cex The size of label for arrows
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
#' @param margin The amount of empty space below, over, at the left and right
#'  of the plot, it is a numeric vector of length four. Usually values between
#'  0 and 0.5 are meaningful, but negative values are also possible, that will
#'  make the plot zoom in to a part of the graph. If it is shorter than four
#'  then it is recycled.
#' @param arrow.width The width of arrows
#' @param arrow.size the size of arrow
# #' @param from,to,bidirection Deprecated. Use `sources.use`,`targets.use`
# #' @param vertex.size Deprecated. Use `vertex.weight`
#' @importFrom igraph graph_from_adjacency_matrix ends E V layout_ in_circle
#' @importFrom grDevices recordPlot
#' @return  an object of class "recordedplot"
#' @export
netVisual_diffInteraction <- function(object, comparison = c(1,2), measure = c("count", "weight", "count.merged", "weight.merged"), color.use = NULL,title.name = NULL, sources.use = NULL, targets.use = NULL, remove.isolate = FALSE, top = 1,
                                      weight.scale = FALSE, vertex.weight = 20, vertex.weight.max = NULL, vertex.size.max = 15, vertex.label.cex=1,vertex.label.color= "black",
                                      edge.weight.max = NULL, edge.width.max=8, alpha.edge = 0.6, label.edge = FALSE,edge.label.color='black',edge.label.cex=0.8,
                                      edge.curved=0.2,shape='circle',layout=in_circle(), margin=0.2,
                                      arrow.width=1,arrow.size = 0.2){
  options(warn = -1)
  measure <- match.arg(measure)
  obj1 <- object@net[[comparison[1]]][[measure]]
  obj2 <- object@net[[comparison[2]]][[measure]]
  net.diff <- obj2 - obj1
  if (measure %in% c("count", "count.merged")) {
    if (is.null(title.name)) {
      title.name = "Differential number of interactions"
    }
  } else if (measure %in% c("weight", "weight.merged")) {
    if (is.null(title.name)) {
      title.name = "Differential interaction strength"
    }
  }
  net <- net.diff
  if ((!is.null(sources.use)) | (!is.null(targets.use))) {
    df.net <- reshape2::melt(net, value.name = "value")
    colnames(df.net)[1:2] <- c("source","target")
    # keep the interactions associated with sources and targets of interest
    if (!is.null(sources.use)){
      if (is.numeric(sources.use)) {
        sources.use <- rownames(net.diff)[sources.use]
      }
      df.net <- subset(df.net, source %in% sources.use)
    }
    if (!is.null(targets.use)){
      if (is.numeric(targets.use)) {
        targets.use <- rownames(net.diff)[targets.use]
      }
      df.net <- subset(df.net, target %in% targets.use)
    }
    cells.level <- rownames(net.diff)
    df.net$source <- factor(df.net$source, levels = cells.level)
    df.net$target <- factor(df.net$target, levels = cells.level)
    df.net$value[is.na(df.net$value)] <- 0
    net <- tapply(df.net[["value"]], list(df.net[["source"]], df.net[["target"]]), sum)
  }

  if (remove.isolate) {
    idx1 <- which(Matrix::rowSums(net) == 0)
    idx2 <- which(Matrix::colSums(net) == 0)
    idx <- intersect(idx1, idx2)
    net <- net[-idx, ]
    net <- net[, -idx]
  }

  net[abs(net) < stats::quantile(abs(net), probs = 1-top)] <- 0

  g <- graph_from_adjacency_matrix(net, mode = "directed", weighted = T)
  edge.start <- igraph::ends(g, es=igraph::E(g), names=FALSE)
  coords<-layout_(g,layout)
  if(nrow(coords)!=1){
    coords_scale=scale(coords)
  }else{
    coords_scale<-coords
  }
  if (is.null(color.use)) {
    color.use = scPalette(length(igraph::V(g)))
  }
  if (is.null(vertex.weight.max)) {
    vertex.weight.max <- max(vertex.weight)
  }
  vertex.weight <- vertex.weight/vertex.weight.max*vertex.size.max+5

  loop.angle<-ifelse(coords_scale[igraph::V(g),1]>0,-atan(coords_scale[igraph::V(g),2]/coords_scale[igraph::V(g),1]),pi-atan(coords_scale[igraph::V(g),2]/coords_scale[igraph::V(g),1]))
  igraph::V(g)$size<-vertex.weight
  igraph::V(g)$color<-color.use[igraph::V(g)]
  igraph::V(g)$frame.color <- color.use[igraph::V(g)]
  igraph::V(g)$label.color <- vertex.label.color
  igraph::V(g)$label.cex<-vertex.label.cex
  if(label.edge){
    igraph::E(g)$label<-igraph::E(g)$weight
    igraph::E(g)$label <- round(igraph::E(g)$label, digits = 1)
  }
  igraph::E(g)$arrow.width<-arrow.width
  igraph::E(g)$arrow.size<-arrow.size
  igraph::E(g)$label.color<-edge.label.color
  igraph::E(g)$label.cex<-edge.label.cex
  #igraph::E(g)$color<- grDevices::adjustcolor(igraph::V(g)$color[edge.start[,1]],alpha.edge)
  igraph::E(g)$color <- ifelse(igraph::E(g)$weight > 0,'#b2182b','#2166ac')
  igraph::E(g)$color <- grDevices::adjustcolor(igraph::E(g)$color, alpha.edge)

  igraph::E(g)$weight <- abs(igraph::E(g)$weight)

  if (is.null(edge.weight.max)) {
    edge.weight.max <- max(igraph::E(g)$weight)
  }
  if (weight.scale == TRUE) {
    #E(g)$width<-0.3+edge.width.max/(max(E(g)$weight)-min(E(g)$weight))*(E(g)$weight-min(E(g)$weight))
    igraph::E(g)$width<- 0.3+igraph::E(g)$weight/edge.weight.max*edge.width.max
  }else{
    igraph::E(g)$width<-0.3+edge.width.max*igraph::E(g)$weight
  }


  if(sum(edge.start[,2]==edge.start[,1])!=0){
    igraph::E(g)$loop.angle[which(edge.start[,2]==edge.start[,1])]<-loop.angle[edge.start[which(edge.start[,2]==edge.start[,1]),1]]
  }
  radian.rescale <- function(x, start=0, direction=1) {
    c.rotate <- function(x) (x + start) %% (2 * pi) * direction
    c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
  }
  label.locs <- radian.rescale(x=1:length(igraph::V(g)), direction=-1, start=0)
  label.dist <- vertex.weight/max(vertex.weight)+2
  plot(g,edge.curved=edge.curved,vertex.shape=shape,layout=coords_scale,margin=margin, vertex.label.dist=label.dist,
       vertex.label.degree=label.locs, vertex.label.family="Helvetica", edge.label.family="Helvetica") # "sans"
  if (!is.null(title.name)) {
    text(0,1.5,title.name, cex = 1.1)
  }
  # https://www.andrewheiss.com/blog/2016/12/08/save-base-graphics-as-pseudo-objects-in-r/
  # grid.echo()
  # gg <-  grid.grab()
  gg <- recordPlot()
  return(gg)
}


#' Visualization of network using heatmap
#'
#' This heatmap can be used to show differential number of interactions or interaction strength in the cell-cell communication network between two datasets;
#' the number of interactions or interaction strength in a single dataset
#' the inferred cell-cell communication network in single dataset, defined by `signaling`
#'
#' When show differential number of interactions or interaction strength in the cell-cell communication network between two datasets, the width of edges represent the relative number of interactions or interaction strength.
#' Red (or blue) colored edges represent increased (or decreased) signaling in the second dataset compared to the first one.
#'
#' The top colored bar plot represents the sum of column of values displayed in the heatmap. The right colored bar plot represents the sum of row of values.
#'
#'
#' @param object A merged CellChat object or a single CellChat object
#' @param comparison a numerical vector giving the datasets for comparison in object.list; e.g., comparison = c(1,2)
#' @param measure "count" or "weight". "count": comparing the number of interactions; "weight": comparing the total interaction weights (strength)
#' @param signaling a character vector giving the name of signaling networks in a single CellChat object
#' @param slot.name the slot name of object. Set is to be "netP" if input signaling is a pathway name; Set is to be "net" if input signaling is a ligand-receptor pair
#' @param color.use the character vector defining the color of each cell group
#' @param color.heatmap A vector of two colors corresponding to max/min values, or a color name in brewer.pal only when the data in the heatmap do not contain negative values
#' @param title.name the name of the title
#' @param width width of heatmap
#' @param height height of heatmap
#' @param font.size fontsize in heatmap
#' @param font.size.title font size of the title
#' @param cluster.rows whether cluster rows
#' @param cluster.cols whether cluster columns
#' @param sources.use a vector giving the index or the name of source cell groups
#' @param targets.use a vector giving the index or the name of target cell groups.
#' @param remove.isolate whether remove the isolate nodes in the communication network
#' @param row.show,col.show a vector giving the index or the name of row or columns to show in the heatmap
#' @importFrom methods slot
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#' @importFrom ComplexHeatmap Heatmap HeatmapAnnotation anno_barplot rowAnnotation
#' @return  an object of ComplexHeatmap
#' @export
netVisual_heatmap <- function(object, comparison = c(1,2), measure = c("count", "weight"), signaling = NULL, slot.name = c("netP", "net"), color.use = NULL, color.heatmap = c("#2166ac","#b2182b"),
                              title.name = NULL, width = NULL, height = NULL, font.size = 8, font.size.title = 10, cluster.rows = FALSE, cluster.cols = FALSE,
                              sources.use = NULL, targets.use = NULL, remove.isolate = FALSE, row.show = NULL, col.show = NULL){
  # obj1 <- object.list[[comparison[1]]]
  # obj2 <- object.list[[comparison[2]]]
  if (!is.null(measure)) {
    measure <- match.arg(measure)
  }
  slot.name <- match.arg(slot.name)
  if (is.list(object@net[[1]])) {
    message("Do heatmap based on a merged object \n")
    obj1 <- object@net[[comparison[1]]][[measure]]
    obj2 <- object@net[[comparison[2]]][[measure]]
    net.diff <- obj2 - obj1
    if (measure == "count") {
      if (is.null(title.name)) {
        title.name = "Differential number of interactions"
      }
    } else if (measure == "weight") {
      if (is.null(title.name)) {
        title.name = "Differential interaction strength"
      }
    }
    legend.name = "Relative values"
  } else {
    message("Do heatmap based on a single object \n")
    if (!is.null(signaling)) {
      net.diff <- slot(object, slot.name)$prob[,,signaling]
      if (is.null(title.name)) {
        title.name = paste0(signaling, " signaling network")
      }
      legend.name <- "Communication Prob."
    } else if (!is.null(measure)) {
      net.diff <- object@net[[measure]]
      if (measure == "count") {
        if (is.null(title.name)) {
          title.name = "Number of interactions"
        }
      } else if (measure == "weight") {
        if (is.null(title.name)) {
          title.name = "Interaction strength"
        }
      }
      legend.name <- title.name
    }
  }

  net <- net.diff


  if ((!is.null(sources.use)) | (!is.null(targets.use))) {
    df.net <- reshape2::melt(net, value.name = "value")
    colnames(df.net)[1:2] <- c("source","target")
    # keep the interactions associated with sources and targets of interest
    if (!is.null(sources.use)){
      if (is.numeric(sources.use)) {
        sources.use <- rownames(net.diff)[sources.use]
      }
      df.net <- subset(df.net, source %in% sources.use)
    }
    if (!is.null(targets.use)){
      if (is.numeric(targets.use)) {
        targets.use <- rownames(net.diff)[targets.use]
      }
      df.net <- subset(df.net, target %in% targets.use)
    }
    cells.level <- rownames(net.diff)
    df.net$source <- factor(df.net$source, levels = cells.level)
    df.net$target <- factor(df.net$target, levels = cells.level)
    df.net$value[is.na(df.net$value)] <- 0
    net <- tapply(df.net[["value"]], list(df.net[["source"]], df.net[["target"]]), sum)
  }
  net[is.na(net)] <- 0

  if (remove.isolate) {
    idx1 <- which(Matrix::rowSums(net) == 0)
    idx2 <- which(Matrix::colSums(net) == 0)
    idx <- intersect(idx1, idx2)
    net <- net[-idx, ]
    net <- net[, -idx]
  }

  mat <- net
  if (is.null(color.use)) {
    color.use <- scPalette(ncol(mat))
  }
  names(color.use) <- colnames(mat)

  if (!is.null(row.show)) {
    mat <- mat[row.show, ]
  }
  if (!is.null(col.show)) {
    mat <- mat[ ,col.show]
    color.use <- color.use[col.show]
  }


  if (min(mat) < 0) {
    color.heatmap.use = colorRamp3(c(min(mat), 0, max(mat)), c(color.heatmap[1], "#f7f7f7", color.heatmap[2]))
    colorbar.break <- c(round(min(mat, na.rm = T), digits = nchar(sub(".*\\.(0*).*","\\1",min(mat, na.rm = T)))+1), 0, round(max(mat, na.rm = T), digits = nchar(sub(".*\\.(0*).*","\\1",max(mat, na.rm = T)))+1))
    # color.heatmap.use = colorRamp3(c(seq(min(mat), -(max(mat)-min(max(mat)))/9, length.out = 4), 0, seq((max(mat)-min(max(mat)))/9, max(mat), length.out = 4)), RColorBrewer::brewer.pal(n = 9, name = color.heatmap))
  } else {
    if (length(color.heatmap) == 3) {
      color.heatmap.use = colorRamp3(c(0, min(mat), max(mat)), color.heatmap)
    } else if (length(color.heatmap) == 2) {
      color.heatmap.use = colorRamp3(c(min(mat), max(mat)), color.heatmap)
    } else if (length(color.heatmap) == 1) {
      color.heatmap.use = grDevices::colorRampPalette((RColorBrewer::brewer.pal(n = 9, name = color.heatmap)))(100)
    }
    colorbar.break <- c(round(min(mat, na.rm = T), digits = nchar(sub(".*\\.(0*).*","\\1",min(mat, na.rm = T)))+1), round(max(mat, na.rm = T), digits = nchar(sub(".*\\.(0*).*","\\1",max(mat, na.rm = T)))+1))
  }
  # col_fun(as.vector(mat))

  df<- data.frame(group = colnames(mat)); rownames(df) <- colnames(mat)
  col_annotation <- HeatmapAnnotation(df = df, col = list(group = color.use),which = "column",
                                      show_legend = FALSE, show_annotation_name = FALSE,
                                      simple_anno_size = grid::unit(0.2, "cm"))
  row_annotation <- HeatmapAnnotation(df = df, col = list(group = color.use), which = "row",
                                      show_legend = FALSE, show_annotation_name = FALSE,
                                      simple_anno_size = grid::unit(0.2, "cm"))

  ha1 = rowAnnotation(Strength = anno_barplot(rowSums(abs(mat)), border = FALSE,gp = gpar(fill = color.use, col=color.use)), show_annotation_name = FALSE)
  ha2 = HeatmapAnnotation(Strength = anno_barplot(colSums(abs(mat)), border = FALSE,gp = gpar(fill = color.use, col=color.use)), show_annotation_name = FALSE)

  if (sum(abs(mat) > 0) == 1) {
    color.heatmap.use = c("white", color.heatmap.use)
  } else {
    mat[mat == 0] <- NA
  }
  ht1 = Heatmap(mat, col = color.heatmap.use, na_col = "white", name = legend.name,
                bottom_annotation = col_annotation, left_annotation =row_annotation, top_annotation = ha2, right_annotation = ha1,
                cluster_rows = cluster.rows,cluster_columns = cluster.rows,
                row_names_side = "left",row_names_rot = 0,row_names_gp = gpar(fontsize = font.size),column_names_gp = gpar(fontsize = font.size),
               # width = unit(width, "cm"), height = unit(height, "cm"),
                column_title = title.name,column_title_gp = gpar(fontsize = font.size.title),column_names_rot = 90,
                row_title = "Sources (Sender)",row_title_gp = gpar(fontsize = font.size.title),row_title_rot = 90,
                heatmap_legend_param = list(title_gp = gpar(fontsize = 8, fontface = "plain"),title_position = "leftcenter-rot",
                                            border = NA, #at = colorbar.break,
                                            legend_height = unit(20, "mm"),labels_gp = gpar(fontsize = 8),grid_width = unit(2, "mm"))
  )
  #  draw(ht1)
  return(ht1)
}

#' Show all the significant interactions (L-R pairs) from some cell groups to other cell groups
#'
#' The dot color and size represent the calculated communication probability and p-values.
#'
#' @param object CellChat object
#' @param sources.use a vector giving the index or the name of source cell groups
#' @param targets.use a vector giving the index or the name of target cell groups.
#' @param signaling a character vector giving the name of signaling pathways of interest
#' @param pairLR.use a data frame consisting of one column named either "interaction_name" or "pathway_name", defining the interactions of interest
#' @param color.heatmap A character string or vector indicating the colormap option to use. It can be the avaibale color palette in viridis_pal() or brewer.pal()
#' @param direction Sets the order of colors in the scale. If 1, the default colors are used. If -1, the order of colors is reversed.
#' @param n.colors number of basic colors to generate from color palette
#' @param thresh threshold of the p-value for determining significant interaction
#' @param comparison a numerical vector giving the datasets for comparison in the merged object; e.g., comparison = c(1,2)
#' @param group a numerical vector giving the group information of different datasets; e.g., group = c(1,2,2)
#' @param remove.isolate whether remove the entire empty column, i.e., communication between certain cell groups
#' @param max.dataset a scale, keep the communications with highest probability in max.dataset (i.e., certrain condition)
#' @param min.dataset a scale, keep the communications with lowest probability in min.dataset (i.e., certrain condition)
#' @param min.quantile,max.quantile minimum and maximum quantile cutoff values for the colorbar, may specify quantile in [0,1]
#' @param line.on whether add vertical line when doing comparison analysis for the merged object
#' @param line.size size of vertical line if added
#' @param color.text.use whether color the xtick labels according to the dataset origin when doing comparison analysis
#' @param color.text the colors for xtick labels according to the dataset origin when doing comparison analysis
#' @param title.name main title of the plot
#' @param font.size,font.size.title font size of all the text and the title name
#' @param show.legend whether show legend
#' @param grid.on,color.grid whether add grid
#' @param angle.x,vjust.x,hjust.x parameters for adjusting the rotation of xtick labels
#' @param return.data whether return the data.frame for replotting
#'
#' @return
#' @export
#'
#' @examples
#'\dontrun{
#' # show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
#' netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11), remove.isolate = FALSE)
#'
#' # show all the significant interactions (L-R pairs) associated with certain signaling pathways
#' netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11), signaling = c("CCL","CXCL"))
#'
#' # show all the significant interactions (L-R pairs) based on user's input (defined by `pairLR.use`)
#' pairLR.use <- extractEnrichedLR(cellchat, signaling = c("CCL","CXCL","FGF"))
#' netVisual_bubble(cellchat, sources.use = c(3,4), targets.use = c(5:8), pairLR.use = pairLR.use, remove.isolate = TRUE)
#'
#'# show all the increased interactions in the second dataset compared to the first dataset
#' netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:8), remove.isolate = TRUE, max.dataset = 2)
#'
#'# show all the decreased interactions in the second dataset compared to the first dataset
#' netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:8), remove.isolate = TRUE, max.dataset = 1)
#'}
netVisual_bubble <- function(object, sources.use = NULL, targets.use = NULL, signaling = NULL, pairLR.use = NULL, color.heatmap = c("Spectral","viridis"), n.colors = 10, direction = -1, thresh = 0.05,
                             comparison = NULL, group = NULL, remove.isolate = FALSE, max.dataset = NULL, min.dataset = NULL,
                             min.quantile = 0, max.quantile = 1, line.on = TRUE, line.size = 0.2, color.text.use = TRUE, color.text = NULL,
                             title.name = NULL, font.size = 10, font.size.title = 10, show.legend = TRUE,
                             grid.on = TRUE, color.grid = "grey90", angle.x = 90, vjust.x = NULL, hjust.x = NULL,
                             return.data = FALSE){
  color.heatmap <- match.arg(color.heatmap)
  if (is.list(object@net[[1]])) {
    message("Comparing communications on a merged object \n")
  } else {
    message("Comparing communications on a single object \n")
  }
  if (is.null(vjust.x) | is.null(hjust.x)) {
    angle=c(0, 45, 90)
    hjust=c(0, 1, 1)
    vjust=c(0, 1, 0.5)
    vjust.x = vjust[angle == angle.x]
    hjust.x = hjust[angle == angle.x]
  }
  if (length(color.heatmap) == 1) {
    color.use <- tryCatch({
      RColorBrewer::brewer.pal(n = n.colors, name = color.heatmap)
    }, error = function(e) {
      scales::viridis_pal(option = color.heatmap, direction = -1)(n.colors)
    })
  } else {
    color.use <- color.heatmap
  }
  if (direction == -1) {
    color.use <- rev(color.use)
  }

  if (is.null(comparison)) {
    cells.level <- levels(object@idents)
    if (is.numeric(sources.use)) {
      sources.use <- cells.level[sources.use]
    }
    if (is.numeric(targets.use)) {
      targets.use <- cells.level[targets.use]
    }
    df.net <- subsetCommunication(object, slot.name = "net",
                                  sources.use = sources.use, targets.use = targets.use,
                                  signaling = signaling,
                                  pairLR.use = pairLR.use,
                                  thresh = thresh)
    df.net$source.target <- paste(df.net$source, df.net$target, sep = " -> ")
    source.target <- paste(rep(sources.use, each = length(targets.use)), targets.use, sep = " -> ")
    source.target.isolate <- setdiff(source.target, unique(df.net$source.target))
    if (length(source.target.isolate) > 0) {
      df.net.isolate <- as.data.frame(matrix(NA, nrow = length(source.target.isolate), ncol = ncol(df.net)))
      colnames(df.net.isolate) <- colnames(df.net)
      df.net.isolate$source.target <- source.target.isolate
      df.net.isolate$interaction_name_2 <- df.net$interaction_name_2[1]
      df.net.isolate$pval <- 1
      a <- stringr::str_split(df.net.isolate$source.target, " -> ", simplify = T)
      df.net.isolate$source <- as.character(a[, 1])
      df.net.isolate$target <- as.character(a[, 2])
      df.net <- rbind(df.net, df.net.isolate)
    }

    df.net$pval[df.net$pval > 0.05] = 1
    df.net$pval[df.net$pval > 0.01 & df.net$pval <= 0.05] = 2
    df.net$pval[df.net$pval <= 0.01] = 3
    df.net$prob[df.net$prob == 0] <- NA
    df.net$prob.original <- df.net$prob
    df.net$prob <- -1/log(df.net$prob)

    idx1 <- which(is.infinite(df.net$prob) | df.net$prob < 0)
    if (sum(idx1) > 0) {
      values.assign <- seq(max(df.net$prob, na.rm = T)*1.1, max(df.net$prob, na.rm = T)*1.5, length.out = length(idx1))
      position <- sort(prob.original[idx1], index.return = TRUE)$ix
      df.net$prob[idx1] <- values.assign[match(1:length(idx1), position)]
    }
    # rownames(df.net) <- df.net$interaction_name_2

    df.net$source <- factor(df.net$source, levels = cells.level[cells.level %in% unique(df.net$source)])
    df.net$target <- factor(df.net$target, levels = cells.level[cells.level %in% unique(df.net$target)])
    group.names <- paste(rep(levels(df.net$source), each = length(levels(df.net$target))), levels(df.net$target), sep = " -> ")

    df.net$interaction_name_2 <- as.character(df.net$interaction_name_2)
    df.net <- with(df.net, df.net[order(interaction_name_2),])
    df.net$interaction_name_2 <- factor(df.net$interaction_name_2, levels = unique(df.net$interaction_name_2))
    cells.order <- group.names
    df.net$source.target <- factor(df.net$source.target, levels = cells.order)
    df <- df.net
  } else {
    dataset.name <- names(object@net)
    df.net.all <- subsetCommunication(object, slot.name = "net",
                                      sources.use = sources.use, targets.use = targets.use,
                                      signaling = signaling,
                                      pairLR.use = pairLR.use,
                                      thresh = thresh)
    df.all <- data.frame()
    for (ii in 1:length(comparison)) {
      cells.level <- levels(object@idents[[comparison[ii]]])
      if (is.numeric(sources.use)) {
        sources.use <- cells.level[sources.use]
      }
      if (is.numeric(targets.use)) {
        targets.use <- cells.level[targets.use]
      }

      df.net <- df.net.all[[comparison[ii]]]
      df.net$interaction_name_2 <- as.character(df.net$interaction_name_2)
      df.net$source.target <- paste(df.net$source, df.net$target, sep = " -> ")
      source.target <- paste(rep(sources.use, each = length(targets.use)), targets.use, sep = " -> ")
      source.target.isolate <- setdiff(source.target, unique(df.net$source.target))
      if (length(source.target.isolate) > 0) {
        df.net.isolate <- as.data.frame(matrix(NA, nrow = length(source.target.isolate), ncol = ncol(df.net)))
        colnames(df.net.isolate) <- colnames(df.net)
        df.net.isolate$source.target <- source.target.isolate
        df.net.isolate$interaction_name_2 <- df.net$interaction_name_2[1]
        df.net.isolate$pval <- 1
        a <- stringr::str_split(df.net.isolate$source.target, " -> ", simplify = T)
        df.net.isolate$source <- as.character(a[, 1])
        df.net.isolate$target <- as.character(a[, 2])
        df.net <- rbind(df.net, df.net.isolate)
      }

      df.net$source <- factor(df.net$source, levels = cells.level[cells.level %in% unique(df.net$source)])
      df.net$target <- factor(df.net$target, levels = cells.level[cells.level %in% unique(df.net$target)])
      group.names <- paste(rep(levels(df.net$source), each = length(levels(df.net$target))), levels(df.net$target), sep = " -> ")
      group.names0 <- group.names
      group.names <- paste0(group.names0, " (", dataset.name[comparison[ii]], ")")

      if (nrow(df.net) > 0) {
        df.net$pval[df.net$pval > 0.05] = 1
        df.net$pval[df.net$pval > 0.01 & df.net$pval <= 0.05] = 2
        df.net$pval[df.net$pval <= 0.01] = 3
        df.net$prob[df.net$prob == 0] <- NA
        df.net$prob.original <- df.net$prob
        df.net$prob <- -1/log(df.net$prob)
      } else {
        df.net <- as.data.frame(matrix(NA, nrow = length(group.names), ncol = 5))
        colnames(df.net) <- c("interaction_name_2","source.target","prob","pval","prob.original")
        df.net$source.target <- group.names0
      }
      # df.net$group.names <- sub(paste0(' \\(',dataset.name[comparison[ii]],'\\)'),'',as.character(df.net$source.target))
      df.net$group.names <- as.character(df.net$source.target)
      df.net$source.target <- paste0(df.net$source.target, " (", dataset.name[comparison[ii]], ")")
      df.net$dataset <- dataset.name[comparison[ii]]
      df.all <- rbind(df.all, df.net)
    }
    if (nrow(df.all) == 0) {
      stop("No interactions are detected. Please consider changing the cell groups for analysis. ")
    }

    idx1 <- which(is.infinite(df.all$prob) | df.all$prob < 0)
    if (sum(idx1) > 0) {
      values.assign <- seq(max(df.all$prob, na.rm = T)*1.1, max(df.all$prob, na.rm = T)*1.5, length.out = length(idx1))
      position <- sort(df.all$prob.original[idx1], index.return = TRUE)$ix
      df.all$prob[idx1] <- values.assign[match(1:length(idx1), position)]
    }

    df.all$interaction_name_2[is.na(df.all$interaction_name_2)] <- df.all$interaction_name_2[!is.na(df.all$interaction_name_2)][1]

    df <- df.all
    df <- with(df, df[order(interaction_name_2),])
    df$interaction_name_2 <- factor(df$interaction_name_2, levels = unique(df$interaction_name_2))

    cells.order <- c()
    dataset.name.order <- c()
    for (i in 1:length(group.names0)) {
      for (j in 1:length(comparison)) {
        cells.order <- c(cells.order, paste0(group.names0[i], " (", dataset.name[comparison[j]], ")"))
        dataset.name.order <- c(dataset.name.order, dataset.name[comparison[j]])
      }
    }
    df$source.target <- factor(df$source.target, levels = cells.order)
  }

  min.cutoff <- quantile(df$prob, min.quantile,na.rm= T)
  max.cutoff <- quantile(df$prob, max.quantile,na.rm= T)
  df$prob[df$prob < min.cutoff] <- min.cutoff
  df$prob[df$prob > max.cutoff] <- max.cutoff


  if (remove.isolate) {
    df <- df[!is.na(df$prob), ]
    line.on <- FALSE
  }
  if (!is.null(max.dataset)) {
    # line.on <- FALSE
    # df <- df[!is.na(df$prob),]
    signaling <- as.character(unique(df$interaction_name_2))
    for (i in signaling) {
      df.i <- df[df$interaction_name_2 == i, ,drop = FALSE]
      cell <- as.character(unique(df.i$group.names))
      for (j in cell) {
        df.i.j <- df.i[df.i$group.names == j, , drop = FALSE]
        values <- df.i.j$prob
        idx.max <- which(values == max(values, na.rm = T))
        idx.min <- which(values == min(values, na.rm = T))
        #idx.na <- c(which(is.na(values)), which(!(dataset.name[comparison] %in% df.i.j$dataset)))
        dataset.na <- c(df.i.j$dataset[is.na(values)], setdiff(dataset.name[comparison], df.i.j$dataset))
        if (length(idx.max) > 0) {
          if (!(df.i.j$dataset[idx.max] %in% dataset.name[max.dataset])) {
            df.i.j$prob <- NA
          } else if ((idx.max != idx.min) & !is.null(min.dataset)) {
            if (!(df.i.j$dataset[idx.min] %in% dataset.name[min.dataset])) {
              df.i.j$prob <- NA
            } else if (length(dataset.na) > 0 & sum(!(dataset.name[min.dataset] %in% dataset.na)) > 0) {
              df.i.j$prob <- NA
            }
          }
        }
        df.i[df.i$group.names == j, "prob"] <- df.i.j$prob
      }
      df[df$interaction_name_2 == i, "prob"] <- df.i$prob
    }
    #df <- df[!is.na(df$prob), ]
  }
  if (remove.isolate) {
    df <- df[!is.na(df$prob), ]
    line.on <- FALSE
  }
  if (nrow(df) == 0) {
    stop("No interactions are detected. Please consider changing the cell groups for analysis. ")
  }
  df$interaction_name_2 <- factor(df$interaction_name_2, levels = unique(df$interaction_name_2))
  df$source.target = droplevels(df$source.target, exclude = setdiff(levels(df$source.target),unique(df$source.target)))

  g <- ggplot(df, aes(x = source.target, y = interaction_name_2, color = prob, size = pval)) +
    geom_point(pch = 16) +
    theme_linedraw() + theme(panel.grid.major = element_blank()) +
    theme(axis.text.x = element_text(angle = angle.x, hjust= hjust.x, vjust = vjust.x),
          axis.title.x = element_blank(),
          axis.title.y = element_blank()) +
    scale_x_discrete(position = "bottom")

  values <- c(1,2,3); names(values) <- c("p > 0.05", "0.01 < p < 0.05","p < 0.01")
  g <- g + scale_radius(range = c(min(df$pval), max(df$pval)), breaks = sort(unique(df$pval)),labels = names(values)[values %in% sort(unique(df$pval))], name = "p-value")
  #g <- g + scale_radius(range = c(1,3), breaks = values,labels = names(values), name = "p-value")
  if (min(df$prob, na.rm = T) != max(df$prob, na.rm = T)) {
    g <- g + scale_colour_gradientn(colors = colorRampPalette(color.use)(99), na.value = "white", limits=c(quantile(df$prob, 0,na.rm= T), quantile(df$prob, 1,na.rm= T)),
                                    breaks = c(quantile(df$prob, 0,na.rm= T), quantile(df$prob, 1,na.rm= T)), labels = c("min","max")) +
      guides(color = guide_colourbar(barwidth = 0.5, title = "Commun. Prob."))
  } else {
    g <- g + scale_colour_gradientn(colors = colorRampPalette(color.use)(99), na.value = "white") +
      guides(color = guide_colourbar(barwidth = 0.5, title = "Commun. Prob."))
  }

  g <- g + theme(text = element_text(size = font.size),plot.title = element_text(size=font.size.title)) +
    theme(legend.title = element_text(size = 8), legend.text = element_text(size = 6))

  if (grid.on) {
    if (length(unique(df$source.target)) > 1) {
      g <- g + geom_vline(xintercept=seq(1.5, length(unique(df$source.target))-0.5, 1),lwd=0.1,colour=color.grid)
    }
    if (length(unique(df$interaction_name_2)) > 1) {
      g <- g + geom_hline(yintercept=seq(1.5, length(unique(df$interaction_name_2))-0.5, 1),lwd=0.1,colour=color.grid)
    }
  }
  if (!is.null(title.name)) {
    g <- g + ggtitle(title.name) + theme(plot.title = element_text(hjust = 0.5))
  }

  if (!is.null(comparison)) {
    if (line.on) {
      xintercept = seq(0.5+length(dataset.name[comparison]), length(group.names0)*length(dataset.name[comparison]), by = length(dataset.name[comparison]))
      g <- g + geom_vline(xintercept = xintercept, linetype="dashed", color = "grey60", size = line.size)
    }
    if (color.text.use) {
      if (is.null(group)) {
        group <- 1:length(comparison)
        names(group) <- dataset.name[comparison]
      }
      if (is.null(color.text)) {
        color <- ggPalette(length(unique(group)))
      } else {
        color <- color.text
      }
      names(color) <- names(group[!duplicated(group)])
      color <- color[group]
      #names(color) <- dataset.name[comparison]
      dataset.name.order <- levels(df$source.target)
      dataset.name.order <- stringr::str_match(dataset.name.order, "\\(.*\\)")
      dataset.name.order <- stringr::str_sub(dataset.name.order, 2, stringr::str_length(dataset.name.order)-1)
      xtick.color <- color[dataset.name.order]
      g <- g + theme(axis.text.x = element_text(colour = xtick.color))
    }
  }
  if (!show.legend) {
    g <- g + theme(legend.position = "none")
  }
  if (return.data) {
    return(list(communication = df, gg.obj = g))
  } else {
    return(g)
  }

}



#' Chord diagram for visualizing cell-cell communication for a signaling pathway
#'
#' Names of cell states will be displayed in this chord diagram
#'
#' @param object CellChat object
#' @param signaling a character vector giving the name of signaling networks
#' @param net a weighted matrix or a data frame with three columns defining the cell-cell communication network
#' @param slot.name the slot name of object: slot.name = "net" when visualizing cell-cell communication network per each ligand-receptor pair associated with a given signaling pathway;
#' slot.name = "netP" when visualizing cell-cell communication network at the level of signaling pathways
#' @param color.use colors for the cell groups
#' @param group A named group labels for making multiple-group Chord diagrams. The sector names should be used as the names in the vector.
#' The order of group controls the sector orders and if group is set as a factor, the order of levels controls the order of groups.
#' @param cell.order a char vector defining the cell type orders (sector orders)
#' @param sources.use a vector giving the index or the name of source cell groups
#' @param targets.use a vector giving the index or the name of target cell groups.
#' @param lab.cex font size for the text
#' @param small.gap Small gap between sectors.
#' @param big.gap Gap between the different sets of sectors, which are defined in the `group` parameter
#' @param annotationTrackHeight annotationTrack Height
#' @param remove.isolate whether remove sectors without any links
#' @param link.visible whether plot the link. The value is logical, if it is set to FALSE, the corresponding link will not plotted, but the space is still ocuppied. The format is a matrix with names or a data frame with three columns
#' @param scale scale each sector to same width; default = FALSE; however, it is set to be TRUE when remove.isolate = TRUE
#' @param link.target.prop If the Chord diagram is directional, for each source sector, whether to draw bars that shows the proportion of target sectors.
#' @param reduce if the ratio of the width of certain grid compared to the whole circle is less than this value, the grid is removed on the plot. Set it to value less than zero if you want to keep all tiny grid.
#' @param directional Whether links have directions. 1 means the direction is from the first column in df to the second column, -1 is the reverse, 0 is no direction, and 2 for two directional.
#' @param transparency Transparency of link colors
#' @param link.border border for links, single scalar or a matrix with names or a data frame with three columns
#' @param title.name title name
#' @param show.legend whether show the figure legend
#' @param legend.pos.x,legend.pos.y adjust the legend position
#' @param nCol number of columns when displaying the figures
#' @param thresh threshold of the p-value for determining significant interaction when visualizing links at the level of ligands/receptors;
#' @param ... other parameters passing to chordDiagram
#' @return an object of class "recordedplot"
#' @export

netVisual_chord_cell <- function(object, signaling = NULL, net = NULL, slot.name = "netP",
                                 color.use = NULL,group = NULL,cell.order = NULL,
                                 sources.use = NULL, targets.use = NULL,
                                 lab.cex = 0.8,small.gap = 1, big.gap = 10, annotationTrackHeight = c(0.03),
                                 remove.isolate = FALSE, link.visible = TRUE, scale = FALSE, directional = 1,link.target.prop = TRUE, reduce = -1,
                                 transparency = 0.4, link.border = NA,
                                 title.name = NULL, show.legend = FALSE, legend.pos.x = 20, legend.pos.y = 20, nCol = NULL,
                                 thresh = 0.05,...){

  if (!is.null(signaling)) {
    pairLR <- searchPair(signaling = signaling, pairLR.use = object@LR$LRsig, key = "pathway_name", matching.exact = T, pair.only = F)
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
      stop(paste0('There is no significant communication of ', signaling))
    } else {
      pairLR <- pairLR[pairLR.name.use,]
    }
    nRow <- length(pairLR.name.use)

    prob <- prob[,,pairLR.name.use]

    if (length(dim(prob)) == 2) {
      prob <- replicate(1, prob, simplify="array")
    }

    if (slot.name == "netP") {
      message("Plot the aggregated cell-cell communication network at the signaling pathway level")
      net <- apply(prob, c(1,2), sum)
      if (is.null(title.name)) {
        title.name <- paste0(signaling, " signaling pathway network")
      }
      # par(mfrow = c(1,1), xpd=TRUE)
      # par(mar = c(5, 4, 4, 2))
      gg <- netVisual_chord_cell_internal(net, color.use = color.use, group = group, cell.order = cell.order, sources.use = sources.use, targets.use = targets.use,
                                          lab.cex = lab.cex,small.gap = small.gap, annotationTrackHeight = annotationTrackHeight,
                                          remove.isolate = remove.isolate, link.visible = link.visible, scale = scale, directional = directional,link.target.prop = link.target.prop, reduce = reduce,
                                          transparency = transparency, link.border = link.border,
                                          title.name = title.name, show.legend = show.legend, legend.pos.x = legend.pos.x, legend.pos.y = legend.pos.y, ...)
    } else if (slot.name == "net") {
      message("Plot the cell-cell communication network per each ligand-receptor pair associated with a given signaling pathway")
      if (is.null(nCol)) {
        nCol <- min(length(pairLR.name.use), 2)
      }
      #   layout(matrix(1:length(pairLR.name.use), ncol = nCol))
      # par(xpd=TRUE)
      # par(mfrow = c(ceiling(length(pairLR.name.use)/nCol), nCol), xpd=TRUE, mar = c(5, 4, 4, 2) +0.1)
      par(mfrow = c(ceiling(length(pairLR.name.use)/nCol), nCol), xpd=TRUE)
      gg <- vector("list", length(pairLR.name.use))
      for (i in 1:length(pairLR.name.use)) {
        #par(mar = c(5, 4, 4, 2))
        title.name <- pairLR$interaction_name_2[i]
        net <- prob[,,i]
        gg[[i]] <- netVisual_chord_cell_internal(net, color.use = color.use, group = group,cell.order = cell.order,sources.use = sources.use, targets.use = targets.use,
                                                 lab.cex = lab.cex,small.gap = small.gap, annotationTrackHeight = annotationTrackHeight,
                                                 remove.isolate = remove.isolate, link.visible = link.visible, scale = scale, directional = directional,link.target.prop = link.target.prop, reduce = reduce,
                                                 transparency = transparency, link.border = link.border,
                                                 title.name = title.name, show.legend = show.legend, legend.pos.x = legend.pos.x, legend.pos.y = legend.pos.y, ...)
      }
    }

  } else if (!is.null(net)) {
    gg <- netVisual_chord_cell_internal(net, color.use = color.use, group = group,cell.order = cell.order,sources.use = sources.use, targets.use = targets.use,
                                        lab.cex = lab.cex,small.gap = small.gap, annotationTrackHeight = annotationTrackHeight,
                                        remove.isolate = remove.isolate, link.visible = link.visible, scale = scale, directional = directional,link.target.prop = link.target.prop, reduce = reduce,
                                        transparency = transparency, link.border = link.border,
                                        title.name = title.name, show.legend = show.legend, legend.pos.x = legend.pos.x,legend.pos.y=legend.pos.y, ...)
  } else {
    stop("Please assign values to either `signaling` or `net`")
  }

  return(gg)
}


#' Chord diagram for visualizing cell-cell communication from a weighted adjacency matrix or a data frame
#'
#' Names of cell states/groups will be displayed in this chord diagram
#'
#' @param net a weighted matrix or a data frame with three columns defining the cell-cell communication network
#' @param color.use colors for the cell groups
#' @param group A named group labels for making multiple-group Chord diagrams. The sector names should be used as the names in the vector.
#' The order of group controls the sector orders and if group is set as a factor, the order of levels controls the order of groups.
#' @param cell.order a char vector defining the cell type orders (sector orders)
#' @param sources.use a vector giving the index or the name of source cell groups
#' @param targets.use a vector giving the index or the name of target cell groups.
#' @param lab.cex font size for the text
#' @param small.gap Small gap between sectors.
#' @param big.gap Gap between the different sets of sectors, which are defined in the `group` parameter
#' @param annotationTrackHeight annotationTrack Height
#' @param remove.isolate whether remove sectors without any links
#' @param link.visible whether plot the link. The value is logical, if it is set to FALSE, the corresponding link will not plotted, but the space is still ocuppied. The format is a matrix with names or a data frame with three columns
#' @param scale scale each sector to same width; default = FALSE; however, it is set to be TRUE when remove.isolate = TRUE
#' @param link.target.prop If the Chord diagram is directional, for each source sector, whether to draw bars that shows the proportion of target sectors.
#' @param reduce if the ratio of the width of certain grid compared to the whole circle is less than this value, the grid is removed on the plot. Set it to value less than zero if you want to keep all tiny grid.
#' @param directional Whether links have directions. 1 means the direction is from the first column in df to the second column, -1 is the reverse, 0 is no direction, and 2 for two directional.
#' @param transparency Transparency of link colors
#' @param link.border border for links, single scalar or a matrix with names or a data frame with three columns
#' @param title.name title name of the plot
#' @param show.legend whether show the figure legend
#' @param legend.pos.x,legend.pos.y adjust the legend position
#' @param ... other parameters passing to chordDiagram
#' @importFrom circlize circos.clear chordDiagram circos.track circos.text get.cell.meta.data
#' @importFrom grDevices recordPlot
#' @importFrom BiocGenerics union
#' @return an object of class "recordedplot"
#' @export

netVisual_chord_cell_internal <- function(net, color.use = NULL, group = NULL, cell.order = NULL,
                                          sources.use = NULL, targets.use = NULL,
                                          lab.cex = 0.8,small.gap = 1, big.gap = 10, annotationTrackHeight = c(0.03),
                                          remove.isolate = FALSE, link.visible = TRUE, scale = FALSE, directional = 1, link.target.prop = TRUE, reduce = -1,
                                          transparency = 0.4, link.border = NA,
                                          title.name = NULL, show.legend = FALSE, legend.pos.x = 20, legend.pos.y = 20,...){
  if (inherits(x = net, what = c("matrix", "Matrix"))) {
    cell.levels <- union(rownames(net), colnames(net))
    net <- reshape2::melt(net, value.name = "prob")
    colnames(net)[1:2] <- c("source","target")
  } else if (is.data.frame(net)) {
    if (all(c("source","target", "prob") %in% colnames(net)) == FALSE) {
      stop("The input data frame must contain three columns named as source, target, prob")
    }
    cell.levels <- as.character(union(net$source,net$target))
  }
  if (!is.null(cell.order)) {
    cell.levels <- cell.order
  }
  net$source <- as.character(net$source)
  net$target <- as.character(net$target)

  # keep the interactions associated with sources and targets of interest
  if (!is.null(sources.use)){
    if (is.numeric(sources.use)) {
      sources.use <- cell.levels[sources.use]
    }
    net <- subset(net, source %in% sources.use)
  }
  if (!is.null(targets.use)){
    if (is.numeric(targets.use)) {
      targets.use <- cell.levels[targets.use]
    }
    net <- subset(net, target %in% targets.use)
  }
  # remove the interactions with zero values
  net <- subset(net, prob > 0)
  # create a fake data if keeping the cell types (i.e., sectors) without any interactions
  if (!remove.isolate) {
    cells.removed <- setdiff(cell.levels, as.character(union(net$source,net$target)))
    if (length(cells.removed) > 0) {
      net.fake <- data.frame(cells.removed, cells.removed, 1e-10*sample(length(cells.removed), length(cells.removed)))
      colnames(net.fake) <- colnames(net)
      net <- rbind(net, net.fake)
      link.visible <- net[, 1:2]
      link.visible$plot <- FALSE
      link.visible$plot[1:(nrow(net) - nrow(net.fake))] <- TRUE
      # directional <- net[, 1:2]
      # directional$plot <- 0
      # directional$plot[1:(nrow(net) - nrow(net.fake))] <- 1
      # link.arr.type = "big.arrow"
      # message("Set scale = TRUE when remove.isolate = FALSE")
      scale = TRUE
    }
  }

  df <- net
  cells.use <- union(df$source,df$target)

  # define grid order
  order.sector <- cell.levels[cell.levels %in% cells.use]

  # define grid color
  if (is.null(color.use)){
    color.use = scPalette(length(cell.levels))
    names(color.use) <- cell.levels
  } else if (is.null(names(color.use))) {
    names(color.use) <- cell.levels
  }
  grid.col <- color.use[order.sector]
  names(grid.col) <- order.sector

  # set grouping information
  if (!is.null(group)) {
    group <- group[names(group) %in% order.sector]
  }

  # define edge color
  edge.color <- color.use[as.character(df$source)]

  if (directional == 0 | directional == 2) {
    link.arr.type = "triangle"
  } else {
    link.arr.type = "big.arrow"
  }

  circos.clear()
  chordDiagram(df,
               order = order.sector,
               col = edge.color,
               grid.col = grid.col,
               transparency = transparency,
               link.border = link.border,
               directional = directional,
               direction.type = c("diffHeight","arrows"),
               link.arr.type = link.arr.type, # link.border = "white",
               annotationTrack = "grid",
               annotationTrackHeight = annotationTrackHeight,
               preAllocateTracks = list(track.height = max(strwidth(order.sector))),
               small.gap = small.gap,
               big.gap = big.gap,
               link.visible = link.visible,
               scale = scale,
               group = group,
               link.target.prop = link.target.prop,
               reduce = reduce,
               ...)
  circos.track(track.index = 1, panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    xplot = get.cell.meta.data("xplot")
    ylim = get.cell.meta.data("ylim")
    sector.name = get.cell.meta.data("sector.index")
    circos.text(mean(xlim), ylim[1], sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5),cex = lab.cex)
  }, bg.border = NA)

  # https://jokergoo.github.io/circlize_book/book/legends.html
  if (show.legend) {
    lgd <- ComplexHeatmap::Legend(at = names(grid.col), type = "grid", legend_gp = grid::gpar(fill = grid.col), title = "Cell State")
    ComplexHeatmap::draw(lgd, x = unit(1, "npc")-unit(legend.pos.x, "mm"), y = unit(legend.pos.y, "mm"), just = c("right", "bottom"))
  }

  if(!is.null(title.name)){
    # title(title.name, cex = 1)
    text(-0, 1.02, title.name, cex=1)
  }
  circos.clear()
  gg <- recordPlot()
  return(gg)
}


#' Chord diagram for visualizing cell-cell communication for a set of ligands/receptors or signaling pathways
#'
#' Names of ligands/receptors or signaling pathways will be displayed in this chord diagram
#'
#' @param object CellChat object
#' @param slot.name the slot name of object: slot.name = "net" when visualizing links at the level of ligands/receptors; slot.name = "netP" when visualizing links at the level of signaling pathways
#' @param signaling a character vector giving the name of signaling networks
#' @param pairLR.use a data frame consisting of one column named either "interaction_name" or "pathway_name", defining the interactions of interest
#' @param net A data frame consisting of the interactions of interest.
#' net should have at least three columns: "source","target" and "interaction_name" when visualizing links at the level of ligands/receptors;
#' "source","target" and "pathway_name" when visualizing links at the level of signaling pathway; "interaction_name" and "pathway_name" must be the matched names in CellChatDB$interaction.
#' @param sources.use a vector giving the index or the name of source cell groups
#' @param targets.use a vector giving the index or the name of target cell groups.
#' @param color.use colors for the cell groups
#' @param lab.cex font size for the text
#' @param small.gap Small gap between sectors.
#' @param big.gap Gap between the different sets of sectors, which are defined in the `group` parameter
#' @param annotationTrackHeight annotationTrack Height
#' @param link.visible whether plot the link. The value is logical, if it is set to FALSE, the corresponding link will not plotted, but the space is still ocuppied. The format is a matrix with names or a data frame with three columns
#' @param scale scale each sector to same width; default = FALSE; however, it is set to be TRUE when remove.isolate = TRUE
#' @param link.target.prop If the Chord diagram is directional, for each source sector, whether to draw bars that shows the proportion of target sectors.
#' @param reduce if the ratio of the width of certain grid compared to the whole circle is less than this value, the grid is removed on the plot. Set it to value less than zero if you want to keep all tiny grid.
#' @param directional Whether links have directions. 1 means the direction is from the first column in df to the second column, -1 is the reverse, 0 is no direction, and 2 for two directional.
#' @param transparency Transparency of link colors
#' @param link.border border for links, single scalar or a matrix with names or a data frame with three columns
#' @param title.name title name of the plot
#' @param show.legend whether show the figure legend
#' @param legend.pos.x,legend.pos.y adjust the legend position
#' @param thresh threshold of the p-value for determining significant interaction when visualizing links at the level of ligands/receptors;
#' @param ... other parameters to chordDiagram
#' @importFrom circlize circos.clear chordDiagram circos.track circos.text get.cell.meta.data
#' @importFrom dplyr select %>% group_by summarize
#' @importFrom grDevices recordPlot
#' @importFrom stringr str_split
#' @return an object of class "recordedplot"
#' @export

netVisual_chord_gene <- function(object, slot.name = "net", color.use = NULL,
                                 signaling = NULL, pairLR.use = NULL, net = NULL,
                                 sources.use = NULL, targets.use = NULL,
                                 lab.cex = 0.8,small.gap = 1, big.gap = 10, annotationTrackHeight = c(0.03),
                                 link.visible = TRUE, scale = FALSE, directional = 1, link.target.prop = TRUE, reduce = -1,
                                 transparency = 0.4, link.border = NA,
                                 title.name = NULL, legend.pos.x = 20, legend.pos.y = 20, show.legend = TRUE,
                                 thresh = 0.05,
                                 ...){
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

  if (is.null(net)) {
    prob <- slot(object, "net")$prob
    pval <- slot(object, "net")$pval
    prob[pval > thresh] <- 0
    net <- reshape2::melt(prob, value.name = "prob")
    colnames(net)[1:3] <- c("source","target","interaction_name")

    pairLR = dplyr::select(object@LR$LRsig, c("interaction_name_2", "pathway_name",  "ligand",  "receptor" ,"annotation","evidence"))
    idx <- match(net$interaction_name, rownames(pairLR))
    temp <- pairLR[idx,]
    net <- cbind(net, temp)
  }

  if (!is.null(signaling)) {
    pairLR.use <- data.frame()
    for (i in 1:length(signaling)) {
      pairLR.use.i <- searchPair(signaling = signaling[i], pairLR.use = object@LR$LRsig, key = "pathway_name", matching.exact = T, pair.only = T)
      pairLR.use <- rbind(pairLR.use, pairLR.use.i)
    }
  }

  if (!is.null(pairLR.use)){
    if ("interaction_name" %in% colnames(pairLR.use)) {
      net <- subset(net,interaction_name %in% pairLR.use$interaction_name)
    } else if ("pathway_name" %in% colnames(pairLR.use)) {
      net <- subset(net, pathway_name %in% as.character(pairLR.use$pathway_name))
    }
  }

  if (slot.name == "netP") {
    net <- dplyr::select(net, c("source","target","pathway_name","prob"))
    net$source_target <- paste(net$source, net$target, sep = "sourceTotarget")
    net <- net %>% dplyr::group_by(source_target, pathway_name) %>% dplyr::summarize(prob = sum(prob))
    a <- stringr::str_split(net$source_target, "sourceTotarget", simplify = T)
    net$source <- as.character(a[, 1])
    net$target <- as.character(a[, 2])
    net$ligand <- net$pathway_name
    net$receptor <- " "
  }

  # keep the interactions associated with sources and targets of interest
  if (!is.null(sources.use)){
    if (is.numeric(sources.use)) {
      sources.use <- levels(object@idents)[sources.use]
    }
    net <- subset(net, source %in% sources.use)
  } else {
    sources.use <- levels(object@idents)
  }
  if (!is.null(targets.use)){
    if (is.numeric(targets.use)) {
      targets.use <- levels(object@idents)[targets.use]
    }
    net <- subset(net, target %in% targets.use)
  } else {
    targets.use <- levels(object@idents)
  }
  # remove the interactions with zero values
  df <- subset(net, prob > 0)

  if (nrow(df) == 0) {
    stop("No signaling links are inferred! ")
  }

  if (length(unique(net$ligand)) == 1) {
    message("You may try the function `netVisual_chord_cell` for visualizing individual signaling pathway")
  }

  df$id <- 1:nrow(df)
  # deal with duplicated sector names
  ligand.uni <- unique(df$ligand)
  for (i in 1:length(ligand.uni)) {
    df.i <- df[df$ligand == ligand.uni[i], ]
    source.uni <- unique(df.i$source)
    for (j in 1:length(source.uni)) {
      df.i.j <- df.i[df.i$source == source.uni[j], ]
      df.i.j$ligand <- paste0(df.i.j$ligand, paste(rep(' ',j-1),collapse = ''))
      df$ligand[df$id %in% df.i.j$id] <- df.i.j$ligand
    }
  }
  receptor.uni <- unique(df$receptor)
  for (i in 1:length(receptor.uni)) {
    df.i <- df[df$receptor == receptor.uni[i], ]
    target.uni <- unique(df.i$target)
    for (j in 1:length(target.uni)) {
      df.i.j <- df.i[df.i$target == target.uni[j], ]
      df.i.j$receptor <- paste0(df.i.j$receptor, paste(rep(' ',j-1),collapse = ''))
      df$receptor[df$id %in% df.i.j$id] <- df.i.j$receptor
    }
  }

  cell.order.sources <- levels(object@idents)[levels(object@idents) %in% sources.use]
  cell.order.targets <- levels(object@idents)[levels(object@idents) %in% targets.use]

  df$source <- factor(df$source, levels = cell.order.sources)
  df$target <- factor(df$target, levels = cell.order.targets)
  # df.ordered.source <- df[with(df, order(source, target, -prob)), ]
  # df.ordered.target <- df[with(df, order(target, source, -prob)), ]
  df.ordered.source <- df[with(df, order(source, -prob)), ]
  df.ordered.target <- df[with(df, order(target, -prob)), ]

  order.source <- unique(df.ordered.source[ ,c('ligand','source')])
  order.target <- unique(df.ordered.target[ ,c('receptor','target')])

  # define sector order
  order.sector <- c(order.source$ligand, order.target$receptor)

  # define cell type color
  if (is.null(color.use)){
    color.use = scPalette(nlevels(object@idents))
    names(color.use) <- levels(object@idents)
    color.use <- color.use[levels(object@idents) %in% as.character(union(df$source,df$target))]
  } else if (is.null(names(color.use))) {
    names(color.use) <- levels(object@idents)
    color.use <- color.use[levels(object@idents) %in% as.character(union(df$source,df$target))]
  }

  # define edge color
  edge.color <- color.use[as.character(df.ordered.source$source)]
  names(edge.color) <- as.character(df.ordered.source$source)

  # define grid colors
  grid.col.ligand <- color.use[as.character(order.source$source)]
  names(grid.col.ligand) <- as.character(order.source$source)
  grid.col.receptor <- color.use[as.character(order.target$target)]
  names(grid.col.receptor) <- as.character(order.target$target)
  grid.col <- c(as.character(grid.col.ligand), as.character(grid.col.receptor))
  names(grid.col) <- order.sector

  df.plot <- df.ordered.source[ ,c('ligand','receptor','prob')]

  if (directional == 2) {
    link.arr.type = "triangle"
  } else {
    link.arr.type = "big.arrow"
  }
  circos.clear()
  chordDiagram(df.plot,
               order = order.sector,
               col = edge.color,
               grid.col = grid.col,
               transparency = transparency,
               link.border = link.border,
               directional = directional,
               direction.type = c("diffHeight","arrows"),
               link.arr.type = link.arr.type,
               annotationTrack = "grid",
               annotationTrackHeight = annotationTrackHeight,
               preAllocateTracks = list(track.height = max(strwidth(order.sector))),
               small.gap = small.gap,
               big.gap = big.gap,
               link.visible = link.visible,
               scale = scale,
               link.target.prop = link.target.prop,
               reduce = reduce,
               ...)

  circos.track(track.index = 1, panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    xplot = get.cell.meta.data("xplot")
    ylim = get.cell.meta.data("ylim")
    sector.name = get.cell.meta.data("sector.index")
    circos.text(mean(xlim), ylim[1], sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5),cex = lab.cex)
  }, bg.border = NA)

  # https://jokergoo.github.io/circlize_book/book/legends.html
  if (show.legend) {
    lgd <- ComplexHeatmap::Legend(at = names(color.use), type = "grid", legend_gp = grid::gpar(fill = color.use), title = "Cell State")
    ComplexHeatmap::draw(lgd, x = unit(1, "npc")-unit(legend.pos.x, "mm"), y = unit(legend.pos.y, "mm"), just = c("right", "bottom"))
  }

  circos.clear()
  if(!is.null(title.name)){
    text(-0, 1.02, title.name, cex=1)
  }
  gg <- recordPlot()
  return(gg)
}




#' River plot showing the associations of latent patterns with cell groups and ligand-receptor pairs or signaling pathways
#'
#' River (alluvial) plot shows the correspondence between the inferred latent patterns and cell groups as well as ligand-receptor pairs or signaling pathways.
#'
#' The thickness of the flow indicates the contribution of the cell group or signaling pathway to each latent pattern. The height of each pattern is proportional to the number of its associated cell groups or signaling pathways.
#'
#' Outgoing patterns reveal how the sender cells coordinate with each other as well as how they coordinate with certain signaling pathways to drive communication.
#'
#' Incoming patterns show how the target cells coordinate with each other as well as how they coordinate with certain signaling pathways to respond to incoming signaling.
#'
#' @param object CellChat object
#' @param slot.name the slot name of object that is used to compute centrality measures of signaling networks
#' @param pattern "outgoing" or "incoming"
#' @param cutoff the threshold for filtering out weak links
#' @param sources.use a vector giving the index or the name of source cell groups of interest
#' @param targets.use a vector giving the index or the name of target cell groups of interest
#' @param signaling a character vector giving the name of signaling pathways of interest
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
                              sources.use = NULL, targets.use = NULL, signaling = NULL,
                              color.use = NULL, color.use.pattern = NULL, color.use.signaling = "grey50",
                              do.order = FALSE, main.title = NULL,
                              font.size = 2.5, font.size.title = 12){
  message("Please make sure you have load `library(ggalluvial)` when running this function")
  requireNamespace("ggalluvial")
  #  suppressMessages(require(ggalluvial))
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
    cells.level = levels(object@idents)
    if (is.null(color.use)) {
      color.use <- scPalette(length(cells.level))[cells.level %in% unique(plot.data$CellGroup)]
    }
    if (is.null(color.use.pattern)){
      color.use.pattern <- ggPalette(nPatterns)
    }
    if (!is.null(sources.use)) {
      if (is.numeric(sources.use)) {
        sources.use <- cells.level[sources.use]
      }
      plot.data <- subset(plot.data, CellGroup %in% sources.use)
    }
    if (!is.null(targets.use)) {
      if (is.numeric(targets.use)) {
        targets.use <- cells.level[targets.use]
      }
      plot.data <- subset(plot.data, CellGroup %in% targets.use)
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

    if (!is.null(signaling)) {
      plot.data <- plot.data[plot.data$Signaling %in% signaling, ]
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

    # data3 = merge(data1, data2, by.x="Pattern", by.y="Pattern")
    # data3$Contribution <- data3$Contribution.x * data3$Contribution.y
    # data3 <- data3[,colnames(data3) %in% c("CellGroup","Signaling","Contribution")]

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
#' Using a contribution score of each cell group to each signaling pathway computed by multiplying W by H obtained from `identifyCommunicationPatterns`, we constructed a dot plot in which the dot size is proportion to the contribution score to show association between cell group and their enriched signaling pathways.
#'
#' @param object CellChat object
#' @param slot.name the slot name of object that is used to compute centrality measures of signaling networks
#' @param pattern "outgoing" or "incoming"
#' @param cutoff the threshold for filtering out weak links. Default is 1/R where R is the number of latent patterns. We set the elements in W and H to be zero if they are less than `cutoff`.
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
  comparison <- "single"
  comparison.name <- paste(comparison, collapse = "-")

  Y <- methods::slot(object, slot.name)$similarity[[type]]$dr[[comparison.name]]
  Groups <- methods::slot(object, slot.name)$similarity[[type]]$group[[comparison.name]]
  prob <- methods::slot(object, slot.name)$prob
  if (is.null(pathway.remove)) {
    similarity <- methods::slot(object, slot.name)$similarity[[type]]$matrix[[comparison.name]]
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
  comparison <- "single"
  comparison.name <- paste(comparison, collapse = "-")
  Y <- methods::slot(object, slot.name)$similarity[[type]]$dr[[comparison.name]]
  clusters <- methods::slot(object, slot.name)$similarity[[type]]$group[[comparison.name]]
  prob <- methods::slot(object, slot.name)$prob
  if (is.null(pathway.remove)) {
    similarity <- methods::slot(object, slot.name)$similarity[[type]]$matrix[[comparison.name]]
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
#' @param comparison a numerical vector giving the datasets for comparison. Default are all datasets when object is a merged object
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
netVisual_embeddingPairwise <- function(object, slot.name = "netP", type = c("functional","structural"), comparison = NULL, color.use = NULL, point.shape = NULL, pathway.remove = NULL, pathway.remove.show = TRUE, dot.size = c(2, 6), label.size = 2.5, dot.alpha = 0.5,
                                        xlabel = "Dim 1", ylabel = "Dim 2", title = NULL,do.label = T, show.legend = T, show.axes = T) {
  type <- match.arg(type)
  if (is.null(comparison)) {
    comparison <- 1:length(unique(object@meta$datasets))
  }
  cat("2D visualization of signaling networks from datasets", as.character(comparison), '\n')
  comparison.name <- paste(comparison, collapse = "-")

  Y <- methods::slot(object, slot.name)$similarity[[type]]$dr[[comparison.name]]
  clusters <- methods::slot(object, slot.name)$similarity[[type]]$group[[comparison.name]]
  object.names <- setdiff(names(methods::slot(object, slot.name)), "similarity")[comparison]
  prob <- list()
  for (i in 1:length(comparison)) {
    object.net <- methods::slot(object, slot.name)[[comparison[i]]]
    prob[[i]] = object.net$prob
  }

  if (is.null(point.shape)) {
    point.shape <- c(21, 0, 24, 23, 25, 10, 12)
  }

  if (is.null(pathway.remove)) {
    similarity <- methods::slot(object, slot.name)$similarity[[type]]$matrix[[comparison.name]]
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
#' @param comparison a numerical vector giving the datasets for comparison. Default are all datasets when object is a merged object
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
netVisual_embeddingPairwiseZoomIn <- function(object, slot.name = "netP", type = c("functional","structural"), comparison = NULL, color.use = NULL, nCol = 1, point.shape = NULL, pathway.remove = NULL, dot.size = c(2, 6), label.size = 2.8, dot.alpha = 0.5,
                                              xlabel = NULL, ylabel = NULL, do.label = T, show.legend = F, show.axes = T) {

  type <- match.arg(type)
  if (is.null(comparison)) {
    comparison <- 1:length(unique(object@meta$datasets))
  }
  cat("2D visualization of signaling networks from datasets", as.character(comparison), '\n')
  comparison.name <- paste(comparison, collapse = "-")

  Y <- methods::slot(object, slot.name)$similarity[[type]]$dr[[comparison.name]]
  clusters <- methods::slot(object, slot.name)$similarity[[type]]$group[[comparison.name]]

  object.names <- setdiff(names(methods::slot(object, slot.name)), "similarity")[comparison]
  prob <- list()
  for (i in 1:length(comparison)) {
    object.net <- methods::slot(object, slot.name)[[comparison[i]]]
    prob[[i]] = object.net$prob
  }

  if (is.null(point.shape)) {
    point.shape <- c(21, 0, 24, 23, 25, 10, 12)
  }

  if (is.null(pathway.remove)) {
    similarity <- methods::slot(object, slot.name)$similarity[[type]]$matrix[[comparison.name]]
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
#' @importFrom dplyr group_by summarise n %>%
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





#' A Seurat wrapper function for plotting gene expression using violin plot or dot plot
#'
#' This function create a Seurat object from an input CellChat object, and then plot gene expression distribution using a modified violin plot or dot plot based on Seurat's function.
#' Please check \code{\link{StackedVlnPlot}} and \code{\link{dotPlot}} for detailed description of the arguments.
#'
#' USER can extract the signaling genes related to the inferred L-R pairs or signaling pathway using \code{\link{extractEnrichedLR}}, and then plot gene expression using Seurat package.
#'
#' @param object seurat object
#' @param features Features to plot gene expression
#' @param signaling a char vector containing signaling pathway names for searching
#' @param enriched.only whether only return the identified enriched signaling genes in the database. Default = TRUE, returning the significantly enriched signaling interactions
#' @param type violin plot or dot plot
#' @param color.use defining the color for each cell group
#' @param group.by Name of one metadata columns to group (color) cells. Default is the defined cell groups in CellChat object
#' @param ... other arguments passing to either VlnPlot or DotPlot from Seurat package
#' @return
#' @export
#'
#' @examples

plotGeneExpression <- function(object, features = NULL, signaling = NULL, enriched.only = TRUE, type = c("violin", "dot"), color.use = NULL, group.by = NULL, ...) {
  type <- match.arg(type)
  meta <- object@meta
  if (is.list(object@idents)) {
    meta$group.cellchat <- object@idents$joint
  } else {
    meta$group.cellchat <- object@idents
  }
  w10x <- Seurat::CreateSeuratObject(counts = object@data.signaling, meta.data = meta)
  if (is.null(group.by)) {
    group.by <- "group.cellchat"
  }
  Seurat::Idents(w10x) <- group.by
  if (!is.null(features) & !is.null(signaling)) {
    warning("`features` will be used when inputing both `features` and `signaling`!")
  }
  if (!is.null(features)) {
    feature.use <- features
  } else if (!is.null(signaling)) {
    res <- extractEnrichedLR(object, signaling = signaling, geneLR.return = TRUE, enriched.only = enriched.only)
    feature.use <- res$geneLR
  }
  if (type == "violin") {
    gg <- StackedVlnPlot(w10x, features = feature.use, color.use = color.use, ...)
  } else if (type == "dot") {
    gg <- dotPlot(w10x, features = feature.use, ...)
  }
  return(gg)
}



#' Dot plot
#'
#'The size of the dot encodes the percentage of cells within a class, while the color encodes the AverageExpression level across all cells within a class
#'
#' @param object seurat object
#' @param features Features to plot (gene expression, metrics)
#' @param rotation whether rotate the plot
#' @param colormap RColorbrewer palette to use (check available palette using RColorBrewer::display.brewer.all()). default will use customed color palette
#' @param color.direction Sets the order of colours in the scale. If 1, the default, colours are as output by RColorBrewer::brewer.pal(). If -1, the order of colours is reversed.
#' @param idents Which classes to include in the plot (default is all)
#' @param group.by Name of one or more metadata columns to group (color) cells by
#' (for example, orig.ident); pass 'ident' to group by identity class
#' @param split.by Name of a metadata column to split plot by;
#' @param legend.width legend width
#' @param scale whther show x-axis text
#' @param col.min Minimum scaled average expression threshold (everything smaller will be set to this)
#' @param col.max Maximum scaled average expression threshold (everything larger will be set to this)
#' @param dot.scale Scale the size of the points, similar to cex
#' @param assay Name of assay to use, defaults to the active assay
#' @param angle.x angle for x-axis text rotation
#' @param hjust.x adjust x axis text
#' @param angle.y angle for y-axis text rotation
#' @param hjust.y adjust y axis text
#' @param show.legend whether show the legend
#' @param ... Extra parameters passed to DotPlot from Seurat package
#' @return ggplot2 object
#' @export
#'
#' @examples
#' @import ggplot2
dotPlot <- function(object, features, rotation = TRUE, colormap = "OrRd", color.direction = 1, scale = TRUE, col.min = -2.5, col.max = 2.5, dot.scale = 6, assay = "RNA",
                    idents = NULL, group.by = NULL, split.by = NULL, legend.width = 0.5,
                    angle.x = 45, hjust.x = 1, angle.y = 0, hjust.y = 0.5, show.legend = TRUE, ...) {

  gg <- Seurat::DotPlot(object, features = features, assay = assay, cols = c("blue", "red"),
                scale = scale, col.min = col.min, col.max = col.max, dot.scale = dot.scale,
                idents = idents, group.by = group.by, split.by = split.by,...)
  gg <- gg + theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
    theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10), axis.line = element_line(colour = 'black')) +
    theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5))+
    theme(axis.text.x = element_text(angle = angle.x, hjust = hjust.x), axis.text.y = element_text(angle = angle.y, hjust = hjust.y))

  gg <- gg + theme(legend.title = element_text(size = 10), legend.text = element_text(size = 8))
  if (is.null(split.by)) {
    gg <- gg + guides(color = guide_colorbar(barwidth = legend.width, title = "Scaled expression"),size = guide_legend(title = 'Percent expressed'))
  }

  if (rotation) {
    gg <- gg + coord_flip()
  }
  if (!is.null(colormap)) {
    if (is.null(split.by)) {
      gg <- gg + scale_color_distiller(palette = colormap, direction = color.direction, guide = guide_colorbar(title = "Scaled Expression", ticks = T, label = T, barwidth = legend.width), na.value = "lightgrey")
    }
  }
  if (!show.legend) {
    gg <- gg + theme(legend.position = "none")
  }
  return(gg)
}


#' Stacked Violin plot
#'
#' @param object seurat object
#' @param features Features to plot (gene expression, metrics)
#' @param color.use defining the color for each cell group
#' @param colors.ggplot whether use ggplot color scheme; default: colors.ggplot = FALSE
#' @param split.by Name of a metadata column to split plot by;
#' @param idents Which classes to include in the plot (default is all)
#' @param show.text.y whther show y-axis text
#' @param line.size line width in the violin plot
#' @param pt.size size of the dots
#' @param plot.margin adjust the white space between each plot
#' @param angle.x angle for x-axis text rotation
#' @param vjust.x adjust x axis text
#' @param hjust.x adjust x axis text
#' @param ... Extra parameters passed to VlnPlot from Seurat package
#' @return ggplot2 object
#' @export
#'
#' @examples
#' @import ggplot2
#' @importFrom  patchwork wrap_plots
StackedVlnPlot<- function(object, features, idents = NULL, split.by = NULL,
                          color.use = NULL, colors.ggplot = FALSE,
                          angle.x = 90, vjust.x = NULL, hjust.x = NULL, show.text.y = TRUE, line.size = NULL,
                          pt.size = 0,
                          plot.margin = margin(0, 0, 0, 0, "cm"),
                          ...) {
  options(warn=-1)
  if (is.null(color.use)) {
    numCluster <- length(levels(Seurat::Idents(object)))
    if (colors.ggplot) {
      color.use <- NULL
    } else {
      color.use <- scPalette(numCluster)
    }
  }
  if (is.null(vjust.x) | is.null(hjust.x)) {
    angle=c(0, 45, 90)
    hjust=c(0, 1, 1)
    vjust=c(0, 1, 0.5)
    vjust.x = vjust[angle == angle.x]
    hjust.x = hjust[angle == angle.x]
  }

  plot_list<- purrr::map(features, function(x) modify_vlnplot(object = object, features = x, idents = idents, split.by = split.by, cols = color.use, pt.size = pt.size,
                                                              show.text.y = show.text.y, line.size = line.size, ...))

  # Add back x-axis title to bottom plot. patchwork is going to support this?
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(), axis.ticks.x = element_line()) +
    theme(axis.text.x = element_text(angle = angle.x, hjust = hjust.x, vjust = vjust.x)) +
    theme(axis.text.x = element_text(size = 10))

  # change the y-axis tick to only max value
  ymaxs<- purrr::map_dbl(plot_list, extract_max)
  plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x +
                            scale_y_continuous(breaks = c(y)) +
                            expand_limits(y = y))

  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}

#' modified vlnplot
#' @param object Seurat object
#' @param features Features to plot (gene expression, metrics)
#' @param split.by Name of a metadata column to split plot by;
#' @param idents Which classes to include in the plot (default is all)
#' @param cols defining the color for each cell group
#' @param show.text.y whther show y-axis text
#' @param line.size line width in the violin plot
#' @param pt.size size of the dots
#' @param plot.margin adjust the white space between each plot
#' @param ... pass any arguments to VlnPlot in Seurat
#' @import ggplot2
#'
modify_vlnplot<- function(object,
                          features,
                          idents = NULL,
                          split.by = NULL,
                          cols = NULL,
                          show.text.y = TRUE,
                          line.size = NULL,
                          pt.size = 0,
                          plot.margin = margin(0, 0, 0, 0, "cm"),
                          ...) {
  options(warn=-1)
  p<- Seurat::VlnPlot(object, features = features, cols = cols, pt.size = pt.size, idents = idents, split.by = split.by,  ... )  +
    xlab("") + ylab(features) + ggtitle("")
  p <- p + theme(text = element_text(size = 10)) + theme(axis.line = element_line(size=line.size)) +
    theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 8), axis.line.x = element_line(colour = 'black', size=line.size),axis.line.y = element_line(colour = 'black', size= line.size))
  # theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5))
  p <- p + theme(legend.position = "none",
                 plot.title= element_blank(),
                 axis.title.x = element_blank(),
                 axis.text.x = element_blank(),
                 axis.ticks.x = element_blank(),
                 axis.title.y = element_text(size = rel(1), angle = 0),
                 axis.text.y = element_text(size = rel(1)),
                 plot.margin = plot.margin ) +
    theme(axis.text.y = element_text(size = 8))
  p <- p + theme(element_line(size=line.size))

  if (!show.text.y) {
    p <- p + theme(axis.ticks.y=element_blank(), axis.text.y=element_blank())
  }
  return(p)
}

#' extract the max value of the y axis
#' @param p ggplot object
#' @importFrom  ggplot2 ggplot_build
extract_max<- function(p){
  ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}


