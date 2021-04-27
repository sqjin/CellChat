
#' Compute the communication probability/strength between any interacting cell groups
#'
#' @param object CellChat object
#' @param type methods for computing the average gene expression per cell group. By default = "triMean", producing fewer but stronger interactions;
#' When setting `type = "truncatedMean"`, a value should be assigned to 'trim',  producing more interactions.
#' @param trim the fraction (0 to 0.25) of observations to be trimmed from each end of x before the mean is computed
#' @param LR.use a subset of ligand-receptor interactions used in inferring communication network
#' @param raw.use whether use the raw data (i.e., `object@data.signaling`) or the projected data (i.e., `object@data.project`).
#' Set raw.use = FALSE to use the projected data when analyzing single-cell data with shallow sequencing depth because the projected data could help to reduce the dropout effects of signaling genes, in particular for possible zero expression of subunits of ligands/receptors.
#' @param population.size whether consider the proportion of cells in each group across all sequenced cells.
#' Set population.size = FALSE if analyzing sorting-enriched single cells, to remove the potential artifact of population size.
#' Set population.size = TRUE if analyzing unsorted single-cell transcriptomes, with the reason that abundant cell populations tend to send collectively stronger signals than the rare cell populations.
#' @param do.fast whether run the code in a fast version. In the fast version, all the calculation is based on average expression per cell group; in the previous version, calculation of ligand and receptor expression is also dependent on the number of cells.
#'
#' To further speed up on large-scale datasets, 1) USER can downsample the data using the function 'subset' from Seurat package (e.g., pbmc.small <- subset(pbmc, downsample = 500)), or using the function `sketchData` from CellChat, in particular for the large cell clusters;
#'
#' 2) we are still looking for methods on faster calculation of average expression per group
#'
#' @param nboot threshold of p-values
#' @param seed.use set a random seed. By default, set the seed to 1.
#' @param Kh parameter in Hill function
#' @param n parameter in Hill function
#' @importFrom future nbrOfWorkers
#' @importFrom future.apply future_sapply
#' @importFrom pbapply pbsapply
#' @importFrom stats aggregate
#' @importFrom Matrix crossprod
#' @importFrom utils txtProgressBar setTxtProgressBar
#'
#' @return A CellChat object with updated slot 'net':
#'
#' object@net$prob is the inferred communication probability (strength) array, where the first, second and third dimensions represent a source, target and ligand-receptor pair, respectively.
#'
#' USER can access all the inferred cell-cell communications using the function 'subsetCommunication(object)', which returns a data frame.
#'
#' object@net$pval is the corresponding p-values of each interaction
#'
#' @export
#'
computeCommunProb <- function(object, type = c("triMean", "truncatedMean", "median"), trim = NULL, LR.use = NULL, raw.use = TRUE, population.size = FALSE, do.fast = TRUE, nboot = 100, seed.use = 1L, Kh = 0.5, n = 1) {
  type <- match.arg(type)
  FunMean <- switch(type,
                    triMean = triMean,
                    truncatedMean = function(x) mean(x, trim = trim, na.rm = TRUE),
                    median = function(x) median(x, na.rm = TRUE))
  if (raw.use) {
    data <- as.matrix(object@data.signaling)
  } else {
    data <- object@data.project
  }
  if (is.null(LR.use)) {
    pairLR.use <- object@LR$LRsig
  } else {
    pairLR.use <- LR.use
  }
  complex_input <- object@DB$complex
  cofactor_input <- object@DB$cofactor
  my.sapply <- ifelse(
    test = future::nbrOfWorkers() == 1,
    yes = sapply,
    no = future.apply::future_sapply
  )

  ptm = Sys.time()

  pairLRsig <- pairLR.use
  group <- object@idents
  geneL <- as.character(pairLRsig$ligand)
  geneR <- as.character(pairLRsig$receptor)
  nLR <- nrow(pairLRsig)
  numCluster <- nlevels(group)
  if (numCluster != length(unique(group))) {
    stop("Please check `unique(object@idents)` and ensure that the factor levels are correct!
         You may need to drop unused levels using 'droplevels' function. e.g.,
         `meta$labels = droplevels(meta$labels, exclude = setdiff(levels(meta$labels),unique(meta$labels)))`")
  }
  if (all(data[1:5, ] == floor(data[1:5, ]))) {
    stop("Please check your input data matrix and ensure that you use the normalized data instead of count data!")
  }


  data.use <- data/max(data)
  nC <- ncol(data.use)

  if (do.fast) {
    # compute the average expression per group
    data.use.avg <- aggregate(t(data.use), list(group), FUN = FunMean)
    data.use.avg <- t(data.use.avg[,-1])
    colnames(data.use.avg) <- levels(group)
    # compute the expression of ligand or receptor
    dataLavg <- computeExpr_LR(geneL, data.use.avg, complex_input)
    dataRavg <- computeExpr_LR(geneR, data.use.avg, complex_input)
    # take account into the effect of co-activation and co-inhibition receptors
    dataRavg.co.A.receptor <- computeExpr_coreceptor(cofactor_input, data.use.avg, pairLRsig, type = "A")
    dataRavg.co.I.receptor <- computeExpr_coreceptor(cofactor_input, data.use.avg, pairLRsig, type = "I")
    dataRavg <- dataRavg * dataRavg.co.A.receptor/dataRavg.co.I.receptor

    dataLavg2 <- t(replicate(nrow(dataLavg), as.numeric(table(group))/nC))
    dataRavg2 <- dataLavg2

    # compute the expression of agonist and antagonist
    index.agonist <- which(!is.na(pairLRsig$agonist) & pairLRsig$agonist != "")
    index.antagonist <- which(!is.na(pairLRsig$antagonist) & pairLRsig$antagonist != "")
    # quantify the communication probability
    Prob <- array(0, dim = c(numCluster,numCluster,nLR))
    Pval <- array(0, dim = c(numCluster,numCluster,nLR))

    set.seed(seed.use)
    permutation <- replicate(nboot, sample.int(nC, size = nC))
    data.use.avg.boot <- my.sapply(
      X = 1:nboot,
      FUN = function(nE) {
        groupboot <- group[permutation[, nE]]
        data.use.avgB <- aggregate(t(data.use), list(groupboot), FUN = FunMean)
        data.use.avgB <- t(data.use.avgB[,-1])
        return(data.use.avgB)
      },
      simplify = FALSE
    )
    pb <- txtProgressBar(min = 0, max = nLR, style = 3, file = stderr())

    for (i in 1:nLR) {
      # ligand/receptor
      dataLR <- Matrix::crossprod(matrix(dataLavg[i,], nrow = 1), matrix(dataRavg[i,], nrow = 1))
      P1 <- dataLR^n/(Kh^n + dataLR^n)
      if (sum(P1) == 0) {
        Pnull = P1
        Prob[ , , i] <- Pnull
        p = 1
        Pval[, , i] <- matrix(p, nrow = numCluster, ncol = numCluster, byrow = FALSE)
      } else {
        # agonist and antagonist
        if (is.element(i, index.agonist)) {
          data.agonist <- computeExpr_agonist(data.use = data.use.avg, pairLRsig, cofactor_input, index.agonist = i, Kh = Kh,  n = n)
          P2 <- Matrix::crossprod(matrix(data.agonist, nrow = 1))
        } else {
          P2 <- matrix(1, nrow = numCluster, ncol = numCluster)
        }
        if (is.element(i, index.antagonist)) {
          data.antagonist <- computeExpr_antagonist(data.use = data.use.avg, pairLRsig, cofactor_input,  index.antagonist = i, Kh = Kh,  n = n)
          P3 <- Matrix::crossprod(matrix(data.antagonist, nrow = 1))
        } else {
          P3 <- matrix(1, nrow = numCluster, ncol = numCluster)
        }
        # number of cells
        if (population.size) {
          P4 <- Matrix::crossprod(matrix(dataLavg2[i,], nrow = 1), matrix(dataRavg2[i,], nrow = 1))
        } else {
          P4 <- matrix(1, nrow = numCluster, ncol = numCluster)
        }

        Pnull = P1*P2*P3*P4
        Prob[ , , i] <- Pnull

        Pnull <- as.vector(Pnull)

        #Pboot <- foreach(nE = 1:nboot) %dopar% {
        Pboot <- sapply(
          X = 1:nboot,
          FUN = function(nE) {
            data.use.avgB <- data.use.avg.boot[[nE]]
            dataLavgB <- computeExpr_LR(geneL[i], data.use.avgB, complex_input)
            dataRavgB <- computeExpr_LR(geneR[i], data.use.avgB, complex_input)
            # take account into the effect of co-activation and co-inhibition receptors
            dataRavgB.co.A.receptor <- computeExpr_coreceptor(cofactor_input, data.use.avgB, pairLRsig[i, , drop = FALSE], type = "A")
            dataRavgB.co.I.receptor <- computeExpr_coreceptor(cofactor_input, data.use.avgB, pairLRsig[i, , drop = FALSE], type = "I")
            dataRavgB <- dataRavgB * dataRavgB.co.A.receptor/dataRavgB.co.I.receptor
            dataLRB = Matrix::crossprod(dataLavgB, dataRavgB)
            P1.boot <- dataLRB^n/(Kh^n + dataLRB^n)
            # agonist and antagonist
            if (is.element(i, index.agonist)) {
              data.agonist <- computeExpr_agonist(data.use = data.use.avgB, pairLRsig, cofactor_input, index.agonist = i, Kh = Kh,  n = n)
              P2.boot <- Matrix::crossprod(matrix(data.agonist, nrow = 1))
            } else {
              P2.boot <- matrix(1, nrow = numCluster, ncol = numCluster)
            }
            if (is.element(i, index.antagonist)) {
              data.antagonist <- computeExpr_antagonist(data.use = data.use.avgB, pairLRsig, cofactor_input, index.antagonist = i, Kh = Kh,  n= n)
              P3.boot <- Matrix::crossprod(matrix(data.antagonist, nrow = 1))
            } else {
              P3.boot <- matrix(1, nrow = numCluster, ncol = numCluster)
            }

            if (population.size) {
              groupboot <- group[permutation[, nE]]
              dataLavg2B <- as.numeric(table(groupboot))/nC
              dataLavg2B <- matrix(dataLavg2B, nrow = 1)
              dataRavg2B <- dataLavg2B
              P4.boot = Matrix::crossprod(dataLavg2B, dataRavg2B)
            } else {
              P4.boot = matrix(1, nrow = numCluster, ncol = numCluster)
            }

            Pboot = P1.boot*P2.boot*P3.boot*P4.boot
            return(as.vector(Pboot))
          }
        )
        Pboot <- matrix(unlist(Pboot), nrow=length(Pnull), ncol = nboot, byrow = FALSE)
        nReject <- rowSums(Pboot - Pnull >= 0)
        p = nReject/nboot
        Pval[, , i] <- matrix(p, nrow = numCluster, ncol = numCluster, byrow = FALSE)
      }
     setTxtProgressBar(pb = pb, value = i)
    }
    close(con = pb)
  } else {
    # compute the expression of ligand and receptor
    dataL <- computeExpr_LR(geneL, data.use, complex_input)
    dataR <- computeExpr_LR(geneR, data.use, complex_input)
    # take account into the effect of co-activation and co-inhibition receptors
    dataR.co.A.receptor <- computeExpr_coreceptor(cofactor_input, data.use, pairLRsig, type = "A")
    dataR.co.I.receptor <- computeExpr_coreceptor(cofactor_input, data.use, pairLRsig, type = "I")
    dataR <- dataR * dataR.co.A.receptor/dataR.co.I.receptor
    # compute the average expression in each cell group
    dataLavg <- aggregate(t(dataL), list(group), FUN = FunMean)
    dataLavg <- t(dataLavg[,-1])
    rownames(dataLavg) <- geneL
    dataRavg <- aggregate(t(dataR), list(group), FUN = FunMean)
    dataRavg <- t(dataRavg[,-1])
    rownames(dataRavg) <- geneR

    dataL.binary = (dataL > 0)*1 ;dataR.binary = (dataR > 0)*1
    dataLavg2 <- aggregate(t(dataL.binary), list(group), FUN = sum)
    dataLavg2 <- t(dataLavg2[,-1])/nC
    dataRavg2 <- aggregate(t(dataR.binary), list(group), FUN = sum)
    dataRavg2 <- t(dataRavg2[,-1])/nC

    # compute the expression of agonist and antagonist
    index.agonist <- which(!is.na(pairLRsig$agonist) & pairLRsig$agonist != "")
    index.antagonist <- which(!is.na(pairLRsig$antagonist) & pairLRsig$antagonist != "")
    # quantify the communication probability
    set.seed(seed.use)
    permutation <- replicate(nboot, sample.int(nC, size = nC))
    Prob <- array(0, dim = c(numCluster,numCluster,nLR))
    Pval <- array(0, dim = c(numCluster,numCluster,nLR))
    pb <- txtProgressBar(min = 0, max = nLR, style = 3, file = stderr())
    for (i in 1:nLR) {
      # ligand/receptor
      dataLR <- Matrix::crossprod(matrix(dataLavg[i,], nrow = 1), matrix(dataRavg[i,], nrow = 1))
      P1 <- dataLR^n/(Kh^n + dataLR^n)
      if (sum(P1) == 0) {
        Pnull = P1
        Prob[ , , i] <- Pnull
        p = 1
        Pval[, , i] <- matrix(p, nrow = numCluster, ncol = numCluster, byrow = FALSE)
      } else {
        # agonist and antagonist
        if (is.element(i, index.agonist)) {
          data.agonist <- computeExprGroup_agonist(data.use = data.use, pairLRsig, cofactor_input, group = group,index.agonist = i, Kh = Kh, FunMean = FunMean, n = n)
          P2 <- Matrix::crossprod(matrix(data.agonist, nrow = 1))
        } else {
          P2 <- matrix(1, nrow = numCluster, ncol = numCluster)
        }
        if (is.element(i, index.antagonist)) {
          data.antagonist <- computeExprGroup_antagonist(data.use = data.use, pairLRsig, cofactor_input, group = group, index.antagonist = i, Kh = Kh, FunMean = FunMean, n = n)
          P3 <- Matrix::crossprod(matrix(data.antagonist, nrow = 1))
        } else {
          P3 <- matrix(1, nrow = numCluster, ncol = numCluster)
        }
        # number of cells
        if (population.size) {
          P4 <- Matrix::crossprod(matrix(dataLavg2[i,], nrow = 1), matrix(dataRavg2[i,], nrow = 1))
        } else {
          P4 <- matrix(1, nrow = numCluster, ncol = numCluster)
        }

        Pnull = P1*P2*P3*P4
        Prob[ , , i] <- Pnull

        Pnull <- as.vector(Pnull)
        dataL.i <- dataL[i,]; dataR.i <- dataR[i,];
        dataL2.i <- dataL.binary[i,]; dataR2.i <- dataR.binary[i,];
        #Pboot <- foreach(nE = 1:nboot) %dopar% {
        Pboot <- my.sapply(
          X = 1:nboot,
          FUN = function(nE) {
            groupboot <- group[permutation[, nE]]
            dataLavgB <- aggregate(matrix(dataL.i, ncol = 1), list(groupboot), FUN = FunMean)
            dataLavgB <- t(dataLavgB[,-1])
            dataLavgB <- matrix(dataLavgB, nrow = 1)

            dataRavgB <- aggregate(matrix(dataR.i, ncol = 1), list(groupboot), FUN = FunMean)
            dataRavgB <- t(dataRavgB[,-1])
            dataRavgB <- matrix(dataRavgB, nrow = 1)
            dataLRB = Matrix::crossprod(dataLavgB, dataRavgB)
            P1.boot <- dataLRB^n/(Kh^n + dataLRB^n)
            # agonist and antagonist
            if (is.element(i, index.agonist)) {
              data.agonist <- computeExprGroup_agonist(data.use = data.use, pairLRsig, cofactor_input, group = groupboot, index.agonist = i, Kh = Kh, FunMean = FunMean, n = n)
              P2.boot <- Matrix::crossprod(matrix(data.agonist, nrow = 1))
            } else {
              P2.boot <- matrix(1, nrow = numCluster, ncol = numCluster)
            }
            if (is.element(i, index.antagonist)) {
              data.antagonist <- computeExprGroup_antagonist(data.use = data.use, pairLRsig, cofactor_input, group = groupboot,index.antagonist = i, Kh = Kh, FunMean = FunMean, n= n)
              P3.boot <- Matrix::crossprod(matrix(data.antagonist, nrow = 1))
            } else {
              P3.boot <- matrix(1, nrow = numCluster, ncol = numCluster)
            }
            dataLavg2B <- by(matrix(dataL2.i, ncol = 1), groupboot, sum)/nC
            dataLavg2B <- matrix(dataLavg2B, nrow = 1)

            dataRavg2B <- by(matrix(dataR2.i, ncol = 1), groupboot, sum)/nC
            dataRavg2B <- matrix(dataRavg2B, nrow = 1)
            if (population.size) {
              P4.boot = Matrix::crossprod(dataLavg2B, dataRavg2B)
            } else {
              P4.boot = matrix(1, nrow = numCluster, ncol = numCluster)
            }

            Pboot = P1.boot*P2.boot*P3.boot*P4.boot
            return(as.vector(Pboot))
          }
        )
        Pboot <- matrix(unlist(Pboot), nrow=length(Pnull), ncol = nboot, byrow = FALSE)
        nReject <- rowSums(Pboot - Pnull >= 0)
        p = nReject/nboot
        Pval[, , i] <- matrix(p, nrow = numCluster, ncol = numCluster, byrow = FALSE)
      }
      setTxtProgressBar(pb = pb, value = i)
    }
    close(con = pb)
  }

  Pval[Prob == 0] <- 1
  dimnames(Prob) <- list(levels(group), levels(group), rownames(pairLRsig))
  dimnames(Pval) <- dimnames(Prob)
  net <- list("prob" = Prob, "pval" = Pval)
  execution.time = Sys.time() - ptm
  object@options$run.time <- as.numeric(execution.time, units = "secs")
  object@options$parameter <- list(type.mean = type, trim = trim, raw.use = raw.use, population.size = population.size,  nboot = nboot, seed.use = seed.use, Kh = Kh, n = n)
  object@net <- net
  return(object)
}


#' Compute the communication probability on signaling pathway level by summarizing all related ligands/receptors
#'
#' @param object CellChat object
#' @param net A list from object@net; If net = NULL, net = object@net
#' @param pairLR.use A dataframe giving the ligand-receptor interactions; If pairLR.use = NULL, pairLR.use = object@LR$LRsig
#' @param thresh threshold of the p-value for determining significant interaction
#'
#' @return A CellChat object with updated slot 'netP':
#'
#' object@netP$prob is the communication probability array on signaling pathway level; USER can convert this array to a data frame using the function 'reshape2::melt()',
#'
#' e.g., `df.netP <- reshape2::melt(object@netP$prob, value.name = "prob"); colnames(df.netP)[1:3] <- c("source","target","pathway_name")` or access all significant interactions using the function \code{\link{subsetCommunication}}
#'
#' object@netP$pathways list all the signaling pathways with significant communications.
#'
#' From version >= 1.1.0, pathways are ordered based on the total communication probabilities. NB: pathways with small total communication probabilities might be also very important since they might be specifically activated between only few cell types.
#'
#' @export
#'
computeCommunProbPathway <- function(object = NULL, net = NULL, pairLR.use = NULL, thresh = 0.05) {
  if (is.null(net)) {
    net <- object@net
  }
  if (is.null(pairLR.use)) {
    pairLR.use <- object@LR$LRsig
  }
  prob <- net$prob
  prob[net$pval > thresh] <- 0
  pathways <- unique(pairLR.use$pathway_name)
  group <- factor(pairLR.use$pathway_name, levels = pathways)
  prob.pathways <- aperm(apply(prob, c(1, 2), by, group, sum), c(2, 3, 1))
  pathways.sig <- pathways[apply(prob.pathways, 3, sum) != 0]
  prob.pathways.sig <- prob.pathways[,,pathways.sig]
  idx <- sort(apply(prob.pathways.sig, 3, sum), decreasing=TRUE, index.return = TRUE)$ix
  pathways.sig <- pathways.sig[idx]
  prob.pathways.sig <- prob.pathways.sig[, , idx]

  if (is.null(object)) {
    netP = list(pathways = pathways.sig, prob = prob.pathways.sig)
    return(netP)
  } else {
    object@netP$pathways <- pathways.sig
    object@netP$prob <- prob.pathways.sig
    return(object)
  }
}


#' Calculate the aggregated network by counting the number of links or summarizing the communication probability
#'
#' @param object CellChat object
#' @param sources.use,targets.use,signaling,pairLR.use Please check the description in function \code{\link{subsetCommunication}}
#' @param remove.isolate whether removing the isolate cell groups without any interactions when applying \code{\link{subsetCommunication}}
#' @param thresh threshold of the p-value for determining significant interaction
#' @param return.object whether return an updated CellChat object
#' @importFrom  dplyr group_by summarize groups
#' @importFrom stringr str_split
#'
#' @return Return an updated CellChat object:
#'
#' `object@net$count` is a matrix: rows and columns are sources and targets respectively, and elements are the number of interactions between any two cell groups. USER can convert a matrix to a data frame using the function `reshape2::melt()`
#'
#' `object@net$weight` is also a matrix containing the interaction weights between any two cell groups
#'
#' `object@net$sum` is deprecated. Use `object@net$weight`
#'
#' @export
#'
aggregateNet <- function(object, sources.use = NULL, targets.use = NULL, signaling = NULL, pairLR.use = NULL, remove.isolate = TRUE, thresh = 0.05, return.object = TRUE) {
  net <- object@net
  if (is.null(sources.use) & is.null(targets.use) & is.null(signaling) & is.null(pairLR.use)) {
    prob <- net$prob
    pval <- net$pval
    pval[prob == 0] <- 1
    prob[pval >= thresh] <- 0
    net$count <- apply(prob > 0, c(1,2), sum)
    net$weight <- apply(prob, c(1,2), sum)
    net$weight[is.na(net$weight)] <- 0
    net$count[is.na(net$count)] <- 0
  } else {
    df.net <- subsetCommunication(object, slot.name = "net",
                                  sources.use = sources.use, targets.use = targets.use,
                                  signaling = signaling,
                                  pairLR.use = pairLR.use,
                                  thresh = thresh)
    df.net$source_target <- paste(df.net$source, df.net$target, sep = "_")
    df.net2 <- df.net %>% group_by(source_target) %>% summarize(count = n(), .groups = 'drop')
    df.net3 <- df.net %>% group_by(source_target) %>% summarize(prob = sum(prob), .groups = 'drop')
    df.net2$prob <- df.net3$prob
    a <- stringr::str_split(df.net2$source_target, "_", simplify = T)
    df.net2$source <- as.character(a[, 1])
    df.net2$target <- as.character(a[, 2])
    cells.level <- levels(object@idents)
    if (remove.isolate) {
      message("Isolate cell groups without any interactions are removed. To block it, set `remove.isolate = FALSE`")
      df.net2$source <- factor(df.net2$source, levels = cells.level[cells.level %in% unique(df.net2$source)])
      df.net2$target <- factor(df.net2$target, levels = cells.level[cells.level %in% unique(df.net2$target)])
    } else {
      df.net2$source <- factor(df.net2$source, levels = cells.level)
      df.net2$target <- factor(df.net2$target, levels = cells.level)
    }

    count <- tapply(df.net2[["count"]], list(df.net2[["source"]], df.net2[["target"]]), sum)
    prob <- tapply(df.net2[["prob"]], list(df.net2[["source"]], df.net2[["target"]]), sum)
    net$count <- count
    net$weight <- prob
    net$weight[is.na(net$weight)] <- 0
    net$count[is.na(net$count)] <- 0
  }
  if (return.object) {
    object@net <- net
    return(object)
  } else {
    return(net)
  }

}


#' Compute averaged expression values for each cell group
#'
#' @param object CellChat object
#' @param features a char vector giving the used features. default use all features
#' @param group.by cell group information; default is `object@idents` when input is a single object and `object@idents$joint` when input is a merged object; otherwise it should be one of the column names of the meta slot
#' @param type methods for computing the average gene expression per cell group.
#'
#' By default = "triMean", defined as a weighted average of the distribution's median and its two quartiles (https://en.wikipedia.org/wiki/Trimean);
#'
#' When setting `type = "truncatedMean"`, a value should be assigned to 'trim'. See the function `base::mean`.
#'
#' @param trim the fraction (0 to 0.25) of observations to be trimmed from each end of x before the mean is computed.
#' @param slot.name the data in the slot.name to use
#' @param data.use a customed data matrix. Default: data.use = NULL and the expression matrix in the 'slot.name' is used
#'
#' @return Returns a matrix with genes as rows, cell groups as columns.

#' @export
#'
computeAveExpr <- function(object, features = NULL, group.by = NULL, type = c("triMean", "truncatedMean", "median"), trim = NULL,
                           slot.name = c("data.signaling", "data"), data.use = NULL) {
  type <- match.arg(type)
  slot.name <- match.arg(slot.name)
  FunMean <- switch(type,
                    triMean = triMean,
                    truncatedMean = function(x) mean(x, trim = trim, na.rm = TRUE),
                    median = function(x) median(x, na.rm = TRUE))

  if (is.null(data.use)) {
    data.use <- slot(object, slot.name)
  }
  if (is.null(features)) {
    features.use <- row.names(data.use)
  } else {
    features.use <- intersect(features, row.names(data.use))
  }
  data.use <- data.use[features.use, , drop = FALSE]
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
  # compute the average expression per group
  data.use.avg <- aggregate(t(data.use), list(labels), FUN = FunMean)
  data.use.avg <- t(data.use.avg[,-1])
  rownames(data.use.avg) <- features.use
  colnames(data.use.avg) <- levels(labels)
  return(data.use.avg)
}



#' Compute the expression of complex in individual cells using geometric mean
#' @param complex_input the complex_input from CellChatDB
#' @param data.use data matrix (row are genes and columns are cells or cell groups)
#' @param complex the names of complex
#' @return
#' @importFrom dplyr select starts_with
#' @importFrom future nbrOfWorkers
#' @importFrom future.apply future_sapply
#' @importFrom pbapply pbsapply
#' @export
computeExpr_complex <- function(complex_input, data.use, complex) {
  Rsubunits <- complex_input[complex,] %>% dplyr::select(starts_with("subunit"))
  my.sapply <- ifelse(
    test = future::nbrOfWorkers() == 1,
    yes = sapply,
    no = future.apply::future_sapply
  )
  data.complex = my.sapply(
    X = 1:nrow(Rsubunits),
    FUN = function(x) {
      RsubunitsV <- unlist(Rsubunits[x,], use.names = F)
      RsubunitsV <- RsubunitsV[RsubunitsV != ""]
      return(geometricMean(data.use[RsubunitsV,]))
    }
  )
  data.complex <- t(data.complex)
  return(data.complex)
}

# Compute the average expression of complex per cell group using geometric mean
# @param complex_input the complex_input from CellChatDB
# @param data.use data matrix (rows are genes and columns are cells)
# @param complex the names of complex
# @param group a factor defining the cell groups
# @param FunMean the function for computing mean expression per group
# @return
# @importFrom dplyr select starts_with
# @importFrom future nbrOfWorkers
# @importFrom future.apply future_sapply
# @importFrom pbapply pbsapply
# #' @export
.computeExprGroup_complex <- function(complex_input, data.use, complex, group, FunMean) {
  Rsubunits <- complex_input[complex,] %>% dplyr::select(starts_with("subunit"))
  my.sapply <- ifelse(
    test = future::nbrOfWorkers() == 1,
    yes = pbapply::pbsapply,
    no = future.apply::future_sapply
  )
  data.complex = my.sapply(
    X = 1:nrow(Rsubunits),
    FUN = function(x) {
      RsubunitsV <- unlist(Rsubunits[x,], use.names = F)
      RsubunitsV <- RsubunitsV[RsubunitsV != ""]
      RsubunitsV <- intersect(RsubunitsV, rownames(data.use))
      if (length(RsubunitsV) > 1) {
        data.avg <- aggregate(t(data.use[RsubunitsV,]), list(group), FUN = FunMean)
        data.avg <- t(data.avg[,-1])
      } else if (length(RsubunitsV) == 1) {
        data.avg <- aggregate(matrix(data.use[RsubunitsV,], ncol = 1), list(group), FUN = FunMean)
        data.avg <- t(data.avg[,-1])
      } else {
        data.avg = matrix(0, nrow = 1, ncol = length(unique(group)))
      }
      return(geometricMean(data.avg))
    }
  )
  data.complex <- t(data.complex)
  return(data.complex)
}

#' Compute the expression of ligands or receptors using geometric mean
#' @param geneLR a char vector giving a set of ligands or receptors
#' @param data.use data matrix (row are genes and columns are cells or cell groups)
#' @param complex_input the complex_input from CellChatDB
# #' @param group a factor defining the cell groups; If NULL, compute the expression of ligands or receptors in individual cells; otherwise, compute the average expression of ligands or receptors per cell group
# #' @param FunMean the function for computing average expression per cell group
#' @return
#' @export
computeExpr_LR <- function(geneLR, data.use, complex_input){
  nLR <- length(geneLR)
  numCluster <- ncol(data.use)
  index.singleL <- which(geneLR %in% rownames(data.use))
  dataL1avg <- data.use[geneLR[index.singleL],]
  dataLavg <- matrix(nrow = nLR, ncol = numCluster)
  dataLavg[index.singleL,] <- dataL1avg
  index.complexL <- setdiff(1:nLR, index.singleL)
  if (length(index.complexL) > 0) {
    complex <- geneLR[index.complexL]
    data.complex <- computeExpr_complex(complex_input, data.use, complex)
    dataLavg[index.complexL,] <- data.complex
  }
  return(dataLavg)
}


#' Modeling the effect of coreceptor on the ligand-receptor interaction
#'
#' @param data.use data matrix
#' @param cofactor_input the cofactor_input from CellChatDB
#' @param pairLRsig a data frame giving ligand-receptor interactions
#' @param type when type == "A", computing expression of co-activation receptor; when type == "I", computing expression of co-inhibition receptor.
#' @return
#' @importFrom future nbrOfWorkers
#' @importFrom future.apply future_sapply
#' @importFrom pbapply pbsapply
#' @export
computeExpr_coreceptor <- function(cofactor_input, data.use, pairLRsig, type = c("A", "I")) {
  type <- match.arg(type)
  if (type == "A") {
    coreceptor.all = pairLRsig$co_A_receptor
  } else if (type == "I"){
    coreceptor.all = pairLRsig$co_I_receptor
  }
  index.coreceptor <- which(!is.na(coreceptor.all) & coreceptor.all != "")
  if (length(index.coreceptor) > 0) {
    my.sapply <- ifelse(
      test = future::nbrOfWorkers() == 1,
      yes = sapply,
      no = future.apply::future_sapply
    )
    coreceptor <- coreceptor.all[index.coreceptor]
    coreceptor.ind <- cofactor_input[coreceptor, grepl("cofactor" , colnames(cofactor_input) )]
    data.coreceptor.ind = my.sapply(
      X = 1:nrow(coreceptor.ind),
      FUN = function(x) {
        coreceptor.indV <- unlist(coreceptor.ind[x,], use.names = F)
        coreceptor.indV <- coreceptor.indV[coreceptor.indV != ""]
        coreceptor.indV <- intersect(coreceptor.indV, rownames(data.use))
        if (length(coreceptor.indV) == 1) {
          return(1 + data.use[coreceptor.indV, ])
        } else if (length(coreceptor.indV) > 1) {
          return(apply(1 + data.use[coreceptor.indV, ], 2, prod))
        } else {
          return(matrix(1, nrow = 1, ncol = ncol(data.use)))
        }
      }
    )
    data.coreceptor.ind <- t(data.coreceptor.ind)
    data.coreceptor <- matrix(1, nrow = length(coreceptor.all), ncol = ncol(data.use))
    data.coreceptor[index.coreceptor,] <- data.coreceptor.ind
  } else {
    data.coreceptor <- matrix(1, nrow = length(coreceptor.all), ncol = ncol(data.use))
  }
  return(data.coreceptor)
}

# Modeling the effect of coreceptor on the ligand-receptor interaction
#
# @param data.use data matrix
# @param cofactor_input the cofactor_input from CellChatDB
# @param pairLRsig a data frame giving ligand-receptor interactions
# @param type when type == "A", computing expression of co-activation receptor; when type == "I", computing expression of co-inhibition receptor.
# @param group a factor defining the cell groups
# @param FunMean the function for computing mean expression per group
# @return
# @importFrom future nbrOfWorkers
# @importFrom future.apply future_sapply
# @importFrom pbapply pbsapply
# #' @export
.computeExprGroup_coreceptor <- function(cofactor_input, data.use, pairLRsig, type = c("A", "I"), group, FunMean) {
  type <- match.arg(type)
  if (type == "A") {
    coreceptor.all = pairLRsig$co_A_receptor
  } else if (type == "I"){
    coreceptor.all = pairLRsig$co_I_receptor
  }
  index.coreceptor <- which(!is.na(coreceptor.all) & coreceptor.all != "")
  if (length(index.coreceptor) > 0) {
    my.sapply <- ifelse(
      test = future::nbrOfWorkers() == 1,
      yes = pbapply::pbsapply,
      no = future.apply::future_sapply
    )
    coreceptor <- coreceptor.all[index.coreceptor]
    coreceptor.ind <- cofactor_input[coreceptor, grepl("cofactor" , colnames(cofactor_input) )]
    data.coreceptor.ind = my.sapply(
      X = 1:nrow(coreceptor.ind),
      FUN = function(x) {
        coreceptor.indV <- unlist(coreceptor.ind[x,], use.names = F)
        coreceptor.indV <- coreceptor.indV[coreceptor.indV != ""]
        coreceptor.indV <- intersect(coreceptor.indV, rownames(data.use))
        if (length(coreceptor.indV) > 1) {
          data.avg <- aggregate(t(data.use[coreceptor.indV,]), list(group), FUN = FunMean)
          data.avg <- t(data.avg[,-1])
          return(apply(1 + data.avg, 2, prod))
          # return(1 + apply(data.avg, 2, mean))
        } else if (length(coreceptor.indV) == 1) {
          data.avg <- aggregate(matrix(data.use[coreceptor.indV,], ncol = 1), list(group), FUN = FunMean)
          data.avg <- t(data.avg[,-1])
          return(1 + data.avg)
        } else {
          return(matrix(1, nrow = 1, ncol = length(unique(group))))
        }
      }
    )
    data.coreceptor.ind <- t(data.coreceptor.ind)
    data.coreceptor <- matrix(1, nrow = length(coreceptor.all), ncol = length(unique(group)))
    data.coreceptor[index.coreceptor,] <- data.coreceptor.ind
  } else {
    data.coreceptor <- matrix(1, nrow = length(coreceptor.all), ncol = length(unique(group)))
  }

  return(data.coreceptor)
}

#' Modeling the effect of agonist on the ligand-receptor interaction
#' @param data.use data matrix
#' @param cofactor_input the cofactor_input from CellChatDB
#' @param pairLRsig the L-R interactions
#' @param group a factor defining the cell groups
#' @param index.agonist the index of agonist in the database
#' @param Kh a parameter in Hill function
#' @param FunMean the function for computing mean expression per group
#' @param n Hill coefficient
#' @return
#' @export
#' @importFrom stats aggregate
computeExprGroup_agonist <- function(data.use, pairLRsig, cofactor_input, group, index.agonist, Kh, FunMean, n) {
  agonist <- pairLRsig$agonist[index.agonist]
  agonist.ind <- cofactor_input[agonist, grepl("cofactor" , colnames(cofactor_input))]
  agonist.indV <- unlist(agonist.ind, use.names = F)
  agonist.indV <- agonist.indV[agonist.indV != ""]
  agonist.indV <- intersect(agonist.indV, rownames(data.use))
  if (length(agonist.indV) == 1) {
    data.avg <- aggregate(matrix(data.use[agonist.indV,], ncol = 1), list(group), FUN = FunMean)
    data.avg <- t(data.avg[,-1])
    data.agonist <- 1 + data.avg^n/(Kh^n + data.avg^n)
  } else if (length(agonist.indV) > 1) {
    data.avg <- aggregate(t(data.use[agonist.indV,]), list(group), FUN = FunMean)
    data.avg <- t(data.avg[,-1])
    data.agonist <- apply(1 + data.avg^n/(Kh^n + data.avg^n), 2, prod)
  } else {
    data.agonist = matrix(1, nrow = 1, ncol = length(unique(group)))
  }
  return(data.agonist)
}

#' Modeling the effect of antagonist on the ligand-receptor interaction
#'
#' @param data.use data matrix
#' @param cofactor_input the cofactor_input from CellChatDB
#' @param pairLRsig the L-R interactions
#' @param group a factor defining the cell groups
#' @param index.antagonist the index of antagonist in the database
#' @param Kh a parameter in Hill function
#' @param n Hill coefficient
#' @param FunMean the function for computing mean expression per group
#' @return
#' @export
#' @importFrom stats aggregate
computeExprGroup_antagonist <- function(data.use, pairLRsig, cofactor_input, group, index.antagonist, Kh, FunMean, n) {
  antagonist <- pairLRsig$antagonist[index.antagonist]
  antagonist.ind <- cofactor_input[antagonist, grepl( "cofactor" , colnames(cofactor_input) )]
  antagonist.indV <- unlist(antagonist.ind, use.names = F)
  antagonist.indV <- antagonist.indV[antagonist.indV != ""]
  antagonist.indV <- intersect(antagonist.indV, rownames(data.use))
  if (length(antagonist.indV) == 1) {
    data.avg <- aggregate(matrix(data.use[antagonist.indV,], ncol = 1), list(group), FUN = FunMean)
    data.avg <- t(data.avg[,-1])
    data.antagonist <- Kh^n/(Kh^n + data.avg^n)
  } else if (length(antagonist.indV) > 1) {
    data.avg <- aggregate(t(data.use[antagonist.indV,]), list(group), FUN = FunMean)
    data.avg <- t(data.avg[,-1])
    data.antagonist <- apply(Kh^n/(Kh^n + data.avg^n), 2, prod)
  } else {
    data.antagonist = matrix(1, nrow = 1, ncol = length(unique(group)))
  }
  return(data.antagonist)
}


#' Modeling the effect of agonist on the ligand-receptor interaction
#' @param data.use data matrix
#' @param cofactor_input the cofactor_input from CellChatDB
#' @param pairLRsig the L-R interactions
# #' @param group a factor defining the cell groups
#' @param index.agonist the index of agonist in the database
#' @param Kh a parameter in Hill function
# #' @param FunMean the function for computing mean expression per group
#' @param n Hill coefficient
#' @return
#' @export
#' @importFrom stats aggregate
computeExpr_agonist <- function(data.use, pairLRsig, cofactor_input, index.agonist, Kh,  n) {
  agonist <- pairLRsig$agonist[index.agonist]
  agonist.ind <- cofactor_input[agonist, grepl("cofactor" , colnames(cofactor_input))]
  agonist.indV <- unlist(agonist.ind, use.names = F)
  agonist.indV <- agonist.indV[agonist.indV != ""]
  agonist.indV <- intersect(agonist.indV, rownames(data.use))
  if (length(agonist.indV) == 1) {
    # data.avg <- aggregate(matrix(data.use[agonist.indV,], ncol = 1), list(group), FUN = FunMean)
    # data.avg <- t(data.avg[,-1])
    data.avg <- data.use[agonist.indV,, drop = FALSE]
    data.agonist <- 1 + data.avg^n/(Kh^n + data.avg^n)
  } else if (length(agonist.indV) > 1) {
    # data.avg <- aggregate(t(data.use[agonist.indV,]), list(group), FUN = FunMean)
    # data.avg <- t(data.avg[,-1])
    data.avg <- data.use[agonist.indV,, drop = FALSE]
    data.agonist <- apply(1 + data.avg^n/(Kh^n + data.avg^n), 2, prod)
  } else {
    # data.agonist = matrix(1, nrow = 1, ncol = length(unique(group)))
    data.agonist = matrix(1, nrow = 1, ncol = ncol(data.use))
  }
  return(data.agonist)
}

#' Modeling the effect of antagonist on the ligand-receptor interaction
#'
#' @param data.use data matrix
#' @param cofactor_input the cofactor_input from CellChatDB
#' @param pairLRsig the L-R interactions
# #' @param group a factor defining the cell groups
#' @param index.antagonist the index of antagonist in the database
#' @param Kh a parameter in Hill function
#' @param n Hill coefficient
# #' @param FunMean the function for computing mean expression per group
#' @return
#' @export
#' @importFrom stats aggregate
computeExpr_antagonist <- function(data.use, pairLRsig, cofactor_input, index.antagonist, Kh, n) {
  antagonist <- pairLRsig$antagonist[index.antagonist]
  antagonist.ind <- cofactor_input[antagonist, grepl( "cofactor" , colnames(cofactor_input) )]
  antagonist.indV <- unlist(antagonist.ind, use.names = F)
  antagonist.indV <- antagonist.indV[antagonist.indV != ""]
  antagonist.indV <- intersect(antagonist.indV, rownames(data.use))
  if (length(antagonist.indV) == 1) {
    # data.avg <- aggregate(matrix(data.use[antagonist.indV,], ncol = 1), list(group), FUN = FunMean)
    # data.avg <- t(data.avg[,-1])
    data.avg <- data.use[antagonist.indV,, drop = FALSE]
    data.antagonist <- Kh^n/(Kh^n + data.avg^n)
  } else if (length(antagonist.indV) > 1) {
    # data.avg <- aggregate(t(data.use[antagonist.indV,]), list(group), FUN = FunMean)
    # data.avg <- t(data.avg[,-1])
    data.avg <- data.use[antagonist.indV,, drop = FALSE]
    data.antagonist <- apply(Kh^n/(Kh^n + data.avg^n), 2, prod)
  } else {
    # data.antagonist = matrix(1, nrow = 1, ncol = length(unique(group)))
    data.antagonist = matrix(1, nrow = 1, ncol = ncol(data.use))
  }
  return(data.antagonist)
}


#' Compute the geometric mean
#' @param x a numeric vector
#' @param na.rm whether remove na
#' @return
#' @export
geometricMean <- function(x,na.rm=TRUE){
  if (is.null(nrow(x))) {
    exp(mean(log(x),na.rm=na.rm))
  } else {
    exp(apply(log(x),2,mean,na.rm=na.rm))
  }
}


#' Compute the Tukey's trimean
#' @param x a numeric vector
#' @param na.rm whether remove na
#' @return
#' @importFrom stats quantile
#' @export
triMean <- function(x, na.rm = TRUE) {
  mean(stats::quantile(x, probs = c(0.25, 0.50, 0.50, 0.75), na.rm = na.rm))
}

#' Compute the average expression per cell group when the percent of expressing cells per cell group larger than a threshold
#' @param x a numeric vector
#' @param trim the percent of expressing cells per cell group to be considered as zero
#' @param na.rm whether remove na
#' @return
#' @importFrom Matrix nnzero
# #' @export
thresholdedMean <- function(x, trim = 0.1, na.rm = TRUE) {
  percent <- Matrix::nnzero(x)/length(x)
  if (percent < trim) {
    return(0)
  } else {
    return(mean(x, na.rm = na.rm))
  }
}

#' Identify all the significant interactions (L-R pairs) from some cell groups to other cell groups
#'
#' @param object CellChat object
#' @param from a vector giving the index or the name of source cell groups
#' @param to a corresponding vector giving the index or the name of target cell groups. Note: The length of 'from' and 'to' must be the same, giving the corresponding pair of cell groups for communication.
#' @param bidirection whether show the bidirectional communication, i.e., both 'from'->'to' and 'to'->'from'.
#' @param pair.only whether only return ligand-receptor pairs without pathway names and communication strength
#' @param pairLR.use0 ligand-receptor pairs to use; default is all the significant interactions
#' @param thresh threshold of the p-value for determining significant interaction
#'
#' @return
#' @export
#'
identifyEnrichedInteractions <- function(object, from, to, bidirection = FALSE, pair.only = TRUE, pairLR.use0 = NULL, thresh = 0.05){
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
    from2 <- c(from, to)
    to <- c(to, from)
    from <- from2
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
    prob_ij[pval_ij ==1] <- 0
    pval[,k] <- pval_ij
    prob[,k] <- prob_ij
    group.names <- c(group.names, paste(group.names.all[from[i]], group.names.all[to[i]], sep = " - "))
  }
  prob[which(prob == 0)] <- NA
  # remove rows that are entirely NA
  pval <- pval[rowSums(is.na(prob)) != ncol(prob), ,drop = FALSE]
  pairLR.use0 <- pairLR.use0[rowSums(is.na(prob)) != ncol(prob), ,drop = FALSE]
  prob <- prob[rowSums(is.na(prob)) != ncol(prob), ,drop = FALSE]
  if (pair.only) {
    pairLR.use0 <- dplyr::select(pairLR.use0, ligand, receptor)
  }
  return(pairLR.use0)
}

