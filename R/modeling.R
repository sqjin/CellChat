
#' Compute the communication probability/strength of any interacting cell groups
#'
#' @param object CellChat object
#' @param LR.use ligand-receptor interactions used in inferring communication network
#' @param trim the fraction (0 to 0.25) of observations to be trimmed from each end of x before the mean is computed; default using the robust Tukey's trimean
#' @param nboot the total number of permutations
#' @param seed.use set a random seed. By default, sets the seed to 1.
#' @param Kh parameter in Hill function
#' @importFrom future nbrOfWorkers
#' @importFrom future.apply future_sapply
#' @importFrom pbapply pbsapply
#' @importFrom stats aggregate
#' @importFrom Matrix crossprod
#'
#' @return
#' @export
#'
#' @examples
computeCommunProb <- function(object, LR.use = NULL, trim = NULL, nboot = 100, seed.use = 1L, Kh = 0.5) {
  data <- object@data.project
  if (is.null(LR.use)) {
    pairLR.use <- object@LR$LRsig
  } else {
    pairLR.use <- LR.use
  }
  complex_input <- object@DB$complex
  cofactor_input <- object@DB$cofactor
  my.sapply <- ifelse(
    test = future::nbrOfWorkers() == 1,
    yes = pbapply::pbsapply,
    no = future.apply::future_sapply
  )
  ptm = Sys.time()

  pairLRsig <- pairLR.use
  group <- object@idents
  geneL <- as.character(pairLRsig$ligand)
  geneR <- as.character(pairLRsig$receptor)
  nLR <- nrow(pairLRsig)
  numCluster <- nlevels(group)
  data.use <- data/max(data)
  nC <- ncol(data.use)

  # compute the expression of ligand
  index.singleL <- which(geneL %in% rownames(data.use))
  dataL1 <- data.use[geneL[index.singleL],]
  dataL <- matrix(nrow = nLR, ncol = nC)

  dataL[index.singleL,] <- dataL1
  index.complexL <- setdiff(1:nLR, index.singleL)
  if (length(index.complexL) > 0) {
    complex <- geneL[index.complexL]
    data.complex <- computeExpr_complex(complex_input, data.use, complex)
    dataL[index.complexL,] <- data.complex
  }

  # compute the expression of receptor
  index.singleR <- which(geneR %in% rownames(data.use))
  dataR1 <- data.use[geneR[index.singleR],]
  dataR <- matrix(nrow = nLR, ncol = nC)

  dataR[index.singleR,] <- dataR1
  index.complexR <- setdiff(1:nLR, index.singleR)
  if (length(index.complexR) > 0) {
    complex <- geneR[index.complexR]
    data.complex <- computeExpr_complex(complex_input, data.use, complex)
    dataR[index.complexR,] <- data.complex
  }

  # take account into the effect of co-activation and co-inhibition receptors
  index.co.A.receptor <- which(!is.na(pairLRsig$co_A_receptor) & pairLRsig$co_A_receptor != "")
  if (length(index.co.A.receptor) > 0) {
    dataR.co.A.receptor <- computeExpr_coreceptor(cofactor_input, data.use, coreceptor.all = pairLRsig$co_A_receptor, index.coreceptor = index.co.A.receptor)
  } else{
    dataR.co.A.receptor <- matrix(1, nrow = nrow(dataR), ncol = ncol(dataR))
  }

  index.co.I.receptor <- which(!is.na(pairLRsig$co_I_receptor) & pairLRsig$co_I_receptor != "")
  if (length(index.co.I.receptor) > 0) {
    dataR.co.I.receptor <- computeExpr_coreceptor(cofactor_input, data.use, coreceptor.all = pairLRsig$co_I_receptor, index.coreceptor = index.co.I.receptor)
  } else {
    dataR.co.I.receptor <- matrix(1, nrow = nrow(dataR), ncol = ncol(dataR))
  }
  dataR <- dataR * dataR.co.A.receptor/dataR.co.I.receptor
  if (is.null(trim)) {
    FunMean <- triMean
  } else {
    FunMean <- function(x) mean(x, trim = trim, na.rm = TRUE)
  }
  # compute the average expression in each cell group
  dataLavg <- aggregate(t(dataL), list(group), FUN = FunMean)
  dataLavg <- t(dataLavg[,-1])
  dataRavg <- aggregate(t(dataR), list(group), FUN = FunMean)
  dataRavg <- t(dataRavg[,-1])

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
  for (i in 1:nLR) {
    # ligand/receptor
    dataLR <- Matrix::crossprod(matrix(dataLavg[i,], nrow = 1), matrix(dataRavg[i,], nrow = 1))
    P1 <- dataLR/(Kh + dataLR)
    if (sum(P1) == 0) {
      Pnull = P1
      Prob[ , , i] <- Pnull
      p = 1
      Pval[, , i] <- matrix(p, nrow = numCluster, ncol = numCluster, byrow = FALSE)
    } else {
      # agonist and antagonist
      if (is.element(i, index.agonist)) {
        data.agonist <- computeExpr_agonist(data.use = data.use, pairLRsig, cofactor_input, group = group,index.agonist = i, Kh = Kh, FunMean = FunMean)
        P2 <- Matrix::crossprod(matrix(data.agonist, nrow = 1))
      } else {
        P2 <- matrix(1, nrow = numCluster, ncol = numCluster)
      }
      if (is.element(i, index.antagonist)) {
        data.antagonist <- computeExpr_antagonist(data.use = data.use, pairLRsig, cofactor_input, group = group, index.antagonist = i, Kh = Kh, FunMean = FunMean)
        P3 <- Matrix::crossprod(matrix(data.antagonist, nrow = 1))
      } else {
        P3 <- matrix(1, nrow = numCluster, ncol = numCluster)
      }
      # number of cells
      P4 <- Matrix::crossprod(matrix(dataLavg2[i,], nrow = 1), matrix(dataRavg2[i,], nrow = 1))

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
          P1.boot <- dataLRB/(Kh + dataLRB)
          # agonist and antagonist
          if (is.element(i, index.agonist)) {
            data.agonist <- computeExpr_agonist(data.use = data.use, pairLRsig, cofactor_input, group = groupboot, index.agonist = i, Kh = Kh, FunMean = FunMean)
            P2.boot <- Matrix::crossprod(matrix(data.agonist, nrow = 1))
          } else {
            P2.boot <- matrix(1, nrow = numCluster, ncol = numCluster)
          }
          if (is.element(i, index.antagonist)) {
            data.antagonist <- computeExpr_antagonist(data.use = data.use, pairLRsig, cofactor_input, group = groupboot,index.antagonist = i, Kh = Kh, FunMean = FunMean)
            P3.boot <- Matrix::crossprod(matrix(data.antagonist, nrow = 1))
          } else {
            P3.boot <- matrix(1, nrow = numCluster, ncol = numCluster)
          }
          dataLavg2B <- by(matrix(dataL2.i, ncol = 1), groupboot, sum)/nC
          dataLavg2B <- matrix(dataLavg2B, nrow = 1)

          dataRavg2B <- by(matrix(dataR2.i, ncol = 1), groupboot, sum)/nC
          dataRavg2B <- matrix(dataRavg2B, nrow = 1)
          P4.boot = Matrix::crossprod(dataLavg2B, dataRavg2B)

          Pboot = P1.boot*P2.boot*P3.boot*P4.boot
          return(as.vector(Pboot))
        }
      )
      Pboot <- matrix(unlist(Pboot), nrow=length(Pnull), ncol = nboot, byrow = FALSE)
      nReject <- rowSums(Pboot - Pnull >= 0)
      p = nReject/nboot
      Pval[, , i] <- matrix(p, nrow = numCluster, ncol = numCluster, byrow = FALSE)
    }

  }
  Pval[Prob == 0] <- 1
  dimnames(Prob) <- list(levels(group), levels(group), rownames(pairLRsig))
  dimnames(Pval) <- dimnames(Prob)
  net <- list("prob" = Prob, "pval" = Pval)
  execution.time = Sys.time() - ptm
  object@options$run.time <- as.numeric(execution.time, units = "secs")
  object@net <- net
  return(object)
}


#' Compute the communication probability on signaling pathway level by summarizing all related ligands/receptors
#'
#' @param object CellChat object
#' @param thresh threshold of the p-value for determining significant interaction
#'
#' @return
#' @export
#'
#' @examples
computeCommunProbPathway <- function(object, thresh = 0.05) {
  net <- object@net
  pairLR.use <- object@LR$LRsig
  prob <- net$prob
  prob[net$pval > thresh] <- 0
  pathways <- unique(pairLR.use$pathway_name)
  group <- factor(pairLR.use$pathway_name, levels = pathways)
  prob.pathways <- aperm(apply(prob, c(1, 2), by, group, sum), c(2, 3, 1))
  pathways.sig <- pathways[apply(prob.pathways, 3, sum) != 0]
  prob.pathways.sig <- prob.pathways[,,pathways.sig]
  object@netP$pathways <- pathways.sig
  object@netP$prob <- prob.pathways.sig
  return(object)
}


#' Calculate the aggregated network by counting the number of links or summarizing the communication probability
#'
#' @param object CellChat object
#' @param thresh threshold of the p-value for determining significant interaction
#'
#' @return
#' @export
#'
#' @examples
aggregateNet <- function(object, thresh = 0.05) {
  net <- object@net
  prob <- net$prob
  pval <- net$pval
  net$count <- apply(pval < thresh, c(1,2), sum)
  prob[pval >= thresh] <- 0
  net$sum <- apply(prob, c(1,2), sum)
  net$sum[is.na(net$sum)] <- 0
  net$count[is.na(net$count)] <- 0
  object@net <- net
  return(object)
}


#' Compute the expression of complex using geometric mean
#' @param complex_input the complex_input from CellChatDB
#' @param data.use data matrix
#' @param complex the names of complex
#' @return
#' @importFrom dplyr select starts_with
#' @importFrom future nbrOfWorkers
#' @importFrom future.apply future_sapply
#' @importFrom pbapply pbsapply
computeExpr_complex <- function(complex_input, data.use, complex) {
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
      return(geometricMean(data.use[RsubunitsV,]))
    }
  )
  data.complex <- t(data.complex)
  return(data.complex)
}


#' Modeling the effect of coreceptor on the ligand-receptor interaction
#'
#' @param data.use data matrix
#' @param cofactor_input the cofactor_input from CellChatDB
#' @param coreceptor.all all the coreceptor in the database
#' @param index.coreceptor the index of coreceptor in the database
#' @return
#' @importFrom future nbrOfWorkers
#' @importFrom future.apply future_sapply
#' @importFrom pbapply pbsapply
computeExpr_coreceptor <- function(cofactor_input, data.use, coreceptor.all, index.coreceptor) {
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
#' @return
#' @importFrom stats aggregate
computeExpr_agonist <- function(data.use, pairLRsig, cofactor_input, group, index.agonist, Kh, FunMean) {
  agonist <- pairLRsig$agonist[index.agonist]
  agonist.ind <- cofactor_input[agonist, grepl("cofactor" , colnames(cofactor_input))]
  agonist.indV <- unlist(agonist.ind, use.names = F)
  agonist.indV <- agonist.indV[agonist.indV != ""]
  agonist.indV <- intersect(agonist.indV, rownames(data.use))
  if (length(agonist.indV) == 1) {
    data.avg <- aggregate(matrix(data.use[agonist.indV,], ncol = 1), list(group), FUN = FunMean)
    data.avg <- t(data.avg[,-1])
    data.agonist <- 1 + data.avg/(Kh + data.avg)
  } else if (length(agonist.indV) > 1) {
    data.avg <- aggregate(t(data.use[agonist.indV,]), list(group), FUN = FunMean)
    data.avg <- t(data.avg[,-1])
    data.agonist <- apply(1 + data.avg/(Kh + data.avg), 2, prod)
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
#' @param FunMean the function for computing mean expression per group
#' @return
#' @importFrom stats aggregate
computeExpr_antagonist <- function(data.use, pairLRsig, cofactor_input, group, index.antagonist, Kh, FunMean) {
  antagonist <- pairLRsig$antagonist[index.antagonist]
  antagonist.ind <- cofactor_input[antagonist, grepl( "cofactor" , colnames(cofactor_input) )]
  antagonist.indV <- unlist(antagonist.ind, use.names = F)
  antagonist.indV <- antagonist.indV[antagonist.indV != ""]
  antagonist.indV <- intersect(antagonist.indV, rownames(data.use))
  if (length(antagonist.indV) == 1) {
    data.avg <- aggregate(matrix(data.use[antagonist.indV,], ncol = 1), list(group), FUN = FunMean)
    data.avg <- t(data.avg[,-1])
    data.antagonist <- Kh/(Kh + data.avg)
  } else if (length(antagonist.indV) > 1) {
    data.avg <- aggregate(t(data.use[antagonist.indV,]), list(group), FUN = FunMean)
    data.avg <- t(data.avg[,-1])
    data.antagonist <- apply(Kh/(Kh + data.avg), 2, prod)
  } else {
    data.antagonist = matrix(1, nrow = 1, ncol = length(unique(group)))
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

