#' Subset the ligand-receptor interactions for given specific signals in CellChatDB
#'
#' @param signaling a character vector
#' @param pairLR.use a dataframe containing ligand-receptor interactions
#' @param key the keyword to match
#' @param matching.exact whether perform exact matching
#' @param pair.only whether only return ligand-receptor pairs without cofactors
#' @importFrom future.apply future_sapply
#' @importFrom dplyr select
#' @return
#' @export
searchPair <- function(signaling = c(), pairLR.use, key = c("pathway_name","ligand"), matching.exact = FALSE, pair.only = TRUE) {
  key <- match.arg(key)
  pairLR = future.apply::future_sapply(
    X = 1:length(signaling),
    FUN = function(x) {
      if (!matching.exact) {
        index <- grep(signaling[x], pairLR.use[[key]])
      } else {
        index <- which(pairLR.use[[key]] %in% signaling[x])
      }
      if (length(index) > 0) {
        if (pair.only) {
          pairLR <- dplyr::select(pairLR.use[index, ], interaction_name, pathway_name, ligand, receptor)
        } else {
          pairLR <- pairLR.use[index, ]
        }
        return(pairLR)
      } else {
        stop(cat(paste("Cannot find ", signaling[x], ".", "Please input a correct name!"),'\n'))
      }
    }
  )
  if (pair.only) {
    pairLR0 <- vector("list", length(signaling))
    for (i in 1:length(signaling)) {
      pairLR0[[i]] <- matrix(unlist(pairLR[c(4*i-3, 4*i-2, 4*i-1, 4*i)]), ncol=4, byrow=F)
    }
    pairLR <- do.call(rbind, pairLR0)
    dimnames(pairLR)[[2]] <- dimnames(pairLR.use)[[2]][1:4]
    rownames(pairLR) <- pairLR[,1]
  } else {
    pairLR0 <- vector("list", length(signaling))
    for (i in 1:length(signaling)) {
      pairLR0[[i]] <- matrix(unlist(pairLR[(i*ncol(pairLR.use)-(ncol(pairLR.use)-1)):(i*ncol(pairLR.use))]), ncol=ncol(pairLR.use), byrow=F)
    }
    pairLR <- do.call(rbind, pairLR0)
    dimnames(pairLR)[[2]] <- dimnames(pairLR.use)[[2]]
    rownames(pairLR) <- pairLR[,1]
  }
  return(as.data.frame(pairLR, stringsAsFactors = FALSE))
}

#' Subset CellChatDB databse by only including interactions of interest
#'
#' @param CellChatDB CellChatDB databse
#' @param search a character
#' @param key the name of the variable in CellChatDB interaction_input
#'
#' @return
#' @export
#'
subsetDB <- function(CellChatDB, search = c(), key = "annotation") {
  interaction_input <- CellChatDB$interaction
  interaction_input <- interaction_input[interaction_input[[key]] == search, ]
  CellChatDB$interaction <- interaction_input
  return(CellChatDB)
}



#' Extract the genes involved in CellChatDB
#'
#' @param CellChatDB CellChatDB databse used in the analysis
#'
#' @return
#' @export
#' @importFrom dplyr select
#'
extractGene <- function(CellChatDB) {
  interaction_input <- CellChatDB$interaction
  complex_input <- CellChatDB$complex
  cofactor_input <- CellChatDB$cofactor
  geneIfo <- CellChatDB$geneInfo
  # check whether all gene names in complex_input and cofactor_input are official gene symbol in geneIfo
  checkGeneSymbol(geneSet = unlist(complex_input), geneIfo)
  checkGeneSymbol(geneSet = unlist(cofactor_input), geneIfo)

  geneL <- unique(interaction_input$ligand)
  geneR <- unique(interaction_input$receptor)
  geneLR <- c(geneL, geneR)
  checkGeneSymbol(geneSet = geneLR[geneLR %in% rownames(complex_input) == "FALSE"], geneIfo)

  geneL <- extractGeneSubset(geneL, complex_input, geneIfo)
  geneR <- extractGeneSubset(geneR, complex_input, geneIfo)
  geneLR <- c(geneL, geneR)

  cofactor <- c(interaction_input$agonist, interaction_input$antagonist, interaction_input$co_A_receptor, interaction_input$co_I_receptor)
  cofactor <- unique(cofactor[cofactor != ""])
  cofactorsubunits <- select(cofactor_input[match(cofactor, rownames(cofactor_input), nomatch=0),], starts_with("cofactor"))
  cofactorsubunitsV <- unlist(cofactorsubunits)
  geneCofactor <- unique(cofactorsubunitsV[cofactorsubunitsV != ""])

  gene.use <- unique(c(geneLR, geneCofactor))
  return(gene.use)

}


#' Extract the gene name
#'
#' @param geneSet gene set
#' @param complex_input complex in CellChatDB databse
#' @param geneIfo official gene symbol
#'
#' @return
#' @importFrom dplyr select starts_with
#' @export
extractGeneSubset <- function(geneSet, complex_input, geneIfo) {
  complex <- geneSet[which(geneSet %in% geneIfo$Symbol == "FALSE")]
  geneSet <- intersect(geneSet, geneIfo$Symbol)
  complexsubunits <- dplyr::select(complex_input[match(complex, rownames(complex_input), nomatch=0),], starts_with("subunit"))
  complex <- intersect(complex, rownames(complexsubunits))
  complexsubunitsV <- unlist(complexsubunits)
  complexsubunitsV <- unique(complexsubunitsV[complexsubunitsV != ""])
  geneSet <- unique(c(geneSet, complexsubunitsV))
  return(geneSet)
}


#' check the official Gene Symbol
#'
#' @param geneSet gene set to check
#' @param geneIfo official Gene Symbol
#' @return
#' @export
#'
checkGeneSymbol <- function(geneSet, geneIfo) {
  geneSet <- unique(geneSet[geneSet != ""])
  genes_notOfficial <- geneSet[geneSet %in% geneIfo$Symbol == "FALSE"]
  if (length(genes_notOfficial) > 0) {
    cat("Issue identified!! Please check the official Gene Symbol of the following genes: ", "\n", genes_notOfficial, "\n")
  }
  return(FALSE)
}
