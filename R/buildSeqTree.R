#' Build a phylogenetic tree
#'
#' @name buildSeqTree
#'
#' @details A wrapper to build phylogenetic tree for ASV sequences
#'          obtained after \code{dada2}.
#'
#' @param x Sequences DNAStringSet
#'
#' @param pml_model See optim.pml \code{phangorn}
#'
#' @param pml_optInv See optim.pml \code{phangorn}
#'
#' @param pml_optGamma See optim.pml \code{phangorn}
#'
#' @param pml_rearrangement See optim.pml \code{phangorn}
#' #'
#' @param pml_control See optim.pml \code{phangorn}
#'
#' @examples
#' library(biomeUtils)
#' data("SprockettTHData")
#' # select few for example
#' seqs <- sample(SprockettTHData@refseq, 5)
#' sp.tree <- buildSeqTree(seqs)
#' sp.tree$tree
#'
#' @importFrom DECIPHER AlignSeqs
#' @importFrom Biostrings DNAStringSet
#' @importFrom phangorn phyDat dist.ml NJ pml optim.pml pml.control
#' @importFrom stats update
#'
#' @export

buildSeqTree <- function(x,
                         pml_model="GTR",
                         pml_optInv=TRUE,
                         pml_optGamma=TRUE,
                         pml_rearrangement = "stochastic",
                         pml_control = pml.control(trace = 0)) {

  # global vars
  alignment <- phangAlign <- dm <- treeNJ <- fit <- fitGTR <- NULL

  if(class(x)=="DNAStringSet"){
    x <- x
  } else {
    x <- DNAStringSet(x)
  }

  alignment <- DECIPHER::AlignSeqs(x, anchor=NA, verbose=T)

  phangAlign <- phangorn::phyDat(as(alignment, "matrix"), type="DNA")

  dm <- phangorn::dist.ml(phangAlign)

  # Note, tip order != sequence order
  treeNJ <- phangorn::NJ(dm)

  fit <- phangorn::pml(treeNJ, data=phangAlign)

  fitGTR <- stats::update(fit, k=4, inv=0.2)

  fitGTR <- phangorn::optim.pml(fitGTR,
                                model = pml_model,
                                optInv = pml_optInv,
                                optGamma = pml_optGamma,
                                rearrangement = pml_rearrangement,
                                control = pml_control)

  return(fitGTR)

}
