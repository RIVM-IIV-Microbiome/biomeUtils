#' Remove Zeros
#'
#' @name removeZeros
#'
#' @details In phyloseq based analysis, prune_samples or subset_samples
#'          does not remove taxa that have zero abundances in remaining
#'          samples. In such cases user may want to exclude these taxa.
#'
#' @param x A phyloseq object or a data.frame with only numerical or matrix
#'
#' @examples
#' library(biomeUtils)
#' data("FuentesIliGutData")
#' ps1 <- subset_samples(FuentesIliGutData, ILI != "L2")
#' ps1
#' ps2 <- removeZeros(ps1)
#' ps2
#' @return If input is phyloseq, returns phyloseq with taxa that
#'         have zeros in all samples.
#'         If input is matrix or data frame, returned object is
#'         without rows and cloumns that sum to zero.
#'
#' @author Sudarshan A. Shetty
#'
#' @importFrom phyloseq prune_taxa taxa_sums
#'
#' @export

removeZeros <- function(x) {
  if (class(x) == "phyloseq") {
    x <- prune_taxa(taxa_sums(x) > 0, x)
    x <- prune_samples(sample_sums(x) > 0, x)
    return(x)
  } else {
    x <- x[rowSums(x) > 0, ]
    x <- x[, colSums(x) > 0]
    return(x)
  }
}
