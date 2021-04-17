#' Calculate Faith's Phylogenetic Diversity
#'
#' @name calculatePD
#'
#' @details A wrapper around \code{picante} for \code{phyloseq} objects
#'          to calculate and add Phylogenetic Diversity and Species Richness
#'          to sample data.
#'
#' @param x A phyloseq object
#'
#' @param include_root Logical if root node be included. Depends whether the
#'                     input phyloseq had rooted tree (default = TRUE).
#'
#' @param justDF Logical if TURE returns only sample data as tibble with
#'               PD and SR values as columns.Default is FALSE
#'
#' @examples
#' library(biomeUtils)
#' data("FuentesIliGutData")
#' # reduce size for example
#' ps1 <- subset_samples(FuentesIliGutData, ILI == "C")
#' ps1 <- prune_taxa(taxa_sums(ps1) > 0, ps1)
#'
#' meta_tib <- calculatePD(ps1, justDF = TRUE)
#' # check
#' meta_tib[c(1, 2, 3), c("PD", "SR")]
#' @return Either a phyloseq object or data.frame with results of PD
#'         added to sample_data()
#' @author Sudarshan A. Shetty
#'
#' @references
#' \itemize{
#' \item{}{Kembel, S.W., Cowan, P.D., et al., 2010. Picante: R tools
#'        for integrating phylogenies and ecology.
#'        \emph{Bioinformatics}, 26(11), pp.1463-1464.
#'        \url{https://doi.org/10.1093/bioinformatics/btq166}
#' }
#' }
#' @importFrom picante pd
#' @importFrom phyloseq sample_data 'sample_data<-'
#' @importFrom microbiome abundances meta
#'
#' @export

calculatePD <- function(x, justDF = FALSE, include_root = TRUE) {
  if (class(x) != "phyloseq") {
    stop("Input is not an object of phyloseq class")
  }

  meta_tb <- meta(x)
  otu_tb <- as.data.frame(abundances(x))
  asv.tree <- x@phy_tree
  df.pd <- pd(t(otu_tb),
    asv.tree,
    include.root = include_root
  )
  # t(ou_table) transposes the table for use in picante and the
  # tre file comes from the first code chunck we used to read
  # tree file (see making a phyloseq object section).

  if (justDF) {
    meta_tb <- cbind(meta_tb, df.pd)
    return(meta_tb)
  } else {
    meta_tb <- cbind(meta_tb, df.pd)
    sample_data(x) <- sample_data(meta_tb)
    # sample_data(x)$Phylogenetic_Diversity <- df.pd$PD
    return(x)
  }
}
