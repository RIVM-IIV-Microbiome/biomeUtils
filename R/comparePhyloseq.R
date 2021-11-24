#' Compare Phyloseq Objects
#'
#' @name comparePhyloseq
#'
#' @details Provided a named list of different \code{phyloseq} objects
#'          `comparePhyloseq` returns a data.frame with basic values in
#'          each of the \code{phyloseq} objects.
#'
#' @param x A named list of \code{phyloseq} objects to compare
#'
#' @examples
#' library(biomeUtils)
#' data("FuentesIliGutData")
#' ps1 <- subset_samples(FuentesIliGutData, ILI == "C")
#' ps1 <- prune_taxa(taxa_sums(ps1) > 0, ps1)
#'
#' ps2 <- subset_samples(FuentesIliGutData, ILI == "L1")
#' ps2 <- prune_taxa(taxa_sums(ps2) > 0, ps2)
#'
#' ps3 <- subset_samples(FuentesIliGutData, ILI == "L2")
#' ps3 <- prune_taxa(taxa_sums(ps3) > 0, ps3)
#'
#' ps.list <- c("C" = ps1, "L1" = ps2, "L2" = ps3)
#'
#' comparePhyloseq(ps.list)
#' @return A tibble
#'
#' @author Sudarshan A. Shetty
#'
#' @references
#' Shetty SA (2021). Utilities for microbiome analytics.
#' \url{https://github.com/RIVM-IIV-Microbiome/biomeUtils}
#'
#' @importFrom phyloseq nsamples ntaxa sample_sums taxa_sums
#' @importFrom microbiome abundances
#' @importFrom tibble rownames_to_column
#'
#' @export

comparePhyloseq <- function(x) {
  ps_df <- NULL
  for (i in x) {
    # print(i)
    # print(names(ps.listi))
    df <- data.frame(
      ntaxa = ntaxa(i),
      nsample = nsamples(i),
      min.reads = min(sample_sums(i)),
      max.reads = max(sample_sums(i)),
      total.read = sum(sample_sums(i)),
      average.reads = round(sum(sample_sums(i)) / nsamples(i), 2),
      singletons = length(taxa_sums(i)[taxa_sums(i) <= 1]),
      sparsity = length(which(abundances(i) == 0)) / length(abundances(i))
    )
    # df$pseq <-
    # print(df)
    # rownames(df) <- names(ps.list$i)
    ps_df <- rbind(ps_df, df)
  }
  # add row names
  rownames(ps_df) <- names(x)
  ps_df <- ps_df %>% rownames_to_column("input")
  return(ps_df)
}
