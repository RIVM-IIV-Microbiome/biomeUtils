#' Get Overview of Taxa Abundances
#'
#' @name getTaxaOverview
#'
#' @param x A phyloseq object
#'
#' @param rank Default is NULL and returns overview matching \code{taxa_names}.
#'             Else user can provide one of the ranks available in 'rank_names'
#'             of phyloseq object.
#'
#' @return A \code{tibble}
#'
#' @examples
#' library(biomeUtils)
#' data("FuentesIliGutData")
#' getProportions(FuentesIliGutData) |>
#'    filterSampleData(ILI == "C") |>
#'    getTaxaOverview(rank="Phylum")
#'
#' @author Sudarshan A. Shetty
#'
#' @seealso getTaxaSummary
#'
#' @references
#' Shetty SA (2020). Utilities for microbiome analytics.
#' \url{https://github.com/RIVM-IIV-Microbiome/biomeUtils}
#'
#' @importFrom microbiome prevalence abundances aggregate_taxa
#' @importFrom phyloseq nsamples ntaxa taxa_sums taxa_sums
#' @importFrom phyloseq merge_phyloseq phy_tree taxa_names<-
#' @importFrom tibble tibble
#' @importFrom stats sd median
#'
#' @export
getTaxaOverview <- function(x, rank = NULL){

  mean_abund <- median_abund <- taxa <- min_abund <- max_abund <- NULL
  stdev_abund <- cv_abund <- dectected_samples <- prevalence <- NULL
  taxa.ov <- n <- NULL

  ## Must be a phyloseq object
  if ( !is(x, "phyloseq") ){
    stop("input must be an phyloseq object.")
  }

  ## the input must have few samples
  if ( nsamples(x) < 1 ){
    stop("input must have at least one sample")
  }

  if ( ntaxa(x) < 1 ) {
    stop("input must have at least one taxa")
  }

  if(!is.null(rank)){
    x <- aggregate_taxa(x, level=rank)
  }

  taxa.ov <- tibble::tibble(taxa = taxa_names(x),
                            # number of taxa per sample
                            #sum_abund = taxa_sums(x),
                            mean_abund = rowMeans(microbiome::abundances(x)),
                            median_abund = apply(microbiome::abundances(x), 1, median),
                            min_abund = apply(microbiome::abundances(x), 1, min),
                            max_abund = apply(microbiome::abundances(x), 1, max),
                            stdev_abund = apply(microbiome::abundances(x), 1, sd),
                            cv_abund = stdev_abund/mean_abund,
                            dectected_samples = rowSums(microbiome::abundances(x) != 0),
                            # percent_of_total = sum_abund/sum(sample_sums(x)) *100,
                            prevalence = microbiome::prevalence(x) *100) |>
    dplyr::arrange(dplyr::desc(mean_abund))

  return(taxa.ov)
}
