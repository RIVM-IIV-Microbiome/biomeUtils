#' Filter Taxa Based on Total Data
#'
#' @name filterTaxaTotal
#'
#' @details Provided a count data \code{phyloseq}, \code{filterTaxaTotal} calculates
#'          the percent of taxon at the rank level is calculated. Those that are
#'          less than the percent_thres are removed. This function works with ASV
#'          level data as the ASVs are merged at specified level for calculation.
#'          Therefore, ASVs that belong to low/rare abundance at a specified rank
#'          are removed.
#'
#' @param x \code{phyloseq} object.
#'
#' @param rank Taxonomic rank to use like Phylum.
#'
#' @param percent_thres Percent cut-off to use. If value is 10, then all phyla with
#'        less that 10 percent of the total counts in all data are removed.
#'
#' @param verbose Logical. Prints a list of removed taxa. Default is TRUE.
#'
#' @examples
#' library(biomeUtils)
#' data('FuentesIliGutData')
#' # below we filter Family that are less than 2% of the total data
#' ps.filt <- filterTaxaTotal(FuentesIliGutData,
#'                            rank = 'Family',
#'                            percent_thres = 2,
#'                            verbose = TRUE)
#' ps.filt
#'
#' @return Filtered phyloseq
#'
#' @author Sudarshan A. Shetty
#'
#' @importFrom phyloseq ntaxa nsamples
#' @importFrom dplyr %>%
#' @export
NULL
filterTaxaTotal <- function(x, rank = "Phylum", percent_thres = 50, verbose = TRUE) {

  Percent <- NULL
  if (!is(x, "phyloseq")) {

    stop("Input is not an object of phyloseq class")
  }

  sum.tb <- getTaxaSummary(x, rank = rank)
  keep_tax <- sum.tb %>%
    dplyr::filter(Percent > percent_thres) %>%
    dplyr::pull(!!rank)

  nw.x <- filterTaxaData(x, get(rank) %in% keep_tax)

  lost <- sum(sample_sums(x)) - sum(sample_sums(nw.x))

  if (lost/sum(sample_sums(x)) * 100 > 20) {
    warning(paste0("!! More than 20% counts in total data !!"))
  }


  if (verbose) {

    remove_tax <- sum.tb %>%
      dplyr::filter(Percent < percent_thres) %>%
      dplyr::pull(!!rank)
    message(paste0("Following taxa were removed at level: ", rank))
    message(paste0(remove_tax, "\n"))
    message(paste0("No. of taxa removed: ", ntaxa(x) - ntaxa(nw.x)))
    message(paste0("No. of samples removed: ", nsamples(x) - nsamples(nw.x)))

  }

  return(nw.x)
}

