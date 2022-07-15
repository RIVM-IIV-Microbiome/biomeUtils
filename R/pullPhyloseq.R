#' Pull Information from Phyloseq Object
#'
#' @name pullPhyloseq
#'
#' @details These are alternative to subset and prune functions in
#'          \code{phyloseq}.
#'
#' The \code{pullSampleAbundance} can be used to get abundances of all taxa.
#'
#' The \code{pullTaxaAbundance} can be used to get abundances of a taxa in all
#' samples.
#'
#' The \code{pullSampleVariable} can be used to get all variables within a
#' column in \code{sample_data}.
#'
#' The \code{pullTaxaRank} can be used to get all taxomic information at a
#' particular rank in \code{tax_table}.
#'
#' The naming is using pull because it uses the \code{dplyr::pull} function.
#'
#' @param x \code{phyloseq} object
#'
#' @param ... Option to pass on to \code{\link[dplyr]{pull}}. Only when using
#'                      \code{filterSampleData} or \code{filterTaxaData}
#'                      sample data or taxonomy table as tibble.
#' @examples
#'
#' library(biomeUtils)
#' data("FuentesIliGutData")
#' pullSampleAbundance(FuentesIliGutData, "sample_1", name= "FeatureID")
#'
#' pullTaxaAbundance(FuentesIliGutData, "ASV1", name = "SampleID")
#'
#' pullTaxaRank(FuentesIliGutData, "Genus")
#'
#' @return A vector of requested value
#'
#' @author Sudarshan A. Shetty
#'
#' @references
#' Shetty SA (2020). Utilities for microbiome analytics.
#' \url{https://github.com/RIVM-IIV-Microbiome/biomeUtils}
#'
#'
NULL

#' @rdname pullPhyloseq
#' @aliases pullSampleAbundance
#' @importFrom dplyr pull
#' @export
#'
pullSampleAbundance <- function(x, ...) {

  if (!is(x, "phyloseq")) {
    stop("Input is not an object of phyloseq class")
  }

  #rownames <- .check_col_sam_data(x,rownames)


  y <- x %>%
    getAbundanceTibble() %>%
    pull(...)

  return(y)

}

#' @rdname pullPhyloseq
#' @aliases pullTaxaAbundance
#' @importFrom dplyr pull
#' @export
#'
pullTaxaAbundance <- function(x, ...) {

  if (!is(x, "phyloseq")) {
    stop("Input is not an object of phyloseq class")
  }

  #rownames <- .check_col_sam_data(x,rownames)

  y <- x %>%
    getAbundanceTibble() %>%
    column_to_rownames(colnames(.)[1]) %>%
    t() %>% as.data.frame() %>%
    rownames_to_column("SampleID") %>%
    pull(...)


  return(y)

}


#' @rdname pullPhyloseq
#' @aliases pullTaxaRank
#' @importFrom dplyr pull
#' @export
#'
pullTaxaRank <- function(x, ...) {

  if (!is(x, "phyloseq")) {
    stop("Input is not an object of phyloseq class")
  }

  #rownames <- .check_col_sam_data(x,rownames)

  y <- x %>%
    getTaxaTibble() %>%
    pull(...)

  return(y)

}
