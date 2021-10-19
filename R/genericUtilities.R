#' Generic utilities
#'
#' @name genericUtilities
#'
#' @details Generic utilities that are simple wrappers for doing common
#'          tasks.
#'
#' \code{meltAbundance} convert abundance data to long format.
#'
#' \code{transposeAbundance} transpose OTU for use with \code{vegan} pkg.
#'
#' @param x \code{phyloseq} object
#'
#' @param ... Arguments other than 'cols' to pass to tidyr pivot_longer function
#'
#' @param transform If transformation is required `transposeAbundance`
#'
#' @examples
#'
#' library(biomeUtils)
#' data("FuentesIliGutData")
#' melt.abund <- meltAbundance(FuentesIliGutData,
#'                             names_to = "SampleID",
#'                             values_to="Abundance")
#' melt.abund
#'
#' t.otu <- transposeAbundance(FuentesIliGutData)
#' dim(t.otu)
#' @return long data format.
#'
#' @author Sudarshan A. Shetty
#'
#' @references
#' Shetty SA (2020). Utilities for microbiome analytics.
#' \url{https://github.com/microsud/biomeUtils}
#'
#'
NULL

#' @rdname genericUtilities
#' @aliases meltAbundance
#' @importFrom dplyr %>%
#' @export
#'
meltAbundance <- function(x, ...){

  ab.tib <- getAbundanceTibble(x) %>%
    tidyr::pivot_longer(-FeatureID, ...)

 return(ab.tib)
}

# For non-phyloseq tools such as vegan where taxa ar columns.
#' @rdname genericUtilities
#' @aliases transposeAbundance
#' @importFrom microbiome abundances
#' @export
transposeAbundance <- function(x, transform = "identity"){

  t.ab <- t(abundances(x, transform = transform))

  return(t.ab)

}


