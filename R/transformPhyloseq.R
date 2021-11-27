#' Common Transformations
#'
#' @name transformation
#'
#' @details Common transformation in microbiome like proportions, clr.
#'
#' \code{getProportions} convert counts to relative abundance.
#'
#' \code{getCLR} convert counts to centred log transform 'clr'
#'
#' @param x \code{phyloseq} object
#'
#' @examples
#' library(biomeUtils)
#' data("FuentesIliGutData")
#' ps.rel <- getProportions(FuentesIliGutData)
#' ps.clr <- getCLR(FuentesIliGutData)
#'
#' @return Transformed phyloseq
#'
#' @author Sudarshan A. Shetty
#'
#' @references
#' Shetty SA (2020). Utilities for microbiome analytics.
#' \url{https://github.com/RIVM-IIV-Microbiome/biomeUtils}
#'
#'
NULL

#' @rdname transformation
#' @aliases getProportions
#' @export
getProportions <- function(x) {

  if (class(x) != "phyloseq") {
    stop("Input is not an object of phyloseq class")
  }
  x <- microbiome::abundances(x, "compositional")
  return(x)

}

#' @rdname transformation
#' @aliases getCLR
#' @export
getCLR <- function(x){
  if (class(x) != "phyloseq") {
    stop("Input is not an object of phyloseq class")
  }
  x <- microbiome::abundances(x, "clr")
  return(x)
}

