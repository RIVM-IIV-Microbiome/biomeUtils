#' Find the Top Taxa
#'
#' @name findTopTaxa
#'
#' @param x A phyloseq object
#'
#' @param top Numeric value, how many top taxa to return. Default return top
#'   five taxa.
#'
#' @param method Specify the method to determine top taxa. Either sum, mean,
#'   median or prevalence. Default is 'mean'.
#'
#' @return A vector of the most \code{top} abundant taxa
#'
#' @examples
#' library(biomeUtils)
#' data("FuentesIliGutData")
#' findTopTaxa(FuentesIliGutData, top= 5L, method="prevalence")
#'
#' @author Sudarshan A. Shetty
#'
#' @references
#' Shetty SA (2020). Utilities for microbiome analytics.
#' \url{https://github.com/microsud/biomeUtils}
#'
#' @importFrom microbiome prevalence abundances
#' @importFrom utils head
#' @importFrom Biobase rowMedians
#'
#' @export

findTopTaxa <- function(x,
                        top= 5L,
                        method = c("mean","sum","median", "prevalence")){

  method <- match.arg(method, c("mean","sum","median","prevalence"))

  if(method == "prevalence"){

    taxs <- microbiome::prevalence(x)


  } else {

    taxs <- switch(method,
                   mean = rowMeans(microbiome::abundances(x),
                                   na.rm = T),
                   sum = rowSums(microbiome::abundances(x),
                                 na.rm = T),
                   median = Biobase::rowMedians(microbiome::abundances(x),
                                                na.rm = T))
    names(taxs) <- taxa_names(x)

  }

  taxs <- sort(taxs,decreasing = TRUE)
  head(names(taxs), n = top)

}
