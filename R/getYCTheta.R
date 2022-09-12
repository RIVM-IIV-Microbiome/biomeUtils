#' Yue and Clayton Measure of Dissimilarity
#'
#' @name getYCTheta
#'
#' @details Calculates Yue and Clayton measure of dissimilarity between samples.
#'          It considers total number of taxa between two samples and the
#'          relative abundance of each taxa.
#'
#' @param x A phyloseq object or a data.frame with taxa-columns and samples-rows
#'
#' @examples
#' library(biomeUtils)
#' data("FuentesIliGutData")
#' yc.dis <- getYCTheta(FuentesIliGutData)
#' #yc.dis[1:3,1:4]
#' @return Object of class dist
#' @author Sudarshan A. Shetty
#'
#' @references
#' \itemize{
#' \item{}{Yue, Jack C., and Murray K. Clayton. 2010. A similarity measure
#'         based on species proportions.
#'        \emph{Communications in Statistics-theory and Methods}, 34(11), pp.2123-2131.
#'        \url{https://doi.org/10.1080/STA-200066418}
#' }
#' {Dixon, P 2003. VEGAN, a package of R functions for community ecology.
#'        \emph{Journal of Vegetation Science}, 14(6), pp.927-930.
#'        \url{https://doi.org/10.1111/j.1654-1103.2003.tb02228.x}
#' }
#' }
#' @importFrom vegan designdist
#' @importFrom microbiome abundances
#'
#' @export
NULL
getYCTheta <- function(x){

  yc.dist <- NULL
  if (is(x) == "phyloseq") {
    yc.dist <- vegan::designdist(t(microbiome::abundances(x)),
                                 method="1-(J/(A+B-J))",
                                 terms = c( "quadratic"),
                                 abcd = FALSE)
  } else {
    yc.dist <- vegan::designdist(x,
                                 method="1-(J/(A+B-J))",
                                 terms = c( "quadratic"),
                                 abcd = FALSE)
  }
  return(yc.dist)
}
