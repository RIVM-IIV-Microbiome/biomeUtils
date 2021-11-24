#' Check Filter Effect
#'
#' @name chekFilterEffect
#'
#' @details Investigate effect of abundance and prevalence thresholds on
#'          taxa numbers.
#'
#' \code{checkAbundanceCutOffs} check how many taxa are lost at different
#' abundance thresholds.
#'
#' \code{checkPrevalenceCutOffs} check how many taxa are lost at different
#' prevalence thresholds.
#'
#' @param x \code{phyloseq} object
#'
#' @param abundance_thres If using 'checkAbundanceCutOffs', can have multiple
#'                          numeric values like c(0.1, 0.01, 0.001, 0.0001).
#'                          If using 'checkPrevalenceCutOffs', a single value.
#'
#' @param prevalence_thres If using 'checkPrevalenceCutOffs', can have multiple
#'                          numeric values like c(0.05, 0.1, 0.15, 0.2).
#'                          If using 'checkAbundanceCutOffs', a single value.
#'
#' @examples
#' library(biomeUtils)
#' library(microbiome)
#' data("FuentesIliGutData")
#' ps.rel <- microbiome::transform(FuentesIliGutData, "compositional")
#' ab.check <- checkAbundanceCutOffs(ps.rel,
#'                                   abundance_thres = c(0.0001, 0.001, 0.01, 0.1),
#'                                   prevalence_thres = 0.1)
#' ab.check
#'
#' pv.check <- checkPrevalenceCutOffs(ps.rel,
#'                                    abundance_thres = 0.01,
#'                                    prevalence_thres = c(0.05, 0.1, 0.15, 0.2))
#' pv.check
#'
#' @return A tibble
#'
#' @author Sudarshan A. Shetty
#'
#' @references
#' Shetty SA (2020). Utilities for microbiome analytics.
#' \url{https://github.com/microsud/biomeUtils}
#'
#'
NULL

#' @rdname chekFilterEffect
#' @aliases checkAbundanceCutOffs
#' @importFrom dplyr as_tibble %>%
#' @importFrom phyloseq ntaxa
#' @export
checkAbundanceCutOffs <- function(x,
                                  abundance_thres = NULL,
                                  prevalence_thres = NULL){

  if (class(x) != "phyloseq") {
    stop("Input is not an object of phyloseq class")
  }
  if(is.null(abundance_thres) || is.null(prevalence_thres)){
    stop("'abundance_thres' and 'prevalence_thres' must be numeric values")
  }

  if(length(abundance_thres) == 1){
    warning("Only one 'abundance_thres' value provided for checking")
  }

  if(length(prevalence_thres) >1){
    warning("Multiple 'prevalence_thres' value provided, using only first value for checking")
    prevalence_thres <- prevalence_thres[1]
  }

  ab.cuts <- lapply(abundance_thres, function(th) {
    length(microbiome::core_members(x,
                                    detection = th,
                                    prevalence = prevalence_thres))
  })
  ab.cuts <- do.call("cbind", ab.cuts)
  colnames(ab.cuts) <- as.character(abundance_thres)
  #ntaxa(x) - ab.cuts
  xdf <- as_tibble(ntaxa(x) - ab.cuts) %>%
    t() %>%
    as.data.frame() %>%
    tibble::rownames_to_column("Abundance")
  colnames(xdf)[2] <- "Taxa_lost"
  xdf$TotalTaxa <- ntaxa(x)
  xdf$PrevalenceCutOff <- prevalence_thres

  return(xdf)

}

#
#' @rdname chekFilterEffect
#' @aliases checkPrevalenceCutOffs
#' @importFrom dplyr as_tibble %>%
#' @importFrom phyloseq ntaxa
#' @export
checkPrevalenceCutOffs <- function(x,
                                   abundance_thres = NULL,
                                   prevalence_thres = NULL){

  if(is.null(abundance_thres) || is.null(prevalence_thres)){
    stop("'abundance_thres' and 'prevalence_thres' must be numeric values")
  }

  if(length(prevalence_thres) == 1){
    warning("Only one 'prevalence_thres' value provided for checking")
  }

  if(length(abundance_thres) >1){
    warning("Multiple 'abundance_thres' value provided, using only first value for checking")
    abundance_thres <- abundance_thres[1]
  }

  prev.cuts <- lapply(prevalence_thres, function(th) {
    length(microbiome::core_members(x,
                                    detection = abundance_thres,
                                    prevalence = th))
  })
  prev.cuts <- do.call("cbind", prev.cuts)
  colnames(prev.cuts) <- as.character(prevalence_thres)

  xdf <- as_tibble(ntaxa(x) - prev.cuts) %>%
    t() %>%
    as.data.frame() %>%
    tibble::rownames_to_column("Prevalence")
  colnames(xdf)[2] <- "Taxa_lost"
  xdf$TotalTaxa <- ntaxa(x)
  xdf$AbundanceCutOff <- abundance_thres

  return(xdf)

}


