#' Microbiota Uniqueness
#'
#' @name uniqueness
#'
#' @details Extracts the minimum value from a \code{dissimilarity} matrix for each
#'          individual. This is the dissimilarity of an individual from their
#'          nearest neighbor. Here, the option of using a one or more reference
#'          samples is provided.
#'
#' @param x A phyloseq object
#'
#' @param dist_mat An object of class \code{dist}
#'
#' @param reference_samples Vector of samples to use as reference.
#'
#' @param ... Option to pass biomeUtils::getSampleTibble like column name and columns
#'            to be included in the output.
#'
#' @return A data frame with uniqueness and sample information
#'
#' @examples
#' library(biomeUtils)
#' data("FuentesIliGutData")
#' # Keep only two groups. controls and L1.
#' ps <- filterSampleData(FuentesIliGutData, ILI != "L2")
#' ps <- getProportions(ps)
#' dist.mat <- phyloseq::distance(ps, "bray")
#' # Define controls as reference samples
#' ref_samples <- rownames(meta(subset_samples(ps, ILI == "C")))
#' muniq <- uniqueness(ps,
#'                     dist_mat=dist.mat,
#'                     reference_samples = ref_samples)
#' head(muniq)
#' @references
#' \itemize{
#' \item{}{Wilmanski T et al. (2021). Gut microbiome pattern reflects
#' healthy ageing and predicts survival in humans.
#' \emph{Nature metabolism}, 3(2), pp.274-286.}
#' }
#'
#' @author Sudarshan A. Shetty
#'
#' @importFrom dplyr left_join rename filter summarise group_by
#'
#' @export

uniqueness <- function(x, dist_mat, reference_samples = NULL, ...){


  S1 <- S2 <- n <- value <- NULL
  dist.melt.sample <- meltDistanceToTable(x,
                                          dist_mat = dist_mat,...)

  # Check if reference samples provided
  if(!is.null(reference_samples)){

    unq.tib <- dist.melt.sample |>
      dplyr::filter(S1 %in% reference_samples) |>
      dplyr::group_by(S2) |>
      dplyr::summarise(uniquness=min(value))
  }
  # If reference samples calculate global uniqueness
  unq.tib <- dist.melt.sample |>
    dplyr::group_by(S2) |>
    dplyr::summarise(uniqueness=min(value))

  sam_tib <- getSampleTibble(x, ...)

  unq.tib <- unq.tib |>
    dplyr::left_join(sam_tib, by = c("S2" = colnames(sam_tib)[1]))

  return(unq.tib)

}
