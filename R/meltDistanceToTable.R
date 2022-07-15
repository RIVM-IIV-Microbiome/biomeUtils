#' Convert Distance matrix to Data Frame
#'
#' @name meltDistanceToTable
#'
#' @details Converts \code{dist} object from \code{phyloseq::distance()}
#'          to a data frame combined with sample data
#'
#' @param x A phyloseq object
#'
#' @param dist_mat An object of class \code{dist}
#'
#' @param name_dist_column A specific name for column where distance values
#'                         are stored. Default="value"
#'
#' @param ... Option to pass biomeUtils::getSampleTibble like column name and columns
#'            to be included in the output.
#'
#' @examples
#' library(biomeUtils)
#' data("FuentesIliGutData")
#' ps <- filterSampleData(FuentesIliGutData, ILI != "L2")
#' ps <- phyloseq::rarefy_even_depth(ps)
#' dist.mat <- phyloseq::distance(ps, "bray")
#' dist.melt.sample <- meltDistanceToTable(ps,
#'                                         dist_mat = dist.mat,
#'                                         name_dist_column = "Bray-Curtis",
#'                                         select_cols = c("participant_id", "ILI"))
#' head(dist.melt.sample)
#'
#' @return A data frame with pairwise distance and sample information
#'
#' @author Sudarshan A. Shetty
#'
#' @importFrom dplyr left_join rename filter
#'
#' @export

meltDistanceToTable <- function(x,
                                dist_mat=NULL,
                                name_dist_column = "value", ...){

  if(!is.character(name_dist_column)){
    stop("name_dist_column must be a single character value")
  }

  if(is(dist_mat)[1] != "dist"){
    stop("dist_mat required an object of class 'dist'")
  }

  # global var
  S1 <- S2 <- NULL

  sam_tib <- getSampleTibble(x, ...)

  dist.tidy <- dist_mat %>%
    as.matrix() %>%
    as.data.frame.table(responseName = name_dist_column) %>%
    dplyr::rename( S1 = "Var1", S2 = "Var2") %>%
    dplyr::left_join(sam_tib, by = c("S1" = colnames(sam_tib)[1])) %>%
    dplyr::left_join(sam_tib, by = c("S2" = colnames(sam_tib)[1]), suffix = c("_S1", "_S2")) %>%
    dplyr::filter(S1!=S2)

  return(dist.tidy)

}
