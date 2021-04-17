#' Unite Genus Species Names
#'
#' @name uniteGenusSpeciesNames
#'
#' @details Unite genus and species names in tax_table \code{phyloseq}
#'          objects. E.g. \emph{Fusicatenibacter} and \emph{saccharivorans} are united
#'          to \emph{Fusicatenibacter.saccharivorans}.
#'
#' @param x A phyloseq object
#'
#' @examples
#' library(biomeUtils)
#' data("FuentesIliGutData")
#' ps <- uniteGenusSpeciesNames(FuentesIliGutData)
#'
#' head(tax_table(ps)[, 7])
#' @return A \code{phyloseq} with united genus and species
#'         in taxa_table .
#'
#' @author Sudarshan A. Shetty
#'
#' @references
#' Shetty SA (2021). Utilities for microbiome analytics.
#' \url{https://github.com/microsud/biomeUtils}
#'
#' @importFrom dplyr select mutate rename
#' @importFrom phyloseq tax_table 'tax_table<-'
#' @importFrom tibble column_to_rownames
#'
#' @export

uniteGenusSpeciesNames <- function(x) {

  # global vars
  Genus.Species <- Genus <- Species <- NULL
  if (class(x) != "phyloseq") {
    stop("Input is not an object of phyloseq class")
  }

  tax_tb <- getTaxaTibble(x) %>%
    mutate(Genus.Species = ifelse(!is.na(Species), paste0(Genus, ".", Species), Species)) %>%
    select(-Species) %>%
    rename(Species = Genus.Species) %>%
    column_to_rownames("FeatureID")

  # tax_tb[1:30, 5:9]
  # rownames(tax_tb) <- tax_tb[,1]
  # tax_tb <- tax_tb[,-1]
  tax_table(x) <- tax_table(as.matrix(tax_tb))
  return(x)
  # cols <- c("Genus", "Species")
  # if("Genus" %in% rank_names(x) && "Species" %in% rank_names(x)) {

  #
  # } else {
  #
  # stop(" 'rank_names(x)' do not include columns with
  #      names 'Genus' and 'Species' ")

  # }
}
