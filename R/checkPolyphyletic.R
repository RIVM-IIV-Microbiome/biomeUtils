#' Check for Polyphyletic Taxa
#'
#' @name checkPolyphyletic
#'
#' @details Check for polyphyletic taxa in \code{tax_table}. Useful to check this
#'          before aggregating at any level.
#'
#' @param x \code{phyloseq} object
#'
#' @param taxa_level Taxonomic level to check. Default is Genus
#'
#' @param return_df Logical. Return data.frame or list
#'
#' @examples
#' library(biomeUtils)
#' library(dplyr)
#' data("FuentesIliGutData")
#' polydf <- checkPolyphyletic(FuentesIliGutData)
#' polydf
#'
#' @return Data frame or list
#'
#' @author Sudarshan A. Shetty
#'
#' @references
#' Shetty SA (2020). Utilities for microbiome analytics.
#' \url{https://github.com/RIVM-IIV-Microbiome/biomeUtils}
#'
#' @importFrom phyloseq tax_table
#' @importFrom dplyr select distinct group_by_at count %>% vars
#' @importFrom rlang sym
#' @export
NULL
checkPolyphyletic <- function(x,
                              taxa_level = "Genus",
                              return_df = TRUE){

  nfeatures <- NULL
  if (!is(x, "phyloseq")) {
    stop("Input is not an object of phyloseq class")
  }

  if(is.null(x@tax_table)){
    stop("tax_table is empty")
    }

  tax <- getTaxaTibble(x)
  num.level <- which(colnames(tax)==taxa_level)
  tax.tb <- .count_tax_tib(tax, level=taxa_level)

  if(!.check_names_dup(tax.tb)){
    stop(paste("No Polyphyly detected"))
  }

  if(!return_df){
    tax.tb <- tax.tb %>%
      dplyr::filter(n != 1) %>%
      dplyr::pull(!!sym(taxa_level))

  } else {

    cols.sel <- 2:num.level
    tax.tb <- tax %>%
      dplyr::select(dplyr::all_of(cols.sel)) %>%
      dplyr::distinct() %>%
      dplyr::group_by(!!sym(taxa_level)) %>%
      dplyr::mutate(nfeatures = dplyr::n()) %>%
      dplyr::filter(!is.na(!!sym(taxa_level)), nfeatures > 1)

  }

  return(tax.tb)

}





#' @importFrom dplyr select distinct group_by_at count %>% vars
#' @importFrom rlang sym
#' @keywords internal
# param tax is output of getTaxaTibble
# param level Taxonomic level to investigate
.count_tax_tib <- function(tax, level="Genus") {

  num.level <- which(colnames(tax)==level)

  cnt.tx <- tax %>%
    dplyr::select(c(2:num.level)) %>% #getTibble returns ASVs as first col
    dplyr::distinct() %>%
    dplyr::group_by_at(dplyr::vars(level)) %>%
    dplyr::count() %>%
    dplyr::filter(!is.na(!!sym(level)))

}


#' @importFrom dplyr %>%
#' @importFrom rlang .data
#' @keywords internal
# param tax.tb is output of .count_tax_tib
.check_names_dup <- function(tax.tb) {

  n <- NULL
  tax.tb %>% {all(n == 1)}

}

