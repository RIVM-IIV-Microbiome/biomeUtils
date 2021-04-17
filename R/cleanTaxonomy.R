#' Clean Taxonomy
#'
#' @name cleanTaxonomy
#'
#' @details Removes, special characters, spaces and prefixes in tax_table.
#'
#' @param x A phyloseq object
#'
#' @examples
#' library(biomeUtils)
#' data("FuentesIliGutData")
#' ps1 <- cleanTaxonomy(FuentesIliGutData)
#' ps1
#' @return phyloseq
#'
#' @author Sudarshan A. Shetty
#'
#' @importFrom phyloseq tax_table
#' @importFrom dplyr mutate_if na_if
#' @importFrom tibble column_to_rownames
#'
#' @export

cleanTaxonomy <- function(x) {
  tax_table(x)[, colnames(tax_table(x))] <- gsub(tax_table(x)[, colnames(tax_table(x))],
    pattern = "[a-z]__", replacement = ""
  )

  tax_table(x)[, colnames(tax_table(x))] <- gsub(tax_table(x)[, colnames(tax_table(x))],
    pattern = "\\[|\\]", replacement = ""
  )

  tax_table(x)[, colnames(tax_table(x))] <- gsub(tax_table(x)[, colnames(tax_table(x))],
    pattern = " ", replacement = "."
  )

  tax_table(x)[, colnames(tax_table(x))] <- gsub(tax_table(x)[, colnames(tax_table(x))],
    pattern = "-", replacement = "."
  )
  tib <- getTaxaTibble(x) %>%
    dplyr::mutate_if(is.character, list(~ na_if(., ""))) %>%
    tibble::column_to_rownames("FeatureID") %>%
    as.matrix()

  tax_table(x) <- tib

  return(x)
}
