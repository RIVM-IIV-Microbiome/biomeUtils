#' Summarize Taxonomic Assignments
#'
#' @name summarizeTaxonomicAssignments
#'
#' @details Percentage of taxa that have been assigned to the different
#'          taxonomic levels. Non assigned levels must be NA. Check
#'          biomeUtils::cleanTaxonomy.
#'
#' @param x A phyloseq object
#'
#' @examples
#' library(biomeUtils)
#' data("SprockettTHData")
#' # run cleanTaxonomy if special characters are
#' # present in tax_table
#' summarizeTaxonomicAssignments(SprockettTHData)
#' @return tibble with overview of taxonomic assignments.
#'
#' @author Sudarshan A. Shetty
#'
#' @importFrom phyloseq tax_table taxa_sums
#' @importFrom tibble tibble
#'
#' @export

summarizeTaxonomicAssignments <- function(x) {
  if (class(x) != "phyloseq") {
    stop("Input is not an object of phyloseq class")
  }
  taxadf <- tax_table(x) %>%
    as.matrix() %>%
    as.data.frame(stringsAsFactors = FALSE)


  whichNA <- apply(taxadf, 2, function(x) {
    sum(is.na(x))
  })


  whichNonNA <- apply(taxadf, 2, function(x) {
    sum(!is.na(x))
  })

  ## dada2 especially species assignments can be
  if (any(colnames(taxadf) == "Species")) {
    whichAmbiguous <- sum(grepl("/", taxadf$Species))
    message(paste0("In Species column- ", whichAmbiguous, " are ambigous assignments"))
  }


  total <- whichNA + whichNonNA

  TaxonomicAssignments <- tibble::tibble(
    Ranks = rank_names(x),
    TaxonomyAssigned = whichNonNA,
    Total = total,
    PrecentAssigned = round(100 * whichNonNA / total, 1)
  )
  return(TaxonomicAssignments)
}
