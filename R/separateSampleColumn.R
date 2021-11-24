#' Separate a Column in Phyloseq Sample Data
#'
#' @name separateSampleColumn
#'
#' @details In some cases user may want to separate a sample column in
#'          \code{phyloseq} object. To ease this step, we make use of the
#'          \code{\link[tidyr:separate]{separate}}. The separation can be done
#'          with regular expression or numeric locations.
#
#' @param x \code{phyloseq} object
#'
#' @param ... Option to pass on to \code{\link[tidyr:separate]{separate}}.
#'
#' @examples
#' library(biomeUtils)
#' library(dplyr)
#' data(FuentesIliGutData)
#' FuentesIliGutData <- mutateSampleData(FuentesIliGutData,
#'                                       column.to.split = paste0(sex,"_",age))
#' FuentesIliGutData <- separateSampleColumn(FuentesIliGutData,
#'                                           column.to.split, c("A", "B"), sep="_")
#'
#' @return Phyloseq with modified sample_data
#'
#' @author Sudarshan A. Shetty
#'
#' @references
#' Shetty SA (2020). Utilities for microbiome analytics.
#' \url{https://github.com/RIVM-IIV-Microbiome/biomeUtils}
#'
#' @importFrom tidyr separate
#' @importFrom dplyr %>%
#' @importFrom phyloseq sample_data
#'
#' @export
#'
separateSampleColumn <- function(x, ...){

  if (class(x) != "phyloseq") {
    stop("Input is not an object of phyloseq class")
  }

  if (is.null(x@sam_data)) {
    stop("separateSampleColumn requires sample data")
  }

  samtb <- getSampleTibble(x, column_id = ".m_holder") %>%
    tidyr::separate(...) %>%
    getSampleTableFromTibble(rownames = ".m_holder")

  phyloseq::sample_data(x) <- samtb
  return(x)

}
