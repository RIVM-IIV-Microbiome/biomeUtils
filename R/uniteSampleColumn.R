#' unite a Column in Phyloseq Sample Data
#'
#' @name uniteSampleColumn
#'
#' @details In some cases user may want to unite a sample column in
#'          \code{phyloseq} object. To ease this step, we make use of the
#'          \code{\link[tidyr:unite]{unite}}. The separation can be done
#'          with regular expression or numeric locations.
#
#' @param x \code{phyloseq} object
#'
#' @param ... Option to pass on to \code{\link[tidyr:unite]{unite}}.
#'
#' @examples
#' library(biomeUtils)
#' library(dplyr)
#' data(FuentesIliGutData)
#' FuentesIliGutData <- uniteSampleColumn(FuentesIliGutData,
#'                                        "united.column", sex:age,
#'                                        remove = FALSE)
#'
#' @return Phyloseq with modified sample_data
#'
#' @author Sudarshan A. Shetty
#'
#' @references
#' Shetty SA (2020). Utilities for microbiome analytics.
#' \url{https://github.com/microsud/biomeUtils}
#'
#' @importFrom tidyr unite
#' @importFrom dplyr %>%
#' @importFrom phyloseq sample_data
#'
#' @export
#'
uniteSampleColumn <- function(x, ...){

  if (class(x) != "phyloseq") {
    stop("Input is not an object of phyloseq class")
  }

  if (is.null(x@sam_data)) {
    stop("uniteSampleColumn requires sample data")
  }

  samtb <- getSampleTibble(x, column_id = ".m_holder") %>%
    tidyr::unite(...) %>%
    getSampleTableFromTibble(rownames = ".m_holder")

  phyloseq::sample_data(x) <- samtb
  return(x)
}
