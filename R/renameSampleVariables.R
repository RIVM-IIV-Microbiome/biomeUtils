#' Rename Sample Variables in Phyloseq Sample Data
#'
#' @name renameSampleVariables
#'
#' @details In some cases user may want to rename one or several sample variables.
#'          With \code{phyloseq}, this needs to be done one by one. To ease this step,
#'          we make use of the \code{\link[dplyr:rename]{rename}} and
#'          \code{\link[dplyr:rename_with]{rename_with}}. The later allows to use
#'          functions such as those that work with string matching.
#
#' @param x \code{phyloseq} object
#'
#' @param ... Option to pass on to \code{\link[dplyr:rename]{rename}} or
#'            \code{\link[dplyr:rename_with]{rename_with}}. If using later,
#'            set with=TRUE.
#'
#' @param with Default is FALSE. Set this to TRUE when using a function. See
#'             \code{\link[dplyr:rename_with]{rename_with}}
#'
#' @examples
#' library(biomeUtils)
#' library(dplyr)
#' data(FuentesIliGutData)
#' FuentesIliGutData <- renameSampleVariables(FuentesIliGutData,
#'                                            with=FALSE,
#'                                            condition = "ILI")
#'
#' FuentesIliGutData <- renameSampleVariables(FuentesIliGutData,
#'                                            with=TRUE,
#'                                            .fn = toupper)
#'
#' @return Phyloseq with modified sample_data
#'
#' @author Sudarshan A. Shetty
#'
#' @references
#' Shetty SA (2020). Utilities for microbiome analytics.
#' \url{https://github.com/RIVM-IIV-Microbiome/biomeUtils}
#'
#' @importFrom dplyr rename_with rename %>%
#' @importFrom phyloseq sample_data
#'
#' @export
#'
renameSampleVariables <- function(x, with=FALSE, ...) {

  s.tb <- holder <- NULL
  if(with){
    s.tb <- getSampleTibble(x, column_id = "holder") %>%
      dplyr::rename_with(...)
    s.tb <- s.tb %>%
      getSampleTableFromTibble(rownames=colnames(s.tb)[1])
    sample_data(x) <- s.tb
    return(x)
  }

  s.tb <- getSampleTibble(x, column_id = "holder") %>%
    dplyr::rename(...)
  s.tb <- s.tb %>%
    getSampleTableFromTibble(rownames="holder")
  sample_data(x) <- s.tb
  return(x)

}
