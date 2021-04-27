#' Filter Phyloseq Object
#'
#' @name filterPhyloseq
#'
#' @details These are alternative to subset and prune functions in
#'          \code{phyloseq}.
#'
#' The \code{filterSampleData} does subsetting similar to
#' \code{phyloseq::subset_samples}.
#' The \code{filterTaxaData} does subsetting similar to
#' \code{phyloseq::subset_taxa}.
#'
#'
#' The \code{filterSampleByNames} does pruning similar to
#' \code{phyloseq::prune_samples}.
#' The \code{filterTaxaByNames} does pruning similar to
#' \code{phyloseq::prune_taxa}.
#'
#' The naming is using filter because it uses the \code{dplyr::filter} function.
#'
#' Note:
#' \code{filterSampleData} will additionally remove any taxa that are
#' zero across all samples i.e. not present is any remaining samples
#' after filtering.
#' \code{filterTaxaData} will additionally remove samples that do not
#' have any of the selected taxa, i.e. samples that sum to zero.
#' \code{filterSampleByNames} will additionally remove any taxa that are not
#' present is any remaining samples
#' after filtering.
#' \code{filterTaxaByNames} will additionally remove samples that do not
#' have any of the selected taxa, i.e. samples that sum to zero.
#'
#'
#'
#'
#' @param x \code{phyloseq} object
#'
#' @param ... Option to pass on to \code{\link[dplyr]{filter}}. Only when using
#'                      \code{filterSampleData} or \code{filterTaxaData}
#'                      sample data or taxonomy table as tibble.
#'
#' @param ids A list of sample ids or taxa ids to filter. Only when using
#'                      \code{filterSampleByNames} or \code{filterTaxaByNames}
#'
#' @param keep Logical. Default is TRUE. Only when using
#'                      \code{filterSampleByNames} or \code{filterTaxaByNames}
#'
#' @examples
#'
#' library(biomeUtils)
#' data("FuentesIliGutData")
#'
#' # Filter from tables subset_*
#' ps.filtered.samples  <- filterSampleData(FuentesIliGutData,
#'                                          ILI == "C" & BMI < 26)
#'
#' ps.filtered.samples
#'
#' ps.filtered.taxa  <- filterTaxaData(FuentesIliGutData,
#'                                     Phylum=="Firmicutes")
#'
#' ps.filtered.taxa
#'
#' # Filter by names like prune_*
#' sams.select <- sample_names(FuentesIliGutData)[1:10]
#' ps.filter.by.sam.names <- filterSampleByNames(FuentesIliGutData,
#'                                               ids = sams.select,
#'                                               keep = TRUE)
#'
#' tax.select <- taxa_names(FuentesIliGutData)[1:10]
#' ps.filter.by.tax.names <- filterTaxaByNames(FuentesIliGutData,
#'                                             ids = tax.select,
#'                                             keep = TRUE)
#'
#' @return Either filtered \code{phyloseq}
#'
#' @author Sudarshan A. Shetty
#'
#' @references
#' Shetty SA (2020). Utilities for microbiome analytics.
#' \url{https://github.com/microsud/biomeUtils}
#'
#'
NULL

#' @rdname filterPhyloseq
#' @aliases filterSampleData
#' @importFrom dplyr filter
#' @importFrom phyloseq sample_data
#' @export
#'
filterSampleData <- function(x, ...) {

  if (class(x) != "phyloseq") {
    stop("Input is not an object of phyloseq class")
  }

  #rownames <- .check_col_sam_data(x,rownames)


  y <- x %>%
    getSampleTibble() %>%
    filter(...)

  # getSampleTibble by default adds the rownames as the first column
  # this will now be again added as rownames to sample_data
  y <-  y %>%
    tibble::column_to_rownames(colnames(y)[1]) %>%
    phyloseq::sample_data()

  sample_data(x) <- y

  return(removeZeros(x))


}


#' @rdname filterPhyloseq
#' @aliases filterTaxaData
#' @importFrom dplyr filter
#' @importFrom phyloseq tax_table
#' @export
#'
filterTaxaData <- function(x, ...) {

  if (class(x) != "phyloseq") {
    stop("Input is not an object of phyloseq class")
  }
  # check for sample names

  y <- x %>%
    getTaxaTibble() %>%
    filter(...)

  y <-  y %>%
    tibble::column_to_rownames(colnames(y)[1]) %>% as.matrix() %>%
    phyloseq::tax_table()

  tax_table(x) <- y

  return(removeZeros(x))


}


#' @rdname filterPhyloseq
#' @aliases filterTaxaByNames
#' @importFrom dplyr filter
#' @importFrom phyloseq tax_table
#' @export
#'
filterTaxaByNames <- function(x,
                              ids = NULL,
                              keep = TRUE) {

  # global var
  sel.tx <- NULL

  if (class(x) != "phyloseq") {
    stop("Input is not an object of phyloseq class")
  }

  if (is.null(ids)) {
    warning("No taxa to discard. Returning original pseq object")
    return(x)
  }


  if (keep){

    y <- x %>%
      getTaxaTibble(column_id = "sel.tx") %>%
      dplyr::filter(sel.tx %in% ids) %>%
      getTaxaTableFromTibble(rownames="sel.tx")

    phyloseq::tax_table(x) <- y

    return(removeZeros(x))

  } else {

    y <- x %>%
      getTaxaTibble(column_id = "sel.tx") %>%
      dplyr::filter(!(sel.tx %in% ids)) %>%
      getTaxaTableFromTibble(rownames="sel.tx")

    phyloseq::tax_table(x) <- y

    return(removeZeros(x))

  }


}



#' @rdname filterPhyloseq
#' @aliases filterSampleByNames
#' @importFrom dplyr filter
#' @importFrom phyloseq tax_table
#' @export
#'
filterSampleByNames <- function(x,
                                ids = NULL,
                                keep = TRUE) {

  # global var
  sel.sm <- NULL
  if (is.null(ids)) {
    warning("No samples to discard. Returning original pseq object")
    return(x)
  }


  if (class(x) != "phyloseq") {
    stop("Input is not an object of phyloseq class")
  }

  if (keep){

    y <- x %>%
      getSampleTibble(column_id = "sel.sm") %>%
      dplyr::filter(sel.sm %in% ids) %>%
      getSampleTableFromTibble(rownames="sel.sm")

    phyloseq::sample_data(x) <- y

    return(removeZeros(x))

  } else {

    y <- x %>%
      getSampleTibble(column_id = "sel.sm") %>%
      dplyr::filter(!(sel.sm %in% ids)) %>%
      getSampleTableFromTibble(rownames="sel.sm")

    phyloseq::sample_data(x) <- y

    return(removeZeros(x))

  }


}




