#' Get Phyloseq slots from Tibbles
#'
#' @name getPhyloseqSlots
#'
#' @details Convert different tibbles into \code{phyloseq} slots.
#'
#' \code{getOtuTableFromTibble} Converts OTU tibble obtained using
#' \code{\link[biomeUtils]{getAbundanceTibble}} to \code{\link[phyloseq]{otu_table}}.
#'
#' \code{getTaxaTableFromTibble} Converts the taxa tibble obtained using
#' \code{\link[biomeUtils]{getTaxaTibble}} to \code{\link[phyloseq]{tax_table}}.
#'
#' \code{getSampleTableFromTibble} Converts the sample data tibble obtained using
#' \code{\link[biomeUtils]{getSampleTibble}} to \code{\link[phyloseq]{sample_data}}.
#'
#' @param x \code{phyloseq} object
#'
#' @param rownames A column in tibble to use are rownames
#'
#' @examples
#' library(biomeUtils)
#' data("FuentesIliGutData")
#' otu_tib  <- getAbundanceTibble(FuentesIliGutData)
#' otu_table <- getOtuTableFromTibble(otu_tib, rownames="FeatureID")
#' class(otu_table)
#'
#' tax_tib <- getTaxaTibble(FuentesIliGutData)
#' taxa_table <- getTaxaTableFromTibble(tax_tib, rownames="FeatureID")
#' class(tax_table)
#'
#' meta_tib <- getSampleTibble(FuentesIliGutData)
#' sample_table <- getSampleTableFromTibble(meta_tib, rownames="SampleID")
#' class(sample_table)
#'
#' @return Either \code{\link[phyloseq]{otu_table}} or
#'                \code{\link[phyloseq]{tax_table}} or
#'                \code{\link[phyloseq]{sample_data}}
#'
#' @author Sudarshan A. Shetty
#'
#' @references
#' Shetty SA (2020). Utilities for microbiome analytics.
#' \url{https://github.com/microsud/biomeUtils}
#'
#'
NULL


#' @rdname getPhyloseqSlots
#' @aliases getOtuTableFromTibble
#' @importFrom tibble column_to_rownames is_tibble
#' @importFrom phyloseq otu_table
#' @export
#'
getOtuTableFromTibble <- function(x, rownames="FeatureID") {

  if (!is_tibble(x)) {
    stop("Input is not a tibble")
  }

  otu_table <- x %>%
    tibble::column_to_rownames(rownames) %>%
    phyloseq::otu_table(taxa_are_rows = T)

  return(otu_table)

}

#' @rdname getPhyloseqSlots
#' @aliases getTaxaTableFromTibble
#' @importFrom tibble column_to_rownames is_tibble
#' @importFrom phyloseq tax_table
#' @export
#'
getTaxaTableFromTibble <- function(x, rownames="FeatureID") {

  if (!is_tibble(x)) {
    stop("Input is not a tibble")
  }

  taxa_table <- x %>%
    tibble::column_to_rownames(rownames) %>%
    as.matrix() %>%
    phyloseq::tax_table()

  return(taxa_table)

}


#' @rdname getPhyloseqSlots
#' @aliases getSampleTableFromTibble
#' @importFrom tibble column_to_rownames is_tibble
#' @importFrom phyloseq sample_data
#' @importFrom dplyr mutate
#' @export
#'
getSampleTableFromTibble <- function(x, rownames= NULL) {

  if (!is_tibble(x)) {
    stop("Input is not a tibble")
  }


  if (is.null(rownames) || is.na(rownames)) {

    x$rowname_2 <- paste0("Sample-", seq(1:nrow(x)))

    sam_table <- x %>%
      tibble::column_to_rownames("rowname_2") %>%
      phyloseq::sample_data()

    return(sam_table)


  } else{
    #rownames <- .check_col_sam_data(vb,rownames)

    sam_table <- x %>%
      tibble::column_to_rownames(rownames) %>%
      phyloseq::sample_data()

    return(sam_table)
  }

}

