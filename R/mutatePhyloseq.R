#' Mutating Slots of Phyloseq
#'
#' @name mutatePhyloseq
#'
#' @details Mutate different \code{phyloseq} slots.
#'
#' \code{mutateAbundanceTable} mutate \code{otu_table}. Will also add the new sample to
#'  \code{sample_data}
#'
#' \code{mutateTaxaTable} mutate \code{tax_table}, Add or modify existing columns.
#'
#' \code{mutateSampleData} mutate \code{sample_data}. Add columns, modify existing
#' columns.
#'
#' @param x \code{phyloseq} object
#'
#' @param ... Arguments that dplyr::mutate can accept
#'
#' @examples
#' library(biomeUtils)
#' library(dplyr)
#' data("FuentesIliGutData")
#' ps_otu <- mutateAbundanceTable(FuentesIliGutData, test = 1:ntaxa(FuentesIliGutData))
#' ps_tax <- mutateTaxaTable(FuentesIliGutData, test = 1:ntaxa(FuentesIliGutData))
#' ps_samdat <- mutateSampleData(FuentesIliGutData, test = paste0(ILI, "-test"))
#' @return A phyloseq
#'
#' @author Sudarshan A. Shetty
#'
#' @references
#' Shetty SA (2020). Utilities for microbiome analytics.
#' \url{https://github.com/RIVM-IIV-Microbiome/biomeUtils}
#'
#'
NULL


#' @rdname mutatePhyloseq
#' @aliases mutateSampleData
#' @importFrom dplyr %>% mutate
#'
#' @export
#'

mutateSampleData <- function(x, ...) {

  # A place holder for storing/preserving sample_names
  # during data organiszation
  .m_holder <- NULL

  if (!is(x, "phyloseq")) {
    stop("Input is not an object of phyloseq class")
  }

  if (is.null(x@sam_data)) {
    stop("mutateSample requires sample data")
  }


  samtb <- getSampleTibble(x, column_id = ".m_holder") %>%
    dplyr::mutate(...) %>%
    getSampleTableFromTibble(rownames = ".m_holder")

  phyloseq::sample_data(x) <- samtb

  return(x)
}


#' @rdname mutatePhyloseq
#' @aliases mutateTaxaTable
#' @importFrom dplyr %>% mutate
#'
#' @export
#'
mutateTaxaTable <- function(x, ...){


  # A place holder for storing/preserving sample_names
  # during data organiszation
  .m_holder <- txtb <- NULL

  if (!is(x, "phyloseq")) {
    stop("Input is not an object of phyloseq class")
  }

  if (is.null(x@tax_table)) {
    stop("mutateTaxa requires tax_table data")
  }


  txtb <- getTaxaTibble(x, column_id=".m_holder") %>%
    dplyr::mutate(...) %>%
    getTaxaTableFromTibble(rownames = ".m_holder")

  phyloseq::tax_table(x) <- txtb

  return(x)
}


#' @rdname mutatePhyloseq
#' @aliases mutateAbundanceTable
#' @importFrom dplyr %>% mutate bind_rows
#' @importFrom phyloseq otu_table sample_names nsamples sample_variables sample_data
#'
#' @export
#'
mutateAbundanceTable <- function(x, ...){


  # A place holder for storing/preserving sample_names
  # during data organiszation
  .m_holder <- NULL
  otutb <- new_sams <- dum.sam <- NULL

  if (!is(x, "phyloseq")) {
    stop("Input is not an object of phyloseq class")
  }

  if (is.null(x@otu_table)) {
    stop("mutateAbundance requires otu_table data")
  }


  otutb <-  getAbundanceTibble(x, column_id=".m_holder") %>%
    dplyr::mutate(...) %>%
    getOtuTableFromTibble(rownames = ".m_holder")

  if (ncol(otutb) > nsamples(x)){

    new_sams <- setdiff(colnames(otutb), sample_names(x))

    # here a dummy sample data needs to be added for new sample
    #make dummy sample_data
    dum.sam <- data.frame(row.names = new_sams,
                          .m.hodler = new_sams)

    dum.sam <- data.frame(matrix(nrow = length(new_sams),
                                 ncol = length(phyloseq::sample_variables(x))))
    colnames(dum.sam) <- phyloseq::sample_variables(x)
    rownames(dum.sam) <- new_sams
    dum.sam <- dplyr::bind_rows(meta(x), dum.sam) %>%
      phyloseq::sample_data()

    # remove old and add new
    x@otu_table <- NULL
    x@sam_data <- NULL
    phyloseq::otu_table(x) <- otutb
    phyloseq::sample_data(x) <- dum.sam
    return(x)
  } else {

    phyloseq::otu_table(x) <- otutb

    return(x)

  }

}




