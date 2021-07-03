#' Get Prevalence and Taxonomy
#'
#' @name getPrevalence
#'
#' @details Get the prevalence of taxa in \code{phyloseq} objects
#'          along with taxonomic classification.
#'
#' @param x A phyloseq object
#'
#' @param return_rank Specify which taxonomic ranks to include in output.
#'                    Must be a character vector \code{phyloseq::rank_names()}
#'
#' @param return_taxa A specific list of taxa for which the values should
#'                    be returned. This can be used if user is not interested
#'                    in all the taxa in input \code{phyloseq}. Default is NULL
#'                    which returns all taxa. This list must match rows names if
#'                    otu_table in phyloseq has taxa_are_rows=TRUE or columns
#'                    names if otu_table in phyloseq has taxa_are_rows=FALSE
#'
#' @param sort Logical. Sort by prevalence value from higher to lower (Default=TRUE)
#'
#' @param ... Option to pass microbiome::prevalence
#'
#' @examples
#' library(biomeUtils)
#' data("FuentesIliGutData")
#' prev_tib <- getPrevalence(FuentesIliGutData,
#'   return_rank = c("Family", "Genus"),
#'   return_taxa = c("ASV4", "ASV17", "ASV85", "ASV83"),
#'   sort = TRUE
#' )
#' head(prev_tib)
#' @return A tibble with prevalence and taxonomy
#'
#' @author Sudarshan A. Shetty
#'
#' @references
#' \itemize{
#' \item{}{A Salonen et al. (2012) The adult intestinal core microbiota
#'         is determined by  analysis depth and health status.
#'        \emph{Clinical Microbiology and Infection}, 18(S4):16 20.
#'        \url{https://doi.org/10.1111/j.1469-0691.2012.03855.x}
#' }
#' \item{}{Lahti L, Shetty S (2012-2019). microbiome R package.
#'         statistical aspects of the study of the microbiome.
#'         \emph{BioConductor}.
#'         \url{https://doi.org/doi:10.18129/B9.bioc.microbiome}
#' }
#' }
#'
#'
#' @importFrom microbiome prevalence
#' @importFrom tibble rownames_to_column
#' @importFrom phyloseq rank_names tax_table
#' @importFrom dplyr arrange left_join select rename filter desc
#' @importFrom rlang syms
#'
#' @export

getPrevalence <- function(x,
                          return_rank = rank_names(x),
                          return_taxa = taxa_names(x),
                          sort = TRUE, ...) {

  # global vars
  Taxa <- NULL

  if (class(x) != "phyloseq") {
    stop("Input is not an object of phyloseq class")
  }

  rank_nm <- .check_ranks(x, return_rank)
  rank_nm <- c("Taxa", rank_nm)

  tax_tib <- getTaxaTibble(x,
    column_id = "Taxa",
    select_cols = rank_names(x),
    select_rows = taxa_names(x)
  ) %>%
    dplyr::select(!!!syms(rank_nm))

  prev_tbl <- prevalence(x, ...) %>%
    as.data.frame() %>%
    rownames_to_column("Taxa") %>%
    rename(prevalence = ".") %>%
    left_join(tax_tib, by = "Taxa")

  # colnames(prev_tbl) <- c("Taxa", "prevalence", rank_nm)

  taxa_rt <- .check_taxa(x, return_taxa)

  prev_tbl <- prev_tbl %>%
    filter(Taxa %in% taxa_rt)

  if (sort) {
    prev_tbl <- prev_tbl %>%
      arrange(desc(prevalence))
    return(prev_tbl)
  } else {
    return(prev_tbl)
  }
}


#' @importFrom phyloseq rank_names
.check_ranks <- function(x, return_rank) {
  if (!is.null(return_rank) || any(return_rank %in% rank_names(x))) {
    return(return_rank)
  } else if (is.null(return_rank) || is.na(return_rank)) {
    return_rank <- rank_names(x)
    return(return_rank)
  } else if (!is.null(return_rank) && !any(return_rank %in% rank_names(x))) {
    stop("Please provide valid taxonomic rank names in 'return_rank' ")
  }
}

#' @importFrom phyloseq taxa_names
.check_taxa <- function(x, return_taxa) {
  if (!is.null(return_taxa) || any(return_taxa %in% taxa_names(x))) {
    return(return_taxa)
  } else if (is.null(return_taxa)) {
    return(taxa_names(x))
  } else if (!is.null(return_taxa) && !any(return_taxa %in% taxa_names(x))) {
    stop("Please provide valid names in 'return_taxa' ")
  }
}
