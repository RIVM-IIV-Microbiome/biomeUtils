#' Get Taxa Summary at Specific Taxonomic Level
#'
#' @name getTaxaSummary
#'
#' @details Get a summary of total abundances of taxa
#'          at specific taxonomic level.
#'
#' @param x A phyloseq object
#'
#' @param rank Taxonomic rank for which abundance summaries
#'             are required.
#'
#' @examples
#' library(biomeUtils)
#' library(dplyr)
#' data("SprockettTHData")
#' x <- subset_samples(SprockettTHData,
#'                     Delivery_Mode=="Vaginal" & Cohort == "Tsimane")
#' x <- removeZeros(x)
#' getTaxaSummary(x, rank="Phylum")
#'
#' @importFrom data.table as.data.table melt
#' @importFrom dplyr left_join group_by ungroup mutate arrange desc summarise
#' @importFrom rlang .data sym
#'
#' @export

getTaxaSummary <-function(x, rank= "Phylum"){

  # global vars
  Abundance <- Counts <- Percent <- variable <- NULL

  suppressMessages(

    #improve speed for large number of rows and samples
    # with data.table
    abund_tb <- getAbundanceTibble(x) %>%
      data.table::as.data.table() %>%
      data.table::melt(id.vars = c("FeatureID"),
                       value.name = "Abundance") %>%
      dplyr::left_join(getTaxaTibble(x))

  )

  rnk <- .check_rank_single(x, rank = rank)

  overview <- abund_tb %>%
    group_by(!!sym(rnk)) %>%
    summarise(Counts=sum(Abundance)) %>%
    ungroup %>%
    mutate(Percent = Counts/sum(Counts)*100) %>%
    arrange(desc(Percent))

  return(overview)

}


#' @importFrom phyloseq rank_names
.check_rank_single <- function(x, rank) {

  if (is.null(rank) || is.na(rank) || length(rank) > 1 || is.numeric(rank)) {

    stop("Please provide valid taxonomic rank name in 'rank'.
         It cannot be multiple ranks or NULL or NA")

  } else if (any(rank %in% rank_names(x)) && length(rank)==1) {
    trank <- rank

    return(trank)
  }

}

