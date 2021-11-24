#' Calculate Fold Difference in Taxa Abundance
#'
#' @name calculateTaxaFoldDifference
#'
#' @details Calculate Log10 fold difference in abundance of taxa between two-
#'          groups. This code is modified from original code
#'          \url{https://github.com/microbiome/microbiome/blob/fa2c0de2fbe000da87be3c185972ed7f0f626591/inst/extdata/check_foldchange.R} Get the prevalence of taxa in \code{phyloseq} objects
#'          along with taxonomic classification and prevalence.
#'
#' @param x A phyloseq object
#'
#' @param group Vector with specifying the groups to compare.
#'              Only two-group comparisons are supported.
#'
#' @param sort Sort the results by descending order of fold difference
#'
#' @param paired Paired comparison (Default: FALSE)
#'
#' @examples
#' library(biomeUtils)
#' data("FuentesIliGutData")
#' # reduce size for example
#' ps1 <- filterSampleData(FuentesIliGutData, ILI != "L2")
#'
#' taxa_fd <- calculateTaxaFoldDifference(ps1, group = "ILI")
#' # check
#' taxa_fd
#' @return A tibble with taxa ids, taxonomic information,
#'         two-group prevalence and fold change values.
#'
#' @author Original Author: Leo Lahti. Adapted by Sudarshan A. Shetty
#'
#' @references
#' Shetty SA (2021). Utilities for microbiome analytics.
#' \url{https://github.com/RIVM-IIV-Microbiome/biomeUtils}
#'
#' @importFrom tibble rownames_to_column
#' @importFrom phyloseq sample_names prune_samples
#' @importFrom microbiome abundances meta
#' @importFrom dplyr left_join relocate rename arrange desc filter mutate
#'
#' @export

calculateTaxaFoldDifference <- function(x, group, sort = FALSE, paired = FALSE) {

  # global vars
  FoldDifference <- NULL
  # Pick the grouping variable from sample metadata
  g <- group
  if (length(g) == 1) {
    g <- sample_data(x)[[g]]
    if (!is.factor(g)) {
      warning(paste("Converting the grouping variable", group, "into a factor."))
      g <- as.factor(g)
    }
    g <- droplevels(g)
    if (!length(levels(g)) == 2) {
      stop(paste(
        "calculateTaxaFoldChange is currently implemented only for
                 two-group comparisons. The selected variable", group, "has",
        length(unique(g)), "levels: ", paste(unique(g), collapse = "/")
      ))
    }
  }
  if (is(x) == "phyloseq") {
    tx.ab <- abundances(x)
  }
  # Calculate fold changes
  if (paired) {
    fc <- apply(tx.ab, 1, function(xi) {
      spl <- split(xi, g)
      log10(mean(spl[[2]] - spl[[1]], na.rm = TRUE))
    })
  } else {
    fc <- apply(tx.ab, 1, function(xi) {
      spl <- split(xi, g)
      log10(mean(spl[[2]], na.rm = TRUE)) - log10(mean(spl[[1]], na.rm = TRUE))
    })
  }

  fcdf <- as.data.frame(fc) %>%
    rename(FoldDifference = fc) %>%
    rownames_to_column("FeatureID")


  lev.g1 <- levels(g)[1]
  lev.g2 <- levels(g)[2]

  g1.sams <- getSampleTibble(x,
    select_rows = sample_names(x),
    select_cols = group
  ) %>%
    filter(.data[[group]] %in% lev.g1)

  g2.sams <- getSampleTibble(x,
    select_rows = sample_names(x),
    select_cols = group
  ) %>%
    filter(.data[[group]] %in% lev.g2)

  lab.prev1 <- paste("Prevalence.", lev.g1, sep = "")
  lab.prev2 <- paste("Prevalence.", lev.g2, sep = "")

  prev.tb.g1 <- getPrevalence(prune_samples(sample_names(x) %in% g1.sams$SampleID, x),
    return_rank = rank_names(x),
    return_taxa = taxa_names(x)
  )
  colnames(prev.tb.g1)[2] <- lab.prev1

  prev.tb.g2 <- getPrevalence(prune_samples(sample_names(x) %in% g2.sams$SampleID, x),
    return_rank = rank_names(x),
    return_taxa = taxa_names(x)
  )
  colnames(prev.tb.g2)[2] <- lab.prev2


  prev.tb <- prev.tb.g1 %>%
    left_join(prev.tb.g2)

  fcdf <- fcdf %>%
    left_join(prev.tb, by = c("FeatureID" = "Taxa")) %>%
    mutate(Enriched = ifelse(FoldDifference > 0, lev.g2,
      ifelse(FoldDifference < 0, lev.g1, "NoChange")
    )) %>%
    dplyr::relocate(rank_names(x), .before = "FoldDifference")

  if (sort) {
    fcdf <- fcdf %>%
      arrange(desc(FoldDifference))
  }

  return(fcdf)
}

