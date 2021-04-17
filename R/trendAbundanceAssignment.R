#' Check Trend in Abundance and Taxonomic Assignment
#'
#' @name trendAbundanceAssignment
#'
#' @details Check if the more abundance-prevalent taxa
#'          have better taxonomic assignments. This is
#'          a pre-check, not corrected for differences
#'          in sequencing depth.
#'
#' @param x A phyloseq object
#'
#' @param quantiles Abundances values to sort. Can be changed to specify
#'                  how many values in a distribution are above or below
#'                  a certain limit.
#'
#' @param plot Logical. Default is TRUE.
#'
#' @examples
#' library(biomeUtils)
#' data("SprockettTHData")
#'
#' p1 <- trendAbundanceAssignment(SprockettTHData,
#'                               quantiles = seq(0, 95, by = 5),
#'                               plot=TRUE)
#' p1 + ggplot2::scale_colour_brewer("", palette = "Spectral") +
#'     ggplot2::theme_minimal() +
#'     ggplot2::theme(axis.text.x = element_text(angle=90, vjust=0.5))
#'
#' @return Either a list with data.frame and plot or just
#'         ggplot object.
#'
#' @author Sudarshan A. Shetty
#'
#' @references
#' \itemize{
#' \item{}{Thorsten Brach, 2018. MicrobiomeX2 Pipeline.
#'        \emph{GitHub}.
#'        \url{https://github.com/TBrach/MicrobiomeX2}
#'        \url{https://github.com/TBrach/Dada_Pipel/blob/master/Generalized_Phyloseq_Analysis_New.Rmd}
#' }
#' }
#' @importFrom stats quantile
#' @importFrom phyloseq taxa_sums 'sample_data<-'
#' @importFrom data.table as.data.table melt
#' @importFrom tidyr separate
#' @importFrom ggplot2 ggplot geom_point scale_x_discrete ylab xlab theme element_text
#' @export

trendAbundanceAssignment <- function(x, quantiles = seq(0, 95, by = 10),
                                     plot=TRUE){

  #globalvars
  Quant <- Var1 <- Var2 <- variable <- aes <- value <- NULL
  taxadf <- tax_table(x) %>%
    as.matrix() %>%
    as.data.frame(stringsAsFactors = FALSE)

  taxa.counts <- phyloseq::taxa_sums(x)

  abQuantiles <- stats::quantile(taxa.counts, probs = quantiles/100)

  which_taxa <- lapply(abQuantiles,
                       function(quant) {taxa.counts >= quant})

  numb_taxa <- sapply(which_taxa, sum)

  filtered.taxas <- lapply(which_taxa, function(indexes){
    taxadf[indexes,]
  })

  assignment_distributions <- lapply(filtered.taxas, .summarizeTaxonomicAssignmentsDF)

  pc_assigned <- sapply(assignment_distributions, function(distri){distri[["PC_assigned"]]})

  rownames(pc_assigned) <- colnames(taxadf)

  colnames(pc_assigned) <- paste("Ab_", quantiles, "_", round(abQuantiles), "_", numb_taxa, sep = "")
  pc_assigned <- as.data.frame(pc_assigned)
  pc_assigned$Var1 <- rownames(pc_assigned)
  pc_assigned <- data.table::as.data.table(pc_assigned)

  pc_assigned_ldf <- data.table::melt(pc_assigned,
                                      id.vars = c("Var1"),
                                      value.name = "value") %>%
    tidyr::separate(variable, c("Type", "Quant", "abQuant", "numb_taxa"), sep="_")

  pc_assigned_ldf$Var1 <- factor(pc_assigned_ldf$Var1, levels = colnames(taxadf), ordered = TRUE)

  pc_assigned_ldf$Quant <- factor(pc_assigned_ldf$Quant, levels = as.character(quantiles), ordered = TRUE)

  trendplot <- ggplot2::ggplot(pc_assigned_ldf, aes(x = Quant, y = value, col = Var1))
  trendplot <- trendplot +
    ggplot2::geom_point(size = 2) +
    ggplot2::scale_x_discrete(breaks = quantiles, labels = paste(quantiles, " (", numb_taxa, ")", sep = "")) +
    ggplot2::ylab("Percentage of taxa assigned") +
    ggplot2::xlab("Total counts quantile (No. of remaining taxa)") +
    ggplot2::theme(axis.text.x = element_text(angle=90, vjust=0.5))

  if(plot){
    return(trendplot)
  }
  return(list(pc_assigned_ldf, trendplot))

}


## Helper

.summarizeTaxonomicAssignmentsDF <- function(taxa){

  # taxa <- as.data.frame(unclass(tax_table(physeq)))

  countNA <- apply(taxa, 2, function(x){sum(is.na(x))})
  countNonNA <- apply(taxa, 2, function(x){sum(!is.na(x))})
  # ambiguous <- apply(taxa, 2, function(x){sum(grepl("/", x))})
  total <- countNA + countNonNA
  assignment_distribution <- data.frame(assigned = countNonNA,
                                        total = total,
                                        PC_assigned = round(100*countNonNA/total, 1))

}






