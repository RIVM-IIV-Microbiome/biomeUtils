#' Read Merged Metaphlan Table
#'
#' @name readMergedMetaphlan
#'
#' @details This function reads Metaphlan3 output that is merged with all samples.
#'          and returns a \code{phyloseq} object with otu_table and tax_table.
#'
#' Note:
#' \code{readMergedMetaphlan} is tested for Metaphlan version 3 and based on the RIVM
#' HPC HUMANN3 pipeline output built by Jeron Frank.
#'
#' @param input_file_path path to where the output of Metaphlan3 is located.
#'
#' @param find_sample_name_pattern Naming convention used for samples within
#'                                 the input file.
#'
#' @param replace_sample_name_pattern specific pattern to replace in sample names.
#'
#' @param return_species Logical whether to return only species level data. Default TRUE.
#'
#' @examples
#' \dontrun{
#'
#' library(biomeUtils)
#'
#' mpa.ps <- readMergedMetaphlan(input_file_path = "../metaphlan_bugs_list_merged.txt",
#'                              find_sample_name_pattern = "_kneaddata_concat_metaphlan_bugs_list",
#'                              replace_sample_name_pattern = "",
#'                              return_species = TRUE)
#' }
#'
#' @return A \code{phyloseq} object
#'
#' @author Sudarshan A. Shetty
#'
#' @references
#' Shetty SA (2020). Utilities for microbiome analytics.
#' \url{https://github.com/microsud/biomeUtils}
#'
#' @importFrom data.table :=
#' @export

readMergedMetaphlan <- function(input_file_path = NULL,
                                find_sample_name_pattern = NULL,
                                replace_sample_name_pattern = "",
                                return_species = TRUE) {

  V1 <- clade_name <- NULL

  if(is.null(input_file_path)){
    stop("Please specific input_file_path ")
  }

  if(is.null(find_sample_name_pattern)){
    stop("Please specific find_sample_name_pattern")
  }

  mpa_tab <- data.table::fread(input_file_path)
  #head(mpa_tab)
  colnames(mpa_tab) <- gsub(find_sample_name_pattern,
                            replace_sample_name_pattern,
                            colnames(mpa_tab))

  # Get lowest taxonomic classification.
  clde.nam <- mpa_tab[,clade_name]

  shortnames <- gsub(paste0(".+\\", "|"), "", clde.nam)
  rownames(mpa_tab) <- shortnames

  tax.lineage <- data.table::as.data.table(mpa_tab$clade_name)
  tax.lineage <- tax.lineage[, c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
                             := data.table::tstrsplit(V1, "|", fixed=T)][]

  tax.lineage <- as.matrix(tax.lineage)

  rownames(tax.lineage) <- gsub(paste0(".+\\", "|"), "", tax.lineage[,1])

  tax.lineage <- tax.lineage[,-1]
  taxmat <- phyloseq::tax_table(as.matrix(tax.lineage))

  otu.tb <- mpa_tab[,-c("clade_name", "clade_taxid")]
  otu.tb[is.na(otu.tb)] <- 0
  rownames(otu.tb) <- shortnames
  #rownames(otu.tb)[1:10]
  otu.tb <- phyloseq::otu_table(as.matrix(otu.tb), taxa_are_rows = T)
  rownames(otu.tb) <- shortnames
  res <- phyloseq::phyloseq(taxmat, otu.tb)
  return(res)
}
