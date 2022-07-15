#' Cluster ASVs into Clusters
#'
#' @name clusterASVs
#'
#' @details clusterASVs clusters ASV sequences
#'          obtained after \code{dada2} in \code{phyloseq} object into clusters
#'          based on user provided sequence identity. If using method="inexact"
#'          this is similar to OTU approach for example clustering of ASVs that
#'          are 97% similar into one operational taxonomic unit. While preference
#'          for ASVs or OTUs is highly debated (see reference section), both
#'          have their pros and cons.The user is advised to make appropriate
#'          choice based on their research questions and properties of data.
#'
#' @param x A phyloseq object
#'
#' @param cluster_method See methods in \code{\link[DECIPHER]{IdClusters}}
#'                       in \code{DECIPHER}
#'                       Default is "inexact"
#'
#' @param cluster_cutoff Cut-off for clustering ASVs.
#'                       Default= 0.01 i.e. 99 percent identity
#'
#' @param threads No. of processors to use
#'
#' @examples
#' \dontrun{
#' library(biomeUtils)
#' library(dplyr)
#' library(data.table)
#' library(tibble)
#'
#' testSprockett <- core(subset_samples(SprockettTHData, Community == "River_4"), 5, 20 / 100)
#' testSprockett <- clusterASVs(testSprockett,
#'                              cluster_cutoff = 0.01,
#'                              threads = 2)
#'}
#' @return A new \code{phyloseq} object with cluster ASVs
#'
#' @author Sudarshan A. Shetty
#'
#' @references
#' \itemize{
#' \item{}{Clustering ASVs after DADA2? generating new OTU table.
#'        \emph{GitHub issue #947:}
#'        \url{https://github.com/benjjneb/dada2/issues/947}
#' }
#' \item{}{Michael McLaren. mikemc/speedyseq: speedyseq v0.2.0.
#'         \emph{Zenodo}, 2020, June 30, (Version v0.2.0).
#'         \url{http://doi.org/10.5281/zenodo.3923184}
#' }
#' \item{}{PD Schloss. Pre-print Review of 97% Identity Threshold Manuscript.
#'         \emph{Blogpost: The Academic Hermit}, 2017, Oct 11.
#'         \url{http://www.academichermit.com/2017/10/11/Review-of-Updating-the-97-identity-threshold.html}
#' }
#' \item{}{Angela Oliverio and Noah Fierer. Intragenomic heterogeneity and its
#'         implications for ESVs.
#'         \emph{Blogpost: Fierer Lab}, 2017, Oct 9.
#'         \url{http://fiererlab.org/2017/10/09/intragenomic-heterogeneity-and-its-implications-for-esvs/}
#' }
#' \item{}{Noah Fierer, Tess Brewer, and Mallory Choudoir. Lumping versus
#'         splitting – is it time for microbial ecologists to abandon OTUs?.
#'         \emph{Blogpost: Fierer Lab}, 2017, May 2.
#'         \url{http://fiererlab.org/2017/05/02/lumping-versus-splitting-is-it-time-for-microbial-ecologists-to-abandon-otus/}
#' }
#' }
#'
#'
#' @importFrom Biostrings DNAStringSet
#' @importFrom DECIPHER AlignSeqs DistanceMatrix IdClusters
#' @importFrom data.table data.table
#' @importFrom phyloseq taxa_names otu_table otu_table<-
#' @importFrom microbiome abundances
#'
#' @export
NULL

clusterASVs <- function(x,
                        cluster_method = "inexact",
                        cluster_cutoff = 0.03,
                        threads = 2) {

  .Deprecated("clusterASVs")
  # global vars
  #utils::globalVariables(c("seqs", "aln", "d.mat", "group", "taxa", "clusters", "."))
  seqs <- aln <- d.mat <- group <- taxa <- clusters <- . <- NULL

  .check_refseq(x)

  if (!is.null(x@refseq)) {
    seqs <- Biostrings::DNAStringSet(refseq(x))
    if (!is.null(names(seqs))) {
      names(seqs) <- NULL
      taxa_names(x) <- seqs
    }

    if(cluster_method=="inexact"){
      clusters.mat <- DECIPHER::IdClusters(myXStringSet = seqs,
                                           method = cluster_method,
                                           cutoff = cluster_cutoff, # cluster.cutoff=0.03
                                           processors = threads)
    } else {
      aln <- DECIPHER::AlignSeqs(seqs, processors = threads)
      d.mat <- DECIPHER::DistanceMatrix(aln, processors = threads)
      clusters.mat <- DECIPHER::IdClusters(d.mat,
                                           #method = cluster_method,
                                           cutoff = cluster_cutoff, # cluster.cutoff=0.03
                                           processors = threads)
    }

    new_names <- data.table::data.table(taxa = taxa_names(x),
                                        prevalence = calculateQC(x)$TaxaQC$taxa_prevalence,
                                        taxa_counts = calculateQC(x)$TaxaQC$taxa_counts,
                                        taxa_mean_counts = calculateQC(x)$TaxaQC$taxa_mean_counts,
                                        group = clusters.mat$cluster)
    new_names <- new_names[, by = group, .(cASV = taxa[which.max(prevalence)])]

    otu <- otu_table(rowsum(microbiome::abundances(x),
                            clusters.mat$cluster,
                            reorder = F),
                     taxa_are_rows = TRUE)

    stopifnot(all.equal(as.character(new_names$group), taxa_names(otu)))
    taxa_names(otu) <- new_names$cASV


  }

  if(is.null(x@phy_tree)){
    #print("yes")
    x <- merge_phyloseq(otu,
                        tax_table(x),
                        sample_data(x),
                        refseq(x))

  } else {
    x <- merge_phyloseq(otu,
                        tax_table(x),
                        sample_data(x),
                        phy_tree(x),
                        refseq(x))
  }


  return(x)
}

# table(clusters$cluster)
.check_refseq <- function(x) {
  if (is.null(x@refseq)) {
    stop("Please provide a valid input
            Either the phyloseq with sequences in refseq slot
         OR Sequences as DNAStringSet")
  }
}
