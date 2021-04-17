#' Calculate QC Metrics for Taxa and Samples
#'
#' @name calculateQC
#'
#' @details Calculate QC metrics for taxa and samples.
#'
#' @param x A phyloseq object
#'
#' @return A list with two tibbles.
#'
#' \emph{SampleQC}
#' \itemize{
#' \item{sample_total_taxa:}{Total taxa in a sample}
#' \item{sample_total_reads:}{Total reads in a sample}
#' \item{number_taxa_account_fifty_percent:}{Number of taxa account for 50 percent}
#' }
#' \emph{TaxaQC}
#' \itemize{
#' \item{taxa_counts:}{Total counts of taxa in all samples}
#' \item{taxa_mean_counts:}{Mean counts of taxa in all samples}
#' \item{taxa_dectected_samples:}{Number of samples in which taxa detected}
#' \item{taxa_prevalence:}{Percent of samples in which taxa detected}
#' }
#'
#' @author Sudarshan A. Shetty
#'
#' @examples
#' library(biomeUtils)
#' data("FuentesIliGutData")
#' qc_Data <- calculateQC(FuentesIliGutData)
#' qc_Data
#'
#' @importFrom microbiome prevalence abundances
#' @importFrom phyloseq nsamples ntaxa sample_sums taxa_sums taxa_sums
#' @importFrom phyloseq merge_phyloseq phy_tree taxa_names<-
#' @importFrom tibble tibble
#' @importFrom stats sd
#'
#' @export

calculateQC <-  function(x) {

  # gloabl vars
  group <- merge_phyloseq <- taxa <-taxa_counts <- taxa_mean_counts <- NULL
  taxa_stdev_counts <- NULL

  ## Must be a phyloseq object
  if ( !is(x, "phyloseq") ){
    stop("input must be an phyloseq object.")
  }

  ## the input must have few samples
  if ( nsamples(x) < 1 ){
    stop("input must have at least one sample")
  }

  if ( ntaxa(x) < 1 ) {
    stop("input must have at least one taxa")
  }

  ## See what versions of the expression data are available in the object
  if(any(sample_sums(x)==1) | any(.check_decimal(sample_sums(x)))){
    stop("Data must be counts")
  }

  # Create a tibbles with QC information
  qcSample <- tibble::tibble(Samples=sample_names(x),
                             # number of taxa per sample
                             sample_total_taxa = colSums(microbiome::abundances(x) != 0),
                             # Total reads/sample
                             sample_total_reads = sample_sums(x),
                             number_taxa_account_fifty_percent = microbiome::coverage(x))

  qcTaxa <- tibble::tibble(taxa = taxa_names(x),
                           # number of taxa per sample
                           taxa_counts = taxa_sums(x),
                           taxa_mean_counts = rowMeans(microbiome::abundances(x)),
                           taxa_stdev_counts = apply(microbiome::abundances(x), 1, sd),
                           taxa_cv = taxa_stdev_counts/taxa_mean_counts,
                           taxa_rank = rank(rowMeans(microbiome::abundances(x))),
                           taxa_dectected_samples = rowSums(microbiome::abundances(x) != 0),
                           taxa_prevalence = microbiome::prevalence(x) *100,
                           percent_of_total = taxa_counts/sum(sample_sums(x)) *100)

  return(list(SampleQC=qcSample, TaxaQC=qcTaxa))


}


#' @keywords internal
#' check decimal
#' \url{https://www.rdocumentation.org/packages/schoolmath/versions/0.4/topics/is.decimal}
.check_decimal <- function(x){

  start <- 1
  end <- length(x)+1
  while(start<end){
    y <- x[start]

    test <- floor(y)
    if(y==test){
      if(start==1){
        result=FALSE
      }else{
        result<- c(result,FALSE)
      }

    }else{
      if(start==1){
        result=TRUE
      }else{
        result <- c(result,TRUE)
      }
    }
    start <- start+1
  }

  return(result)
}


