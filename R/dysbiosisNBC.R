#' Neighborhood Classification Dysbiosis Scores
#'
#' @name dysbiosisScore
#'
#' @details Calculates neighborhood classification based microbiome dysbiosis scores.
#'
#' Dysbiosis score or index are often used to identify the extent to which
#' clinically classified diseased samples are different from non diseased
#' reference samples. There are several approaches large-scale bacterial marker
#' profiling, relevant taxon-based methods, neighborhood classification,
#' random forest prediction, and combined alpha-beta diversity reviewed by Wei 2021.
#' Here we provide two scores based on neighborhood classification.
#'
#' \itemize{
#'
#' \item{'median-BC' }{The median Bray-Curtis distance between the test
#' sample and the reference samples (Lloyd-Price J, Arze C, Ananthakrishnan AN
#' et al. 2021).}
#'
#' \item{'BC-ED' }{The Euclidian distance between test sample and the centroid (median)
#' of PCoA Axis 1 and 2 (based on Bray-Curtis dissimilarity) of the reference
#' samples. (Mahapatra S et al 2021)}
#'
#' }
#'
#' @param x A phyloseq object
#'
#' @param reference_samples Vector of samples to use as reference.
#'
#' @param method Either 'median-BC' or 'BC-ED'
#'
#' @return A tibble data frame with dysbiosis score and sample information
#'
#' @examples
#' library(biomeUtils)
#' data("FuentesIliGutData")
#' # Keep only two groups. controls and L1.
#' ps <- filterSampleData(FuentesIliGutData, ILI != "L2")
#' ps <- getProportions(ps)
#' # Define controls as reference samples
#' ref_samples <- rownames(meta(subset_samples(ps, ILI == "C")))
#'
#' dys_bc <- dysbiosisScore(ps,
#'                          reference_samples = ref_samples,
#'                          method="median-BC")
#' dys_bc_ed <- dysbiosisScore(ps,
#'                             reference_samples = ref_samples,
#'                             method="BC-ED")
#'
#' @references
#' \itemize{
#' \item{}{Lloyd-Price J, Arze C, Ananthakrishnan AN et al. (2019).
#' Multi-omics of the gut microbial ecosystem in inflammatory bowel diseases.
#' \emph{Nature}, 569(7758), pp.655-662.}
#'
#' \item{}{Mahapatra S et al. (2021). Nanopore 16S rRNA sequencing reveals
#' alterations in nasopharyngeal microbiome and enrichment of Mycobacterium
#' and Mycoplasma in patients with COVID 19.
#' \emph{medRxiv preprint}, 2021.11.10.21266147.}
#'
#' }
#'
#' @author Sudarshan A. Shetty
#'
#' @importFrom dplyr left_join rename filter summarise group_by mutate
#' @importFrom stats median
#' @export
dysbiosisScore <- function(x,
                           reference_samples = NULL,
                           method=c("median-BC", "BC-ED")){

  S1 <- S2 <- n <- dysbiosis.score <- value <- pcoa.tab <- dys.tib <- NULL
  dist.melt.sample <- sam_tib <- centroids <- hold <- NULL
  cent.1 <- cent.2 <- a.two <- a.one <- is.dysbiotic <- NULL
  Axis.1 <- Axis.2 <- S2 <- dysbiosis.score <- NULL

  if(method=="median-BC"){

    if(is.null(reference_samples)){
      stop("Provide a reference set of samples")
    }

    dist_mat <- phyloseq::distance(x, "bray")

    dist.melt.sample <- meltDistanceToTable(x,
                                            dist_mat = dist_mat)

    # check
    dys.tib <- dist.melt.sample |>
      dplyr::filter(S1 %in% reference_samples) |>
      dplyr::group_by(S2) |>
      dplyr::summarise(dysbiosis.score =median(value))

    sam_tib <- getSampleTibble(x)

    dys.tib <- dys.tib |>
      dplyr::left_join(sam_tib, by = c("S2" = colnames(sam_tib)[1]))

    return(dys.tib)
  }
  if(method =="BC-ED"){
    pcoa.tab <- phyloseq::plot_ordination(x,
                                          ordination =
                                            phyloseq::ordinate(x, "PCoA", "bray"),
                                          justDF = T)
    centroids <- pcoa.tab |>
      dplyr::mutate(hold = rownames(pcoa.tab)) |>
      dplyr::filter(hold %in% reference_samples) |>
      dplyr::summarise(cent.1 = median(Axis.1),
                       cent.2 = median(Axis.2))
    pcoa.tab <- pcoa.tab |>
      dplyr::mutate(a.one = (Axis.1 - centroids[[1]]),
                    a.two = (Axis.2 - centroids[[2]])) |>
      dplyr::mutate(dysbiosis.score = sqrt((a.one)^2 + (a.two)^2)) |>
      dplyr::select(-c(a.one, a.two))
    return(pcoa.tab)
  }

}

# Helper classification
#' @importFrom stats quantile
.is_dysbiotic <- function(dys.tib, reference_samples, percentile=.90){
  is.dysbiotic <- S2 <- dysbiosis.score <- NULL
  ref.dys <- dys.tib |>
    dplyr::filter(S2 %in% reference_samples) |>
    dplyr::pull(dysbiosis.score)
  dys.tib <- dys.tib |>
    dplyr::mutate(is.dysbiotic = ifelse(dysbiosis.score >= quantile(ref.dys,probs = percentile),
                                        "yes", "no"))
}

