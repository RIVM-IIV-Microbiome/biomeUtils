#' Save Sequence
#'
#' @name saveSequence
#'
#' @details Saves sequences either from a DNAStringSet or refseq slot
#'          \code{phyloseq} object. If taxa_names is sequences, use
#'          'add_refseq' function from \code{microbiomeutilties} before
#'          using this function.
#'
#' @param x Either a phyloseq object with refseq information or
#'          a object of class DNAStringSet.
#'
#' @param file File name with extension ".fasta" or ".fa", e.g.
#'             'myseqs.fasta' or file path where you wish to save the
#'             output fasta e.g. 'mypath/myseqs.fasta'.
#'             As a fail safe, default is NULL which returns
#'             an error.
#' @examples
#'
#' # data("SprockettTHData")
#' # ps <- SprockettTHData
#' # seq_id <- c("ASV_30_Lachnospiraceae", "ASV_451_UC5-1-2E3", "ASV_586_CAG-56")
#'
#' # saveSequence(ps,
#' #             file="inst/extras/sequences.fasta")
#' @return A fasta file saved in user specified location
#'
#' @author Sudarshan A. Shetty
#'
#' @importFrom phyloseq refseq
#' @importFrom Biostrings writeXStringSet
#'
#' @export

saveSequence <- function(x, file = NULL) {
  if (is.null(file) || is.na(file)) {
    stop("Provide a file name e.g.
         'mypath/myseqs.fasta'")
  }
  if (is(x)[1] == "DNAStringSet") {
    writeXStringSet(x, file)
  } else {
    .check_refseq(x)
    seqs <- refseq(x)
    writeXStringSet(seqs, file)
  }
}


# when taxa_names is seqs
.check_refseq <- function(x) {
  if (is.null(refseq(x))) {
    stop("refseq() not found. Please provide a phyloseq with
         sequences stored in refseq()")
  }
}
