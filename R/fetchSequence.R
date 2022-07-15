#' Retrieve Sequence from Phyloseq
#'
#' @name fetchSequence
#'
#' @details Retrieve sequences from refseq slot in phyloseq
#'          \code{phyloseq} object.
#'
#' @param x A phyloseq object
#'
#' @param seq_id A list of IDs matching taxa_names(phyloseq) that also have
#'               sequence information stored in refseq slot of
#'               \code{phyloseq} object.
#'
#' @param as_DNAStringSet Logical (Default=TRUE). If FALSE sequences are
#'                        returned as character strings.
#'
#' @examples
#' library(biomeUtils)
#' data("SprockettTHData")
#' ps <- SprockettTHData
#' seq_id <- c("ASV_30_Lachnospiraceae", "ASV_451_UC5-1-2E3", "ASV_586_CAG-56")
#' myseqs <- fetchSequence(ps, seq_id = seq_id, as_DNAStringSet = TRUE)
#' myseqs
#' @return Sequences for user species IDs either as vector
#'         or DNAStringSet object
#'
#' @author Sudarshan A. Shetty
#'
#' @importFrom phyloseq refseq
#'
#' @export

fetchSequence <- function(x,
                          seq_id = NULL,
                          as_DNAStringSet = TRUE) {
  if (!is(x, "phyloseq")) {
    stop("Input is not an object of phyloseq class")
  }

  if (is.null(x)) {
    stop("refseq(pseq) is NULL. Please provide reference sequences
         in phyloseq refseq slot.")
  }

  seqs <- refseq(x)

  # names(seqs)[1:4]
  list.df <- data.frame(seq_id = seq_id)

  user.seqs <- seqs[list.df$seq_id]

  if (as_DNAStringSet) {
    user.seqs <- seqs[list.df$seq_id]

    return(user.seqs)
  } else {
    return(as.vector(user.seqs))
  }
}
