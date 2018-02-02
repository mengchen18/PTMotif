#' Find consecutive but variable position (CVP) motif in a list of sequences
#' @author Chen Meng
#' @param seqs a character vector of sequences
#' @param min.seqs minimum number of motif
#' @param ncores the number of cores to be used, passed to \code{mclapply}.
#' @note the algorithm uses an exhaustive approach to find all consecutive
#'   motifs, so it could be slow when the number of sequences is large.
#' @return a names integer vector. The names are the sequence motifs and
#'   the integers indicate the frequency of each motifs in the input sequences.
#' @importFrom parallel mclapply
#' @examples
#'   seqs <- c(paste(LETTERS[1:6], collapse = ""),
#'   paste(LETTERS[2:7], collapse = ""),
#'   paste(LETTERS[3:8], collapse = ""),
#'   paste(LETTERS[4:9], collapse = ""),
#'   paste(LETTERS[5:10], collapse = ""))
#'   motif_cvp(seqs)
#' @export
#'
motif_cvp <- function(seqs, min.seqs=1, ncores = 1) {

  if (length(seqs) < min.seqs) {
    warning("length of seqs < min.seqs")
    return()
  }

  seqs <- strsplit(seqs, "|")
  aseq <- mclapply(seqs, subSeqs, mc.cores = ncores)
  aseq <- unlist(aseq)
  count <- c(table(aseq))
  count <- count[count >= min.seqs]
  sort(count, decreasing = TRUE)
}
