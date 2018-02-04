#' Find fixed position motifs (inconsecutive or consecutive) in a list of
#'   sequences using an exhausitive approach
#' @author Chen Meng
#' @param seqs a character vector of sequences
#' @param min.seqs minimum number of motif
#' @param ncores the number of cores to be used, passed to \code{mclapply}.
#' @note the algorithm uses an exhaustive approach to find all (long)
#'   fixed position motifs, so it could be slow when the number of sequences is large.
#' @return a names integer vector. The names are the sequence motifs and
#'   the integers indicate the frequency of each motifs in the input sequences.
#' @importFrom parallel mclapply
#' @importFrom stringr str_detect
#' @export
#'
motif_exh <- function(seqs, min.seqs=1, ncores = 1) {

  if (length(seqs) < min.seqs) {
    warning("length of seqs < min.seqs")
    return()
  }

  seq0 <- seqs
  seqs <- strsplit(seqs, "|")
  dots <- rep(".", length(seqs[[1]]))
  ls <- combn(seq_along(seqs), m = 2)
  motifs <- mclapply(1:ncol(ls), function(x, mot) {
    i <- seqs[[ls[1, x]]] == seqs[[ls[2, x]]] & seqs[[ls[1, x]]] != "_"
    if (sum(i) <= 1)
      return()
    mot[i] <- seqs[[ls[1, x]]][i]
    paste0(mot, collapse = "")
  }, mot = dots, mc.cores = ncores)
  motifs <- unique(unlist(motifs))

  count <- mclapply(motifs, function(x) sum(str_detect(seq0, x)),
                    mc.cores = ncores)
  count <- unlist(count)
  names(count) <- motifs
  count <- count[count >= min.seqs]
  sort(count, decreasing = TRUE)
}
