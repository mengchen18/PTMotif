#' Find consecutive but variable position (CVP) motif in a list of sequences
#' @author Chen Meng
#' @param seqs a character vector of sequences
#' @param min.seqs minimum number of motif
#' @param genes the gene name from where a sequences is discovered, a character vector
#'   has the same length as \code{seqs}
#' @param ncores the number of cores to be used, passed to \code{mclapply}.
#' @note the algorithm uses an exhaustive approach to find all consecutive
#'   motifs, so it could be slow when the number of sequences is large.
#' @return a list consists of:
#'   $mw - motif wise count, a names integer vector. The names are the sequence motifs and
#'   the integers indicate the frequency of each motifs in the input sequences.
#'   $gw - gene wise count, a list of two elements: 1) unique gene and 2) motif genes, i.e.
#'   genes include a specific type of motif
##' @importFrom parallel mclapply
#' @examples
#'   seqs <- c(paste(LETTERS[1:6], collapse = ""),
#'   paste(LETTERS[2:7], collapse = ""),
#'   paste(LETTERS[3:8], collapse = ""),
#'   paste(LETTERS[4:9], collapse = ""),
#'   paste(LETTERS[5:10], collapse = ""))
#'   motif_cvp(seqs)
#' @export
#'
motif_cvp <- function(seqs, min.seqs=1, genes=NULL, ncores = 1) {

  if (length(seqs) < min.seqs) {
    warning("length of seqs < min.seqs")
    return()
  }
  
  seq0 <- seqs
  seqs <- strsplit(seqs, "|")
  aseq <- mclapply(seqs, subSeqs, mc.cores = ncores, trim.char=TRUE)
  aseq <- unlist(aseq)
  count <- c(table(aseq))
  count <- count[count >= min.seqs]
  res <- sort(count, decreasing = TRUE)
  
  if (!is.null(genes)) {
    gn <- lapply(names(res), function(x) 
      unique(genes[grep(x, seq0)])
    )
    names(gn) <- names(res)
    gw <- list(motif.gene = gn, fg.gene = unique(genes))
  } else 
    gw <- NULL
  
  list(mw = res, gw = gw)
}
