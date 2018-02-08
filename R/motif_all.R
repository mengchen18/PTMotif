
#' Find ALL motifs in a list of sequences using motif-all algorithm
#' @author Chen Meng
#' @param seqs a character vector of sequences
#' @param min.seqs minimum number of motif
#' @param ncores the number of cores to be used, passed to \code{mclapply}.
#' @param verbose logical, whether print detailed information
#' @references He, Zengyou, Can Yang, Guangyu Guo, Ning Li, and Weichuan Yu. 2011.
#'   "Motif-All: Discovering All Phosphorylation Motifs." BMC Bioinformatics
#'   12 Suppl 1 (February):S22.
#' @return a names integer vector. The names are the sequence motifs and
#'   the integers indicate the frequency of each motifs in the input sequences.
#' @importFrom parallel mclapply
#' @importFrom stringr str_detect
#' @export
#'

motif_all <- function(seqs, min.seqs, ncores = 1, verbose = FALSE) {

  # AA letters
  aaa <- c("A","C","D","E","F","G","H","I","K","L",
           "M","N","P","Q","R","S","T","V","W","Y")

  # def function
  toAAFreq <- function(x, aaa, n, excludeCol = NULL) {
    # x - a matrix
    # amino acid character vector
    # minimum number of frequency
    # excludeCol
    freqmat <- apply(x, 2, function(x) table(x)[aaa])
    icc <- freqmat >= n
    i <- which(icc, arr.ind = TRUE)
    df <- data.frame(AA = aaa[i[, 1]], col = i[, 2], freq = freqmat[which(icc)],
                     stringsAsFactors = FALSE)
    if (!is.null(df))
      df <- df[!df$col %in% excludeCol, ]
    if (nrow(df) == 0)
      return(list(df))
    split(df, 1:nrow(df))
  }

  #
  seq0 <- seqs
  seqs <- strsplit(seqs, "|")
  seqmat <- do.call(rbind, seqs)
  # seqmat[seqmat == "_"] <- NA
  nc <- ncol(seqmat)

  ##
  ll <- list()
  v <- toAAFreq(seqmat, aaa = aaa, n = min.seqs)
  if (nrow(v[[1]]) == 0) {
    message(paste("No motif presents more than", min.seqs, "times."))
    return()
  }
  ll[[1]] <- v
  len <- ncol(seqmat)
  st <- rep(TRUE, nrow(seqmat))
  for (i in 2:nc) {
    if (verbose)
      cat(paste("Trying to find motifs having", i, "amino acids ...\n"))
    v <- mclapply(v, function(x) {
      rindex <- st
      for (ir in 1:nrow(x)) {
        rindex <- rindex & seqmat[, x$col[ir]] == x$AA[ir]
      }
      smat <- seqmat[rindex, ]
      vv <- toAAFreq(smat, aaa = aaa, n = min.seqs, excludeCol = 1:max(x$col))
      lapply(vv, rbind, x)
    }, mc.cores = ncores)
    v <- unlist(v, recursive = FALSE, use.names = FALSE)
    v <- setdiff(unique(v), ll[[i-1]])
    if (length(v) == 0)
      break()
    ll[[i]] <- v
  }
  res <- unlist(ll[-1], recursive = FALSE)
  
  if (length(res) == 0) {
    message(paste("No motif presents more than", min.seqs, "times."))
    return()
  }
  if (verbose)
    cat("Summarizing results...\n")
  sv <- rep(".", nc)
  res <- sapply(res, function(x) {
    pp <- sv
    pp[x$col] <- x$AA
    c(paste(pp, collapse = ""), x$freq[1])
  })
  res <- structure(as.integer(res[2, ]), names = res[1, ])
  sort(res, decreasing = TRUE)
}

