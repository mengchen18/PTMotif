#' Find over-represented motifs in a foreground list using hypergeometric test
#' @author Chen Meng
#' @param fg.count a names integer vector. The names are the sequence motifs and
#'   the integers indicate the frequency of each motifs in the input sequences,
#'   usually return by \code{\link{motif_all}} or \code{\link{motif_cvp}}.
#' @param n.fg.seqs the number of foreground sequences
#' @param bg.seqs a character vector of background sequences
#' @param ncores the number of cores to be used, passed to \code{mclapply}.
#' @param max.fdr the maximum FDR to be reported
#' @return a \code{data.frame} of the over-representation data analysis
#' @importFrom parallel mclapply
#' @export
#'

motifor <- function(fg.count, n.fg.seqs, bg.seqs, ncores = 1, max.fdr = 1e-2) {
  
  motifs <- names(fg.count)
  count <- mclapply(motifs, function(x) sum(str_detect(bg.seqs, x)),
                    mc.cores = ncores)
  count <- unlist(count)

  nbg <- length(bg.seqs)
  pv <- phyper(q = fg.count, m = count, n = nbg-count,
               k = n.fg.seqs, lower.tail = FALSE, log.p = FALSE)
  fdr <- p.adjust(pv, method = "fdr")
  or <- (fg.count/n.fg.seqs)/(count/nbg)
  df <- data.frame(motif = names(pv),
                   OR = or,
                   pvalue = pv,
                   FDR = fdr,
                   fg.count = fg.count,
                   fg.total = n.fg.seqs,
                   bg.count = count,
                   bg.totol = nbg,
                   stringsAsFactors = FALSE,
                   row.names = NULL)
  df <- df[df$FDR < max.fdr, ]
  df[order(df$FDR, decreasing = FALSE), ]
}
