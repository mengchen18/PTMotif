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
#' @importFrom stats p.adjust phyper
#' @export
#'

motifor <- function(fg.count, n.fg.seqs, bg.seqs, ncores = 1, max.fdr = 1e-2) {

  motifs <- names(fg.count)
  count <- mclapply(motifs, function(x) sum(str_detect(bg.seqs, x)),
                    mc.cores = ncores)
  count <- unlist(count)

  nbg <- length(bg.seqs)
  
  pv <- phyper(q = fg.count-1, m = count, n = nbg-count,
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
                   bg.total = nbg,
                   stringsAsFactors = FALSE,
                   row.names = NULL)
  df <- df[df$FDR < max.fdr, ]
  df[order(df$FDR, decreasing = FALSE), ]
}



#' Find over-represented motifs in a foreground list using hypergeometric test
#' @author Chen Meng
#' @param motif a motif
#' @param fg.genes unique foreground genes
#' @param fg.genes.motif unique foreground genes containing the motif
#' @param bg.seqs the background sequences
#' @param bg.genes a character vector of genes names of the bg.seqs, it
#'   should have the same length as \code{bg.seqs}
#' @importFrom parallel mclapply
#' @importFrom stringr str_detect
#' @importFrom stats p.adjust phyper
#' @export

motifgeneor <- function(motif, fg.genes, fg.genes.motif, bg.seqs, bg.genes)  {

  x <- length(fg.genes.motif)
  i <- str_detect(bg.seqs, motif)
  bg.genes.nomotif <- bg.genes[!i]
  n <- length(unique(bg.genes.nomotif))
  bg.genes.motif <- bg.genes[i]
  m <- length(unique(bg.genes.motif))
  k <- length(fg.genes)
  pv <- phyper(q = x-1, m = m, n = n, k = k, lower.tail = FALSE, log.p = FALSE)

  nbg <- m+n
  count <- m
  or <- (x/k)/(count/nbg)
  data.frame(motif = motif,
             OR = or,
             pvalue = pv,
             fg.count = x,
             fg.total = k,
             bg.count = count,
             bg.total = nbg,
             stringsAsFactors = FALSE,
             row.names = NULL)
  
}



