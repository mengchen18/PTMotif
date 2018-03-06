#' Detect over-represented motifs
#' @author Chen Meng
#' @param fg.seqs a character vector of foreground sequences
#' @param bg.seqs a character vector of background sequences
#' @param fg.genes character vector of gene names of foreground sequences, it
#'   should have the same length as \code{fg.seqs}
#' @param bg.genes character vector of gene names of background sequences, it
#'   should have the same length as \code{bg.seqs}
#' @param method the algorithm used to detect potential motifs, multiple are allowed
#' @param min.seqs the minimum frequency of a motif should be considered
#' @param ncores the number of cores to be used, passed to \code{mclapply}.
#' @param max.fdr the maximum FDR to be reported
#' @param annotate whether the discovered motif should be annotated by known motifs
#' @details more details here
#' @return a \code{data.frame} of the over-representation data analysis
#' @importFrom parallel mclapply
#' @importFrom stringr str_detect
#' @export
#'

detectMotifs <- function(fg.seqs, bg.seqs, fg.genes=NULL, bg.genes=NULL,
                         method = c("cvp", "all")[1], min.seqs = 5,
                         max.fdr = 1e-2, ncores = 1, verbose = TRUE, annotate = FALSE) {

  method <- match.arg(method, c("cvp", "all"), several.ok = TRUE)

  # check for fg./bg.genes, length, NA, avaliability, etc

  if (is.null(bg.genes) && !is.null(fg.genes)) {
    message("bg.genes is NULL, fg.genes is ignored. ")
    fg.genes <- NULL
  }

  if (verbose)
    message("Checking input sequences ...")
  seqs <- checkSeqs(fg.seqs, bg.seqs, option = "extend")

  if (length(fg.seqs) < min.seqs) {
    message("No potential motif discovered, try to lower the min.seqs.")
    return (NULL)
  }

  motif_cvp <- motif_all <- NULL
  if ("cvp" %in% method) {
    if (verbose)
      message("Detecting potential motifs using 'cvp' algorithm ...")
    motif_cvp <- motif_cvp(seqs$fg.seqs, min.seqs=min.seqs, ncores = ncores, genes = fg.genes)
  }
  if ("all" %in% method) {
    if (verbose)
      message("Detecting potential motifs using 'all' algorithm ...")
    motif_all <- motif_all(seqs$fg.seqs, min.seqs=min.seqs, ncores = ncores, verbose=verbose, genes = fg.genes)
  }
  motifs <- unique(c(motif_cvp$mw, motif_all$mw))

  if (length(motifs) == 0)
    return()

  if (verbose)
    message("Evaluating the significance of detected motifs ... ")
  r <- motifor(fg.count = motifs, n.fg.seqs = length(seqs$fg.seqs),
          bg.seqs = seqs$bg.seqs, ncores = ncores, max.fdr = max.fdr)

  if (annotate) {
    if (verbose)
      message("Annotating detected motifs ...")
    r <- annotKnownMotifs(r)
  }

  if (!is.null(bg.genes)) {
    gw <- unique(c(motif_cvp$gw$motif.gene, motif_all$gw$motif.gene))
    fg.genes.unique <- gw$fg.gene

    l <- lapply(r$motif, function(m) {
      motifgeneor(motif = m, fg.genes = fg.genes.unique, fg.genes.motif = gw[[m]],
                  bg.seqs = bg.seqs, bg.genes = bg.genes)
    })
    l <- do.call(rbind, l)
  }


  cbibd(r, l)
}
