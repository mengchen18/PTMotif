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
#' @param max.pvalue the maximum FDR to be reported
#' @param max.gw.pvalue the maximum p value of on gene wise to be reported, only useful
#'   when both \code{fg.genes} and \code{bg.genes} are specified
#' @param center the amino acid centered at the sequences, sequences with other center AA would be 
#'   removed from the list. To disable this function, set \code{center = NULL}. Only used in 
#'   \code{checkSeqs} and \code{motif_all}.
#' @param verbose logical, whether print message
#' @param annotate whether the discovered motif should be annotated by known motifs
#' @details more details here
#' @return a \code{data.frame} of the over-representation data analysis
#' @importFrom parallel mclapply
#' @importFrom stringr str_detect
#' @export
#'

detectMotifs <- function(fg.seqs, bg.seqs, fg.genes=NULL, bg.genes=NULL,
                         method = c("cvp", "all")[1], min.seqs = 5,
                         max.pvalue = 1e-6, max.gw.pvalue = NULL, center = "STY",
                         ncores = 1, verbose = TRUE, annotate = FALSE) {
  
  method <- match.arg(method, c("cvp", "all"), several.ok = TRUE)
  
  # check for fg./bg.genes, length, NA, avaliability, etc
  if (is.null(bg.genes) && !is.null(fg.genes)) {
    message("bg.genes is NULL, fg.genes is ignored. ")
    fg.genes <- NULL
  }
  if (is.null(fg.genes) && !is.null(bg.genes)) {
    message("fg.genes is NULL, bg.genes is ignored. ")
    bg.genes <- NULL
  }
  if (!is.null(fg.genes) && !is.null(bg.genes)) {
    if (length(fg.seqs) != length(fg.genes))
      stop("fg.genes should have the same length as fg.seqs")
    if (length(bg.seqs) != length(bg.genes))
      stop("bg.genes should have the same length as bg.seqs")
    if (any(is.na(fg.genes)))
      fg.genes[is.na(fg.genes)] <- "unnamedGene"
    if (any(is.na(bg.genes)))
      bg.genes[is.na(bg.genes)] <- "unnamedGene"
  }
  
  if (verbose)
    message("Checking input sequences ...")
  seqs <- checkSeqs(fg.seqs, bg.seqs, option = "extend", center = center)
  
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
    motif_all <- motif_all(seqs$fg.seqs, min.seqs=min.seqs, ncores = ncores, 
                           verbose=verbose, genes = fg.genes, center = center)
  }
  
  im <- setdiff(names(motif_cvp$mw), names(motif_all$mw))
  motifs <- c(motif_cvp$mw[im], motif_all$mw)
  
  if (length(motifs) == 0)
    return()
  
  if (verbose)
    message("Evaluating the significance of detected motifs ... ")
  
  r <- motifor(fg.count = motifs, n.fg.seqs = length(seqs$fg.seqs),
               bg.seqs = seqs$bg.seqs, ncores = ncores, max.pvalue = max.pvalue)
  
  if (annotate) {
    if (verbose)
      message("Annotating detected motifs ...")
    r <- annotKnownMotifs(r)
  }
  
  l <- NULL
  gw.keeprow <- NULL
  if (!is.null(bg.genes) && !is.null(fg.genes)) {
    if (verbose)
      message("Evaluating the significance of detected motifs on gene level ... ")
    
    if (nrow(r) == 0) {
      l <- r[, c("OR", "pvalue", "fg.count", "fg.total", "bg.count", "bg.total")]
    } else {
      gw <- c(motif_cvp$gw$motif.gene[im], motif_all$gw$motif.gene)
      fg.genes.unique <- unique(c(motif_all$gw$fg.gene, motif_cvp$gw$fg.gene))
      
      l <- mclapply(r$motif, function(m) {
        motifgeneor(motif = m, fg.genes = fg.genes.unique, fg.genes.motif = gw[[m]],
                    bg.seqs = bg.seqs, bg.genes = bg.genes)
      }, mc.cores = ncores)
      l <- do.call(rbind, l)
      l <- l[, -1]
      
      if (!is.null(max.gw.pvalue))
        gw.keeprow <- l$pvalue < max.gw.pvalue
    }
    colnames(l) <- paste("gw", colnames(l), sep = ".")
  }
  
  if (!is.null(l)) {
    r <- cbind(r, l)
    if (!is.null(gw.keeprow))
      r <- r[gw.keeprow, ]
  }
  
  return(r)
}
