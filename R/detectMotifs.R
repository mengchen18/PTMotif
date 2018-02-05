#' Detect over-represented motifs
#' @author Chen Meng
#' @param fg.seqs a character vector of foreground sequences
#' @param bg.seqs a character vector of background sequences
#' @param method the algorithm used to detect potential motifs, multiple are allowed
#' @param min.seqs the minimum frequency of a motif should be considered
#' @param ncores the number of cores to be used, passed to \code{mclapply}.
#' @param max.fdr the maximum FDR to be reported
#' @details more details here
#' @return a \code{data.frame} of the over-representation data analysis
#' @importFrom parallel mclapply
#' @importFrom stringr str_detect
#' @export
#'

detectMotifs <- function(fg.seqs, bg.seqs, method = c("cvp", "exh")[1], min.seqs = 5,
                         max.fdr = 1e-2, ncores = 1, verbose = TRUE) {

  method <- match.arg(method, c("cvp", "exh"), several.ok = TRUE)

  if (verbose)
    message("Checking input sequences ...")
  seqs <- checkSeqs(fg.seqs, bg.seqs, option = "extend")

  motif_cvp <- motif_exh <- NULL
  if ("cvp" %in% method) {
    if (verbose)
      message("Detecting potential motifs using 'cvp' algorithm ...")
    motif_cvp <- motif_cvp(seqs$fg.seqs, min.seqs=min.seqs, ncores = ncores)
  }
  if ("exh" %in% method) {
    if (verbose)
      message("Detecting potential motifs using 'exh' algorithm ...")
    motif_exh <- motif_exh(seqs$fg.seqs, min.seqs=min.seqs, ncores = ncores)
  }
  motifs <- c(motif_cvp, motif_exh)
  
  if (length(motifs) == 0) {
    message("No potential motif discovered, try to lower the min.seqs.")
    return (NULL)
  }


  if (verbose)
    message("Evaluating the significance of detected motifs")
  motifor(fg.count = motifs, n.fg.seqs = length(seqs$fg.seqs),
          bg.seqs = seqs$bg.seqs, ncores = ncores, max.fdr = max.fdr)

}
