#' Annotate discovered motif with known motif
#' @author Chen Meng
#' @references http://hprd.org/PhosphoMotif_finder
#' @param x a \code{data.frame} has at least a column names as "motif"
#' @return a \code{data.frame} has two other columns attached at the end 
#'   of origin input \code{x} for known motifs and the corresponding kinases
#' @export

annotKnownMotifs <- function(x) {
  emptvec <- rep(NA, nrow(x))
  s <- mapply(function(m, k) {
    i <- grepl(m, x$motif)
    
    mm <- emptvec
    kk <- emptvec
    
    mm[i] <- m
    kk[i] <- k
    
    list(motif = mm, kinase = kk)
  }, m = knownMotifs$motif, k = knownMotifs$kinase,
  SIMPLIFY = FALSE)
  
  mot <- lapply(s, "[[", "motif")
  kin <- lapply(s, "[[", "kinase")
  
  motvec <- do.call(paste, c(mot, sep = "  "))
  kinvec <- do.call(paste, c(kin, sep = "||"))
  
  x$knownMotifs <- gsub("  NA|NA  |NA", "", motvec)
  x$knownKinases <- gsub("\\|\\|NA|NA\\|\\||NA", "", kinvec)
  x
}

