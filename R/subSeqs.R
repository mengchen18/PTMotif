#' Find all subsequences of a sequence in give length(s)
#' @author Chen Meng
#' @param x a character vector represents the sequence, each elment should be a single character
#' @param k the length of subsequences
#' @param trim.char the character to be excluded from the sequences
#' @return a character vector contains all the subsequences
#' @examples
#'   subSeqs(x=LETTERS[1:5])
#'   subSeqs(x=sample(LETTERS, 400, replace = TRUE), k = c(2, 4))
#' @export
#'
subSeqs <- function(x, k=2:length(x), trim.char=NULL) {

  if (!is.null(trim.char))
    x <- x[x != "_"]

  n <- length(x)
  k <- k[k<=n]
  st <- rep("", n)
  xl <- lapply(1:n, function(i) {
    c(rep(" ", i-1), x[1:(n-i+1)])
  } )

  xlc <- unlist(lapply(k, function(i) {
    do.call(paste0, xl[i:1])
  }))
  xlc <- trimws(xlc)
  xlc <- xlc[nchar(xlc) %in% k]
  unique(xlc)
}
