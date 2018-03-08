#' Internal function used to check the input foreground and background sequence vectors
#' @author Chen Meng
#' @param fg.seqs a character vector of foreground sequences
#' @param bg.seqs a character vector of background sequences
#' @param option what should be done when some sequences in fg.seqs do not exist in the bg.seqs.
#'   acceptable options are "extend" - extend the background list; "trim" - trim the foreground list,
#'   i.e. remove the ones in foreground but not in background; "error" - stop and returns
#'   error message.
#' @param center the amino acid centered at the sequences, sequences with other center AA would be 
#'   removed from the list. To disable this function, set \code{center = NULL}.
#' @return a list has two components: fg.seqs and bg.seqs
#' @export
#'
checkSeqs <- function(fg.seqs, bg.seqs, option = c("extend", "trim", "error")[1], center = "STY") {

  # sequence length
  n1 <- unique(nchar(fg.seqs))
  n2 <- unique(nchar(bg.seqs))
  if (length(n1) > 1)
    stop("Foreground sequences need to have the same length.")
  if (length(n2) > 1)
    stop("Background sequences need to have the same length.")
  if (n1 != n2)
    stop("Foreground and background sequneces should have a same length.")

  # duplicated seqs
  if (any(duplicated(fg.seqs))) {
    fg.seqs <- unique(fg.seqs)
    warning("Duplicated sequences detected in foreground list and removed. ")
  }
  if (any(duplicated(bg.seqs))) {
    bg.seqs <- unique(bg.seqs)
    warning("Duplicated sequences detected in background list and removed. ")
  }

  #sequence character, only AA letters and _
  aaa <- c("_","A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
  j <- unique(unlist(strsplit(c(fg.seqs, bg.seqs), "|"))) %in% aaa
  if (!all(j)) 
    stop(paste("fg.seqs and bg.seqs should be amino acid sequences, allowed letters include:", 
      paste(aaa, collapse = ",")))
  
  # if center is not NULL
  if (!is.null(center)) {
    center  <- center[1]
    center <- strsplit(center, "|")[[1]]
    fg.seqs <- fg.seqs[substr(fg.seqs, (n1+1)/2, (n1+1)/2) %in% center]
    bg.seqs <- bg.seqs[substr(bg.seqs, (n2+1)/2, (n2+1)/2) %in% center]
  }
    
  # all fg should be in bg
  options <- match.arg(option, choices = c("extend", "trim", "error"))
  ei <- fg.seqs %in% bg.seqs
  if (!all(ei)) {
    message("Not all foreground sequences exist in background sequences.")
    if (option == "extend") {
      message("Background sequences are extended.")
      bg.seqs <- c(bg.seqs, fg.seqs[!ei])
    } else if (option == "trim") {
      message("foreground sequences are trimmed.")
      fg.seqs <- fg.seqs[ei]
    } else if (option == 'error') {
      stop("Double check your input sequences!")
    } else
      stop("Unknown option!")
  }
  
  

  list(fg.seqs = fg.seqs, bg.seqs = bg.seqs)
}
