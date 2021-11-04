#' @description now function returns the time when something happens

now <- function(verbose=TRUE) { 
  if(verbose) {   
    return(as.character(Sys.time()))
  }
}

#' @description msg function returns a message if verbose = TRUE

msg <- function(..., verbose=TRUE) {
  if(verbose) {
    mm <- paste0(...)
    message(mm)
  }
}

#' @description now.msg function returns a message and the time
#'  that something happens when verbose = TRUE

now.msg <- function(...,verbose=TRUE) {
  msg(now(), " - ", ..., verbose = verbose)
} 

#Examples

#now.msg("DeSeq2 done.", verbose = FALSE)
# now.msg ("DeSeq2 done.", verbose = FALSE)
# now(verbose = TRUE)
# now(verbose = FALSE)
# msg("DeSeq2 done.", verbose = TRUE)
# msg("DeSeq2 done.", verbose = FALSE)
