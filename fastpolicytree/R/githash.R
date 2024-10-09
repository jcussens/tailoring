#' Returns git hash for compiled C code
#'
#' @return A string which is the relevant githash
#'
#' @export
githash <- function() {
    githash_rcpp()
}
