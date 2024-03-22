#' @description A package for constructing optimal policy trees from
#'     covariate and reward data. It provides a subset of the
#'     functionality of the policytree package but constructs the
#'     trees more quickly (particularly for discrete data) and thus
#'     allows e.g. deeper trees to be constructed.
#'
#' @useDynLib fastpolicytree
#' @importFrom Rcpp evalCpp
#' @keywords internal
"_PACKAGE"

