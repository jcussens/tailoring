#' Construct an optimal policy tree from covariate and reward data
#'
#' This function accepts almost the same input and generates the same
#' type of output as the policy_tree function in the policytree
#' package. The only difference is that this function has no
#' 'split.step' argument (since it is effectively hard-coded to the
#' value 1).
#' 
#' @param X The covariates used. Dimension \eqn{N*p} where \eqn{p} is the number of features.
#' @param Gamma The rewards for each action. Dimension \eqn{N*d} where \eqn{d} is the number of actions.
#' @param depth The depth of the fitted tree. Default is 3.
#' @param min.node.size An integer indicating the smallest terminal node size permitted. Default is 1.
#' @param strategy.datatype If set to 0 policytree style sorted sets are used to represent datasets during solving. If set to 1 then unsorted sets are used which are sorted 'on demand'. If set to to 2 then the choice of representation is decided automatically. Default is 2 (choice is automatically made).
#' @param strategy.find.reward.ub If TRUE upper bounds on rewards are computed. Default is FALSE
#' @param strategy.find.dummy.split.reward If TRUE then the reward for a dummy split (where the left split has no datapoints) is computed. Default is FALSE.
#' @param strategy.use.last.rewards If TRUE an upper bound on the reward for a split is computed from the reward for the most recent split value for the current covariate. Default is TRUE
#' @param strategy.use.cutoffs If TRUE then tree finding is aborted if it can be deduced that the reward for the tree is beaten by some existing tree. Default is FALSE
#' @param strategy.use.cache If TRUE a cache of optimal trees for (sub-)datasets is used. Default is TRUE
#' @param strategy.exploitbinaryvars If TRUE then covariates with only 2 values are treated specially. Default is TRUE
#' 
#' @return A policy_tree object.
#'
#' @examples
#' X <-  data.frame(
#'  X1=c(-0.32, 0.16, 0.34, 1.24, 0.22, 0.45, 1.48, 0.65,-0.93,-1.11),
#'  X2=c(-0.58, 0.90,-0.22, 1.54,-0.57,-1.08,-1.42,-1.98,-0.02, 0.05),
#'  X3=c(0.70,-1.49, 0.36,-0.05,-0.14, 1.57,-0.18,-1.98,-1.77,-1.25),
#'  X4=c(0.21, 0.34, 0.60,-0.05,-0.66,-0.69, 0.52, 0.31,-0.03, 1.09),
#'  X5=c(0.16, 0.96,-1.07,-0.97, 2.02,-0.43,-0.79,-2.08, 1.21, 0.39))
#' Gamma  <- data.frame(
#'  control=c(0.8502363,-1.4950411,1.9608062,0.7487925,2.9718517,
#'   0.8952429,-0.2563680,5.9945581,-1.8485703,-1.2840477),
#'  treat=c(-2.91607259,-2.25464535, 0.28214637,-0.17284650,-0.09480810,
#'   1.48786125,2.08600119,-2.05283394,0.72903608,-0.04416392))
#' tree3  <- fastpolicytree(X,Gamma)
#' tree3
#' tree2  <- fastpolicytree(X,Gamma,depth=2)
#' tree2
#' # to get a human-readable display of the trees use the
#' # policytree package...
#' #library(policytree)
#' #print(tree3)
#' #print(tree2)
#' @export

## This R code was produced by copying and editing the file policy_tree.R from
## the policytree package

fastpolicytree <- function(X, Gamma, depth = 3, min.node.size = 1,
                           strategy.datatype  = 2,
                           strategy.find.reward.ub = FALSE,
                           strategy.find.dummy.split.reward = FALSE,
                           strategy.use.last.rewards = TRUE,
                           strategy.use.cutoffs = FALSE,
                           strategy.use.cache = TRUE,
                           strategy.exploitbinaryvars = TRUE
                           ) {
  n.features <- ncol(X)
  n.actions <- ncol(Gamma)
  n.obs <- nrow(X)
  valid.classes <- c("matrix", "data.frame")

  
  if (!inherits(X, valid.classes) || !inherits(Gamma, valid.classes)) {
    stop(paste("Currently the only supported data input types are:",
               "`matrix`, `data.frame`"))
  }
  if (!is.numeric(as.matrix(X)) || any(dim(X) == 0)) {
    stop("The feature matrix X must be numeric")
  }
  if (!is.numeric(as.matrix(Gamma)) || any(dim(Gamma) == 0)) {
    stop("The reward matrix Gamma must be numeric")
  }
  if (anyNA(X)) {
    stop("Covariate matrix X contains missing values.")
  }
  if (anyNA(Gamma)) {
    stop("Gamma matrix contains missing values.")
  }
  if (depth < 0 ) {
    stop("`depth` cannot be negative.")
  }
  if (n.obs != nrow(Gamma)) {
    stop("X and Gamma do not have the same number of rows")
  }
  if (as.integer(min.node.size) != min.node.size || min.node.size < 1) {
    stop("min.node.size should be an integer greater than or equal to 1.")
  }
  if( as.integer(strategy.datatype) != strategy.datatype || strategy.datatype < 0 || strategy.datatype > 2 )
      stop("strategy.datatype should be either 0, 1 or 2.")

  
  action.names <- colnames(Gamma)
  if (is.null(action.names)) {
    action.names <- as.character(1:ncol(Gamma))
  }
  action.names <- utils::type.convert(action.names, as.is = TRUE) # TRUE to not convert character to factor
  columns <- colnames(X)
  if (is.null(columns)) {
    columns <- make.names(1:ncol(X))
  }

  result <- tree_search_rcpp(as.matrix(X), as.matrix(Gamma), depth, min.node.size,
                             strategy.datatype,
                           as.integer(strategy.find.reward.ub),
                           as.integer(strategy.find.dummy.split.reward),
                           as.integer(strategy.use.last.rewards),
                           as.integer(strategy.use.cutoffs),
                           as.integer(strategy.use.cache),
                           as.integer(strategy.exploitbinaryvars))
  
  tree <- list(nodes = result[[1]])

  tree[["_tree_array"]] <- result[[2]]
  tree[["depth"]] <- depth
  tree[["n.actions"]] <- n.actions
  tree[["n.features"]] <- n.features
  tree[["action.names"]] <- action.names
  tree[["columns"]] <- columns
  class(tree) <- "policy_tree"

  tree
}
