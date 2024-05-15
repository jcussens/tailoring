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
#' @export
fastpolicytree <- function(X, Gamma, depth = 2, min.node.size = 1,
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
    stop("X and Gamma does not have the same number of rows")
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
