#include <Rcpp.h>
using namespace Rcpp;

#include "simple_opttree.h"
#include "tree.h"

// [[Rcpp::export]]
Rcpp::List tree_search_rcpp(
   const Rcpp::NumericMatrix& X,             /**< covariates */
   const Rcpp::NumericMatrix& Y,             /**< rewards */
   int                   depth,              /**< depth of tree */
   int                   min_node_size       /**< minimal node size */
   )
{
   int num_rows = X.rows();
   int num_cols_x = X.cols();
   int num_cols_y = Y.cols();
   /* '.begin()' returns an iterator pointing to the first element */
   /* R matrices are column-major */
   NODE* root = tree_search_simple(depth, min_node_size, X.begin(), Y.begin(),
      num_rows, num_cols_x, num_cols_y);
   Rcpp::NumericMatrix tree_array(num_nodes, 4);
   Rcpp::List nodes;

   /* construct nodes and array representation */
   

   

