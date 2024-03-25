#include <Rcpp.h>
using namespace Rcpp;

#include "simple_opttree.h"
#include "tree.h"

// This code uses/adapts C++ source contained in the policytree R package 

// [[Rcpp::export]]
Rcpp::List tree_search_rcpp(
   const Rcpp::NumericMatrix& X,
   const Rcpp::NumericMatrix& Y,
   int                   depth,
   int                   min_node_size
   )
{
   int num_rows = X.rows();
   int num_cols_x = X.cols();
   int num_cols_y = Y.cols();
   double reward;
   int perfect;

   NODE* root = tree_search_simple(depth, min_node_size, X.begin(), Y.begin(),
				   num_rows, num_cols_x, num_cols_y, &reward,
				   &perfect);
   int num_nodes;
   NODE** treenodes = breadth_first_nodes(root, depth, &num_nodes);
   
   Rcpp::NumericMatrix tree_array(num_nodes, 4);
   Rcpp::List nodes;

   int i = 1;
   int j = 0;
   int idx;
   
   // construct node list and array representation 
   for( idx = 0; idx < num_nodes; idx++)
   {
      NODE* node = treenodes[idx];

      // just ignore non-nodes
      if( node == NULL )
	continue;
      
      if( is_leaf(node) )
      {
         int action = get_action(node);
         auto list_node = Rcpp::List::create(Rcpp::Named("is_leaf") = true,
            Rcpp::Named("action") = action + 1); // C++ index
         nodes.push_back(list_node);
         tree_array(j, 0) = -1;
         tree_array(j, 1) = action + 1;
      }
      else
      {
         int index = get_index(node);
         double value = get_value(node);
         auto list_node = Rcpp::List::create(Rcpp::Named("is_leaf") = false,
            Rcpp::Named("split_variable") = index + 1, // C++ index
            Rcpp::Named("split_value") = value,
            Rcpp::Named("left_child") = i + 1,
            Rcpp::Named("right_child") = i + 2);
         nodes.push_back(list_node);
         tree_array(j, 0) = index + 1;
         tree_array(j, 1) = value;
         tree_array(j, 2) = i + 1; // left child
         tree_array(j, 3) = i + 2; // right child
         i += 2;
      }
      j++;
   }
   free(treenodes);
   
   Rcpp::List result;
   result.push_back(nodes);
   result.push_back(tree_array);

   return result;
}
   
   

