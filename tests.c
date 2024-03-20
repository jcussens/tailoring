/** @file tests.c
 *  @brief Tests
 *  @author James Cussens
 */
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "simple_opttree.h"
#include "tree.h"


/** re-arrange input data */
static
void make_data(
   double                input[],            /**< input */
   int                   num_rows,           /**< number of rows in input */
   int                   num_cols_x,         /**< number of covariates in input */
   int                   num_cols_y,         /**< number of actions in input */
   double                data_x[],           /**< (output) covariate data */
   double                data_y[]            /**< (output) reward data */
   )
{
   int i = 0;

   int row;
   int col;

   for(row = 0; row < num_rows; row++)
   {
      for(col = 0; col < num_cols_x; col++)
         data_x[col*num_rows+row] = input[i++];
      for(col = 0; col < num_cols_y; col++)
         data_y[col*num_rows+row] = input[i++];
   }
}

/** test 1 */
static
void test_1(
   int                   depth               /**< depth of tree */
   )
{
   int num_rows = 1;
   int num_cols_x = 1;
   int num_cols_y = 2;
   
   double input[1*(1+2)] = {1.0,1.0,2.0};
   double data_x[1*1];
   double data_y[1*2];
   char* covnames[1] = {"X1"};
   char* actionnames[2] = {"A1","A2"};

   make_data(input, num_rows, num_cols_x, num_cols_y, data_x, data_y); 
   
   NODE* tree = tree_search_simple(depth, 1, data_x, data_y, num_rows, num_cols_x, num_cols_y);

   print_tree_policytree(tree, covnames, depth, num_cols_y, actionnames);

   tree_free(tree);
}

/** main */
int main(
   int                   argc,               /**< number of command line arguments */
   char*                 argv[]              /**< command line arguments */
  )
{
   int depth;

   for( depth = 0; depth < 5; depth++)
      test_1(depth);
}
