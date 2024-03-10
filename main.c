#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "reading.h"
#include "simple_opttree.h"
#include "tree.h"

#define DEFAULT_MIN_NODE_SIZE 1
#define DEFAULT_DEPTH 3

void freedata(
   int                   num_cols_x,
   int                   num_cols_y,
   char**                covnames,
   char**                actionnames,
   double*               data_x,
   double*               data_y,
   NODE*                 tree
   )
{
   int i;
   
   for(i = 0; i < num_cols_x; i++)
      free(covnames[i]);
   free(covnames);

   for(i = 0; i < num_cols_y; i++)
      free(actionnames[i]);
   free(actionnames);

   tree_free(tree);
   
   free(data_x);
   free(data_y);
}

int process_commandline(
   int                   argc,               /**< number of command line arguments */
   char*                 argv[],             /**< command line arguments */
   char**                filename,           /**< *filename will be the name of the file with the data */
   int*                  num_cols_y,         /**< *num_cols_y will be number of actions */
   int*                  depth               /**< *depth will be required depth of tree */
   )
{

   if( argc < 2 )
   {
      printf("Need to supply at least a filename and the number of actions.\n");
      return 1;
   }

   *filename = argv[1];
   
   *num_cols_y = atoi(argv[2]);

   if( *num_cols_y == 0 )
   {
      printf("Number of actions either not a valid number of set to 0.\n");
      return 1;
   }
   
   if( argc > 3)
      *depth = atoi(argv[3]);
   else
      *depth = DEFAULT_DEPTH;
   
   if( *depth == 0 )
      printf("Warning: tree depth set to 0.\n");

   return 0;
}

/** 
 * call like ./a.out filename nactions [depth]
 * file should have a single header line and then each subsequent line is covariates followed by nactions reward values.
 * assume white space separator with no white space in names of covariates
 */
int main(
   int                   argc,               /**< number of command line arguments */
   char*                 argv[]              /**< command line arguments */
  )
{

   char* filename;
   int num_cols_y;
   int depth;
   int num_rows;
   int num_cols_x;
   double* data_x;
   double* data_y;
   char** covnames;
   char** actionnames;

   int status;
   NODE* tree;
  
   int min_node_size = DEFAULT_MIN_NODE_SIZE;

   status = process_commandline(argc, argv, &filename, &num_cols_y, &depth); 

   assert( status == 0 || status == 1);
   if( status == 1 )
      return 1;

  status = readfile(filename, num_cols_y, &data_x, &data_y, &num_rows, &num_cols_x, &covnames, &actionnames);

  assert( status == 0 || status == 1);
  if( status == 1 )
     return 1;
     
  tree = tree_search_simple(depth, min_node_size, data_x, data_y, num_rows, num_cols_x, num_cols_y);
     
  print_tree_policytree(tree, covnames, depth, num_cols_y, actionnames);

  printf("Reward: %g\n", get_reward(tree));

  freedata(num_cols_x, num_cols_y, covnames, actionnames, data_x, data_y, tree);

  return 0;
}
