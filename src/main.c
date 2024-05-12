/** @file main.c
 *  @brief main
 *  @author James Cussens
 */
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "reading.h"
#include "simple_opttree.h"
#include "tree.h"
#include "strategy.h"

#define DEFAULT_MIN_NODE_SIZE 1 /**< default size for minimum number of units to be in a leaf */
#define DEFAULT_DEPTH 3         /**< default depth limit for a policy tree */

/** free all data */
static
void freedata(
   int                   num_cols_x,         /**< number of covariates */
   int                   num_cols_y,         /**< number of actions */
   char**                covnames,           /**< covnames[j] is the name of covariate j */
   char**                actionnames,        /**< actionnames[j] is the name of action j */
   double*               data_x,             /**< covariates, data_x+(j*num_rows) points to values for covariate j */
   double*               data_y,             /**< data_y[d*num_rows+elt] is the reward for action d for unit elt */
   NODE*                 tree,               /**< policy tree */
   STRATEGY*             strategy            /**< strategy */
   )
{
   int i;

   assert(num_cols_x == 0 || covnames != NULL );
   assert(num_cols_y == 0 || actionnames != NULL );
   
   for(i = 0; i < num_cols_x; i++)
      free(covnames[i]);
   if( covnames != NULL )
      free(covnames);

   for(i = 0; i < num_cols_y; i++)
      free(actionnames[i]);
   if( actionnames != NULL )
      free(actionnames);

   if( tree != NULL )
      tree_free(tree);

   if( data_x != NULL )
      free(data_x);

   if( data_y != NULL )
      free(data_y);
   
   if( strategy != NULL )
      free(strategy);
}

/** process command line arguments
   @return 1 if all is well, otherwise 0 */
static
int process_commandline(
   int                   argc,               /**< number of command line arguments */
   char**                argv,               /**< command line arguments */
   char**                filename,           /**< *filename will be the name of the file with the data */
   int*                  num_cols_y,         /**< *num_cols_y will be number of actions */
   int*                  depth,              /**< *depth will be required depth of tree */
   int*                  optargs,            /**< for optional arguments supplied on command line */
   int                   noptargs            /**< number of values for optargs */
   )
{
   int i;
   
   if( argc < 3 )
   {
      printf("Need to supply at least a filename and the number of actions.\n");

      printf("Usage: fpt FILENAME NUM_OF_ACTIONS [DEPTH] [DATA_REP] [COMPUTE_UBS] [COMPUTE_DUMMY] [USE_LAST_REWARDS] [USE_CUTOFFS] [USE_CACHE] [EXPLOIT_BINARY]\n");
      printf("Default value of DEPTH is 3, for DATA_REP, 0 means delayed sorting, 1 means policytree style, any other value (the default) means use data to decide.\n");
      printf("All other optional arguments are 0 (off) or 1 (on), with 1 being the default.\n");

      return 1;
   }

   *filename = argv[1];
   
   *num_cols_y = atoi(argv[2]);

   if( *num_cols_y == 0 )
   {
      printf("Number of actions either not a valid number of set to 0.\n");
      return 1;
   }
   else if( *num_cols_y == 1 )
   {
      printf("Warning: number of actions set to 1.\n");
   }
   
   if( argc > 3)
      *depth = atoi(argv[3]);
   else
      *depth = DEFAULT_DEPTH;
   
   if( *depth == 0 )
      printf("Warning: tree depth set to 0.\n");

   
   for( i = 4; i < 4 + noptargs; i++ )
   {
      int* argptr = optargs + i - 4;
      if( argc > i)
      {
         *argptr = atoi(argv[i]);
         if(*argptr != 0 && *argptr != 1)
            printf("Warning: Argument %d ignored.\n", i);
      }
      else
         *argptr = -1;
   }
   
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

   char* filename = NULL;
   int num_cols_y = -1;
   int depth = -1;
   int num_rows = -1;
   int num_cols_x = -1;
   double* data_x = NULL;
   double* data_y = NULL;
   char** covnames = NULL;
   char** actionnames = NULL;

   int status;
   NODE* tree = NULL;
   double reward;
   
   int min_node_size = DEFAULT_MIN_NODE_SIZE;

   int optargs[6];
   STRATEGY* strategy = get_unint_strategy();
   
   /* NODE** nodes; */
   /* int num_nodes; */
   /* int i; */

   status = process_commandline(argc, argv, &filename, &num_cols_y, &depth, optargs, 6); 

   assert( status == 0 || status == 1);
   if( status == 1 )
      return 1;

   assert(filename != NULL);
   assert(num_cols_y >= 1);
   assert(depth >= 0);
   assert(min_node_size >= 1);
   
   status = readfile(filename, num_cols_y, &data_x, &data_y, &num_rows, &num_cols_x, &covnames, &actionnames);
   
   assert( status == 0 || status == 1);
   if( status == 1 )
      return 1;

   assert( num_rows >= 0 );
   assert( num_rows == 0 || data_x != NULL );
   assert( num_rows == 0 || data_y != NULL );
   assert( num_cols_x >= 0 );
   assert( actionnames != NULL );
   assert( num_cols_x == 0 || covnames != NULL );
   
   if( num_cols_x == 0 )
   {
      printf("Warning: No covariates supplied.\n");     
   }

   if( num_rows == 0 )
   {
      printf("Warning: No data supplied.\n");     
   }

   /* either accept the user's choice of dataset representation or decide based on nature of input data */
   if( optargs[0] == 0 )
      use_simple_sets(strategy);
   else if( optargs[0] == 1 )
      use_sorted_sets(strategy);
   else
      decide_datatype(strategy, data_x, num_rows, num_cols_x);

   /* default is to compute reward upper bounds */
   set_find_reward_ub(strategy, (optargs[1] == 0) ? 0 : 1);

   /* default is to compute dummy split rewards */
   set_find_dummy_split_reward(strategy, (optargs[2] == 0) ? 0 : 1);

   /* default is to use last rewards (to avoid considering some splits) */
   set_use_last_rewards(strategy, (optargs[3] == 0) ? 0 : 1);

   /* default is to use cutoffs (to avoid considering some splits) */
   set_use_cutoffs(strategy, (optargs[4] == 0) ? 0 : 1);

   /* default is to use cache (to avoid considering some splits) */
   set_use_cache(strategy, (optargs[5] == 0) ? 0 : 1);

   /* set exploit binary variables */
   set_exploit_binaryvars(strategy, 1);
   
   if( num_rows > 0 )
   {
      tree = tree_search_simple(strategy, depth, min_node_size, data_x, data_y, num_rows, num_cols_x, num_cols_y, &reward);
      
      print_tree_policytree(tree, covnames, depth, num_cols_y, actionnames);
      
      printf("Reward: %g\n", reward);
   }
   else
   {
      printf("Did not build tree since no data supplied.\n");
   }

   printf("Git commit is " GITHASH ".\n");
      
   /* nodes = breadth_first_nodes(tree, depth, &num_nodes); */
   /* for( i = 0; i < num_nodes; i++ ) */
   /* { */
   /*    printf("%d/%d::\n",i,num_nodes); */
   /*    print_tree(nodes[i],NULL); */
   /* } */
   
   freedata(num_cols_x, num_cols_y, covnames, actionnames, data_x, data_y, tree, strategy);
   
   return 0;
}
