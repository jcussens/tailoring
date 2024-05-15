/** @file main.c
 *  @brief main
 *  @author James Cussens
 */
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include "reading.h"
#include "simple_opttree.h"
#include "tree.h"
#include "strategy.h"


#define NOPTARGS 10              /**< number of optional arguments the user can put on the command line */

#define MINNODESIZE_PARAMNAME          "min_node_size"
#define DEPTH_PARAMNAME                "depth"
#define DATATYPE_PARAMNAME             "datatype"
#define VERBOSITY_PARAMNAME            "verbosity"
#define FINDREWARDUB_PARAMNAME         "find_reward_ub"
#define FINDUMMYSPLITREWARD_PARAMNAME  "find_dummy_split_reward"
#define USELASTREWARDS_PARAMNAME       "use_last_rewards"
#define USECUTOFFS_PARAMNAME           "use_cutoffs"
#define USECACHE_PARAMNAME             "use_cache"
#define EXPLOITBINARYVARS_PARAMNAME    "exploit_binaryvars"

#define MINNODESIZE_DEFAULTVAL          1
#define DEPTH_DEFAULTVAL                3
#define DATATYPE_DEFAULTVAL             2
#define VERBOSITY_DEFAULTVAL            0
#define FINDREWARDUB_DEFAULTVAL         0
#define FINDUMMYSPLITREWARD_DEFAULTVAL  0
#define USELASTREWARDS_DEFAULTVAL       1
#define USECUTOFFS_DEFAULTVAL           0
#define USECACHE_DEFAULTVAL             1
#define EXPLOITBINARYVARS_DEFAULTVAL    1

#define MINNODESIZE_MINVAL          1
#define DEPTH_MINVAL                0
#define DATATYPE_MINVAL             0
#define VERBOSITY_MINVAL            0
#define FINDREWARDUB_MINVAL         0
#define FINDUMMYSPLITREWARD_MINVAL  0
#define USELASTREWARDS_MINVAL       0
#define USECUTOFFS_MINVAL           0
#define USECACHE_MINVAL             0
#define EXPLOITBINARYVARS_MINVAL    0

#define MINNODESIZE_MAXVAL          -1
#define DEPTH_MAXVAL                -1
#define DATATYPE_MAXVAL             2
#define VERBOSITY_MAXVAL            2
#define FINDREWARDUB_MAXVAL         1
#define FINDUMMYSPLITREWARD_MAXVAL  1
#define USELASTREWARDS_MAXVAL       1
#define USECUTOFFS_MAXVAL           1
#define USECACHE_MAXVAL             1
#define EXPLOITBINARYVARS_MAXVAL    1


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
   const char**          optnames,           /**< names for optional parameteres */
   int*                  optvals,            /**< values for optional parameters */
   const int*            optmins,            /**< min values for optional parameters */
   const int*            optmaxs             /**< max values for optional parameters */
   )
{
   int i;
   
   if( argc < 3 )
   {
      printf("Need to supply at least a filename and the number of actions.\n");

      printf("Usage: fpt FILENAME NUM_OF_ACTIONS [[OPTNAME OPTVAL] ... ]\nOptional arguments (with default values) are:\n");
      for(i = 0; i < NOPTARGS; i++)
      {
         printf("--%s %d\n", optnames[i], optvals[i]);
      }
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

   for(i = 3; i < argc; i+=2 )
   {
      int j;
      char* argvi;
      int found_arg;

      if(!(argv[i][0] == '-' && argv[i][1] == '-'))
      {
         printf("All optional arguments must start with '--'.\n");
         return 1;
      }

      argvi = argv[i]+2;
      found_arg = 0;
      for(j = 0; j < NOPTARGS; j++)
      {
         if(!strcmp(argvi,optnames[j]))
         {
            int val = atoi(argv[i+1]);
            if( !((optmins[j] == -1 || val >= optmins[j]) && (optmaxs[j] == -1 || val <= optmaxs[j]) ) )
            {
               printf("Value for --%s must be: %d <= value", optnames[j], optmins[j]);
               if( optmaxs[j] == -1 )
                  printf("\n");
               else
                  printf(" <= %d\n", optmaxs[j]);
               return 1;
            }
            optvals[j] = val;
            found_arg = 1;
            break;
         }
      }

      if( !found_arg )
      {
         printf("Unrecognised optional argument: %s\n", argv[i]);
         return 1;
      }
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
   int min_node_size = -1;
   int verbosity = -1;
   int num_rows = -1;
   int num_cols_x = -1;
   double* data_x = NULL;
   double* data_y = NULL;
   char** covnames = NULL;
   char** actionnames = NULL;

   int status;
   NODE* tree = NULL;
   double reward;


   const char* optnames[NOPTARGS] = {
      MINNODESIZE_PARAMNAME, DEPTH_PARAMNAME, DATATYPE_PARAMNAME,
      VERBOSITY_PARAMNAME, FINDREWARDUB_PARAMNAME, FINDUMMYSPLITREWARD_PARAMNAME,
      USELASTREWARDS_PARAMNAME, USECUTOFFS_PARAMNAME, USECACHE_PARAMNAME, EXPLOITBINARYVARS_PARAMNAME };

   /* store default values for all optional arguments */
   int optvals[NOPTARGS] = {
      MINNODESIZE_DEFAULTVAL, DEPTH_DEFAULTVAL, DATATYPE_DEFAULTVAL,
      VERBOSITY_DEFAULTVAL, FINDREWARDUB_DEFAULTVAL, FINDUMMYSPLITREWARD_DEFAULTVAL,
      USELASTREWARDS_DEFAULTVAL, USECUTOFFS_DEFAULTVAL, USECACHE_DEFAULTVAL, EXPLOITBINARYVARS_DEFAULTVAL };

   const int optmins[NOPTARGS] = {
      MINNODESIZE_MINVAL, DEPTH_MINVAL, DATATYPE_MINVAL,
      VERBOSITY_MINVAL, FINDREWARDUB_MINVAL, FINDUMMYSPLITREWARD_MINVAL,
      USELASTREWARDS_MINVAL, USECUTOFFS_MINVAL, USECACHE_MINVAL, EXPLOITBINARYVARS_MINVAL };

   const int optmaxs[NOPTARGS] = {
      MINNODESIZE_MAXVAL, DEPTH_MAXVAL, DATATYPE_MAXVAL,
      VERBOSITY_MAXVAL, FINDREWARDUB_MAXVAL, FINDUMMYSPLITREWARD_MAXVAL,
      USELASTREWARDS_MAXVAL, USECUTOFFS_MAXVAL, USECACHE_MAXVAL, EXPLOITBINARYVARS_MAXVAL };

   STRATEGY* strategy = get_unint_strategy();
   
   status = process_commandline(argc, argv, &filename, &num_cols_y, optnames, optvals, optmins, optmaxs); 

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

   min_node_size = optvals[0];

   depth = optvals[1];
   
   /* either accept the user's choice of dataset representation or decide based on nature of input data */
   if( optvals[2] == 0 )
      use_simple_sets(strategy);
   else if( optvals[2] == 1 )
      use_sorted_sets(strategy);
   else
      decide_datatype(strategy, data_x, num_rows, num_cols_x);

   verbosity = optvals[3];
   
   /* whether to reward upper bounds */
   set_find_reward_ub(strategy, optvals[4]);

   /* whether to compute dummy split rewards */
   set_find_dummy_split_reward(strategy, optvals[5]);

   /* whether to use last rewards (to avoid considering some splits) */
   set_use_last_rewards(strategy, optvals[6]);

   /* whether to use cutoffs (to avoid considering some splits) */
   set_use_cutoffs(strategy, optvals[7]);

   /* whether to use cache (to avoid considering some splits) */
   set_use_cache(strategy, optvals[8]);

   /* whether to exploit binary variables */
   set_exploit_binaryvars(strategy, optvals[9]);
   
   if( num_rows > 0 )
   {
      tree = tree_search_simple(strategy, verbosity, depth, min_node_size, data_x, data_y, num_rows, num_cols_x, num_cols_y, &reward);
      
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
