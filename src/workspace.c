/** @file workspace.c
 *  @brief Workspace functions
 *  @author James Cussens
 */

#include "workspace.h"
#include "units.h"
#include "tree.h"
#include <assert.h>
#include <stdlib.h>

#define LEFT 0                               /**< index for left children */
#define RIGHT 1                              /**< index for right children */

/** pre-allocated workspace used by many functions */
struct workspace
{
   int                   num_cols_y;         /**< number of rewards */
   double*               rewards;            /**< size is number of actions */
   double*               rewards2;           /**< size is number of actions */
   UNITS**               sets;               /**< units for LEFT and RIGHT and each depth and each covariate, 
                                                each with space = number of units */
   UNITS                 tmpunits;           /**< additional space = number of units */
   int*                  tmp2;               /**< additional space = number of units */
   NODE**                trees;              /**< a tree for each possible depth */
};

/** make workspace to provide pre-allocated space for various functions */
WORKSPACE* make_workspace(
   const STRATEGY*       strategy,           /**< tree-building strategy */
   int                   depth,              /**< (maximum) depth of returned tree */
   CONST_UNITS           initial_units,      /**< initial units */
   int                   num_rows,           /**< number of units in full dataset */
   int                   num_cols_x,         /**< number of covariates */
   int                   num_cols_y          /**< number of rewards/actions */
   )
{
   WORKSPACE* workspace;
   int d;
   int i;
   /* for tree of, say, depth=3, we have tree depths of 0,1,2 and 3 to consider */
   int ndepths = depth+1;
   
   workspace = (WORKSPACE*) malloc(sizeof(WORKSPACE));

   workspace->num_cols_y = num_cols_y;
   
   workspace->rewards = (double*) malloc(num_cols_y*sizeof(double));
   workspace->rewards2 = (double*) malloc(num_cols_y*sizeof(double));

   workspace->sets = (UNITS**) malloc(2*sizeof(UNITS*));
   for( i = 0; i < 2; i++)
   {
      workspace->sets[i] = (UNITS*) malloc(ndepths*sizeof(UNITS));
      for( d = 0; d < ndepths; d++)
      {
         workspace->sets[i][d] = shallow_copy_units(strategy, initial_units, num_cols_x);
      }
   }

   workspace->tmpunits = shallow_copy_units(strategy, initial_units, num_cols_x);
   workspace->tmp2 = (int*) malloc(num_rows*sizeof(int));

   workspace->trees = (NODE**) malloc(ndepths*sizeof(NODE*));
   for(d = 0; d < ndepths; d++)
   {
      workspace->trees[d] = make_tree(d);
   }

   return workspace;
}

/** free the workspace */
void free_workspace(
   const STRATEGY*       strategy,           /**< tree-building strategy */
   WORKSPACE*            workspace,          /**< workspace */
   int                   depth,              /**< (maximum) depth of returned tree */
   int                   num_cols_x          /**< number of covariates */
   )
{
   int i;
   int d;
   /* for tree of, say, depth=3, we have tree depths of 0,1,2 and 3 to consider */
   int ndepths = depth+1;
   
   free(workspace->rewards);
   free(workspace->rewards2);

   for( i = 0; i < 2; i++)
   {
      for( d = 0; d < ndepths; d++)
      {
         shallow_free_units(strategy, workspace->sets[i][d], num_cols_x);
      }
      free(workspace->sets[i]);
   }
   free(workspace->sets);

   shallow_free_units(strategy, workspace->tmpunits, num_cols_x);
   free(workspace->tmp2);
   
   for( d = 0; d < ndepths; d++)
      tree_free(workspace->trees[d]);
   free(workspace->trees);

   free(workspace);
   
}

/** return array of doubles, one double for each action */
double* get_rewards_space(
   const WORKSPACE*      workspace           /**< workspace */
   )
{
   return workspace->rewards;
}

/** return tmp units */
UNITS get_tmpunits(
   const WORKSPACE*      workspace           /**< workspace */
   )
{
   return workspace->tmpunits;
}

/** return secondary tmp units */
int* get_tmp2(
   const WORKSPACE*      workspace           /**< workspace */
   )
{
   return workspace->tmp2;
}

   
/** return array of zeroes, one zero for each reward */
double* get_rewards_space_zeroed(
   WORKSPACE*            workspace           /**< workspace */
   )
{
   int d;
   for(d = 0; d < workspace->num_cols_y; d++)
      workspace->rewards2[d] = 0.0;
   return workspace->rewards2;
}


/** get left sorted sets associated with a given depth */
UNITS get_left_sorted_sets(
   const WORKSPACE*      workspace,          /**< workspace */
   int                   depth               /**< depth */
   )
{
   assert( workspace != NULL );
   assert( depth >= 0 );
   assert( workspace->sets != NULL );
   assert( workspace->sets[LEFT] != NULL );
   assert( workspace->sets[LEFT][depth] != NULL );
   
   return workspace->sets[LEFT][depth];
}

/** get right sorted sets associated with a given depth */
UNITS get_right_sorted_sets(
   const WORKSPACE*      workspace,          /**< workspace */
   int                   depth               /**< depth */
   )
{
   assert( workspace != NULL );
   assert( depth >= 0 );

   return workspace->sets[RIGHT][depth];
}

/** record a tree of given depth in workspace as the best tree */
void record_best_tree(
   WORKSPACE*            workspace,          /**< workspace */
   NODE*                 tree,               /**< tree */
   int                   depth               /**< depth of tree */
   )
{
   tree_copy(tree, workspace->trees[depth]);
}

/** retrieve the best tree of given depth from workspace */
void retrieve_best_tree(
   const WORKSPACE*      workspace,          /**< workspace */
   NODE*                 tree,               /**< tree */
   int                   depth               /**< depth of tree */
   )
{
   tree_copy(workspace->trees[depth], tree);
}


