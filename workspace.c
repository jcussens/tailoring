#include "workspace.h"
#include <stdlib.h>

#define LEFT 0
#define RIGHT 1

struct workspace
{
   int                   num_cols_y;         /**< number of rewards */
   double*               rewards;            /**< size is number of actions */
   double*               rewards2;           /**< size is number of actions */
   SORTED_SET****        sets;               /**< a sorted set for LEFT and RIGHT and each depth and each covariate, 
                                                each with space = number of units */
   NODE**                trees;              /**< a tree for each possible depth */
   SORTED_SET*           sorted_set;         /**< an uninitialised sorted set */
};

/** make workspace to provide pre-allocated space for various functions */
WORKSPACE* make_workspace(
   int                   depth,              /**< (maximum) depth of returned tree */
   SORTED_SET**          initial_sorted_sets, /**< initial sorted sets */
   int                   num_rows,           /**< number of units in full dataset */
   int                   num_cols_x,         /**< number of covariates */
   int                   num_cols_y          /**< number of rewards/actions */
   )
{
   WORKSPACE* workspace;
   int d;
   int i;

   workspace = (WORKSPACE*) malloc(sizeof(WORKSPACE));

   workspace->num_cols_y = num_cols_y;
   
   workspace->rewards = (double*) malloc(num_cols_y*sizeof(double));
   workspace->rewards2 = (double*) malloc(num_cols_y*sizeof(double));

   workspace->sets = (SORTED_SET****) malloc(2*sizeof(SORTED_SET***));
   for( i = 0; i < 2; i++)
   {
      workspace->sets[i] = (SORTED_SET***) malloc(depth*sizeof(SORTED_SET**));
      for( d = 0; d < depth; d++)
      {
         workspace->sets[i][d] = shallow_copy_sorted_sets(
            initial_sorted_sets, num_cols_x);
      }
   }

   workspace->trees = (NODE**) malloc(depth*sizeof(NODE*));
   for(d = 0; d < depth; d++)
   {
      workspace->trees[d] = make_tree(d);
   }

   workspace->sorted_set = make_uninitialised_sorted_set();
   
   return workspace;
}

/** free the workspace */
void free_workspace(
   WORKSPACE*            workspace,          /**< workspace */
   int                   depth,              /**< (maximum) depth of returned tree */
   int                   num_cols_x          /**< number of covariates */
   )
{
   int i;
   int d;
   
   free(workspace->rewards);
   free(workspace->rewards2);

   for( i = 0; i < 2; i++)
   {
      for( d = 0; d < depth; d++)
      {
         shallow_free_sorted_sets(workspace->sets[i][d], num_cols_x);
      }
      free(workspace->sets[i]);
   }
   free(workspace->sets);

   for( d = 0; d < depth; d++)
      tree_free(workspace->trees[d]);
   free(workspace->trees);

   free(workspace->sorted_set);
   
}

/** return array of doubles, one double for each action */
double* get_rewards(
   WORKSPACE*            workspace           /**< workspace */
   )
{
   return workspace->rewards;
}

/** return array of zeroes, one zero for each reward */
double* get_rewards2(
   WORKSPACE*            workspace           /**< workspace */
   )
{
   int d;
   for(d = 0; d < workspace->num_cols_y; d++)
      workspace->rewards2[d] = 0.0;
   return workspace->rewards2;
}


/** get left sorted sets associated with a given depth */
SORTED_SET** get_left_sorted_sets(
   WORKSPACE*            workspace,          /**< workspace */
   int                   depth               /**< depth */
   )
{
   return workspace->sets[LEFT][depth];
}

/** get right sorted sets associated with a given depth */
SORTED_SET** get_right_sorted_sets(
   WORKSPACE*            workspace,          /**< workspace */
   int                   depth               /**< depth */
   )
{
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
   WORKSPACE*            workspace,          /**< workspace */
   NODE*                 tree,               /**< tree */
   int                   depth               /**< depth of tree */
   )
{
   tree_copy(workspace->trees[depth], tree);
}

/** return an uninitialised sorted set */
SORTED_SET* get_sorted_set(
   WORKSPACE*            workspace           /**< workspace */
   )
{
   return workspace->sorted_set;
}

