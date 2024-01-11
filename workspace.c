#include "sorted_set.h"
#include "tree.h"

#define LEFT 0
#define RIGHT 1

struct workspace
{
   double*               rewards;            /**< size is number of actions */
   SORTED_SET****        sets;               /**< a sorted set for LEFT and RIGHT and each depth and each covariate, 
                                                each with space = number of units */
   NODE**                trees;              /**< a tree for each possible depth */
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
   
   workspace->rewards = (double*) malloc(num_cols_y*sizeof(double));

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

   
   return workspace;
}

void free_workspace(
   WORKSPACE*            workspace,          /**< workspace */
   int                   depth,              /**< (maximum) depth of returned tree */
   int                   num_cols_x          /**< number of covariates */
   )
{
   free(workspace->double_array);

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
   
}

/** return array of doubles, one double for each action */
double* get_rewards(
   WORKSPACE*            workspace           /**< workspace */
   )
{
   return workspace->rewards;
}

SORTED_SET** get_left_sorted_sets(
   WORKSPACE*            workspace,          /**< workspace */
   int                   depth
   )
{
   return workspace->sets[LEFT][depth];
}

SORTED_SET** get_right_sorted_sets(
   WORKSPACE*            workspace,          /**< workspace */
   int                   depth
   )
{
   return workspace->sets[RIGHT][depth];
}

/** record a tree in workspace as the best tree */
void record_best_tree(
   WORKSPACE*            workspace,          /**< workspace */
   NODE*                 tree,
   int                   depth
   )
{
   tree_copy(tree, workspace->trees[depth]);
}

/** retrieve the best tree from workspace */
void retrieve_best_tree(
   WORKSPACE*            workspace,          /**< workspace */
   NODE*                 tree,
   int                   depth
   )
{
   tree_copy(workspace->trees[depth], tree);
}

