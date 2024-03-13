#include "simple_opttree.h"
#include "sorted_set.h"
#include "workspace.h"
#include <assert.h>
#include <stdlib.h>

/* #define VERYVERBOSE */

#ifdef VERYVERBOSE
#define VERBOSE
#endif

#ifdef VERBOSE
#include <stdio.h>
#endif


static
void update_left_rewards(
   double*               left_rewards,
   int*                  elts,
   int                   nelts,
   const double*         data_y,             /**< rewards (column major) */
   int                   num_rows,           /**< number of units */
   int                   num_cols_y          /**< number of rewards/actions */
   )
{
   int i;
   
   for( i = 0; i < nelts; i++ )
   {
      int elt = elts[i];
      int d;
      
      if( num_cols_y == 2 )
      {
        left_rewards[0] += data_y[elt];
        left_rewards[1] += data_y[num_rows+elt];
      }
      else
      {
        for( d = 0; d < num_cols_y; d++ )
          left_rewards[d] += data_y[d*num_rows+elt];
      }
   }
}

static
void update_best_right_reward
(
   double*               left_rewards,
   double*               nosplit_rewards,
   int                   num_cols_y,         /**< number of rewards/actions */
   double*               best_right_reward,
   int*                  best_right_action
   )
{
   int d;
   
   *best_right_reward = nosplit_rewards[0] - left_rewards[0];
   *best_right_action = 0;
      
   for( d = 1; d < num_cols_y; d++ )
   {
      double right_reward = nosplit_rewards[d] - left_rewards[d];
      if( right_reward > *best_right_reward )
      {
         *best_right_reward = right_reward;
         *best_right_action = d;
      }
   }   
}


static
void update_best_left_reward(
   double*               left_rewards,
   int                   num_cols_y,         /**< number of rewards/actions */
   double*               best_left_reward,
   int*                  best_left_action
   )
{
   int d;

   *best_left_reward = left_rewards[0];
   *best_left_action = 0;
      
   for( d = 1; d < num_cols_y; d++ )
   {
      if( left_rewards[d] > *best_left_reward )
      {
         *best_left_reward = left_rewards[d];
         *best_left_action = d;
      }
   }   
}

/**
 * For each unit find and record the best and worst actions for that unit
 */
static
void store_best_worst_actions(
   const double*         data_y,             /**< rewards (column major) */
   int                   num_rows,           /**< number of units */
   int                   num_cols_y,         /**< number of rewards */
   int**                 best_actions,       /**< (pointer to) best actions */
   int**                 worst_actions       /**< (pointer to) worst actions */
  )
{
   int i;
   int best_action;
   int worst_action;
   int d;
   double best_reward;
   double worst_reward;

   assert(num_cols_y > 1);

   *best_actions = (int*) malloc(num_rows*sizeof(int));
   *worst_actions = (int*) malloc(num_rows*sizeof(int));  

   for( i = 0; i < num_rows; i++ )
   {
      /* there are always at least two actions
         initialise best to action 0
         and worst to action 1, so that even
         if we get equal rewards the best and worst action must differ
         (it is just simpler to assume that best and worst differ
      */
      
      best_action = 0;
      best_reward = data_y[i];
      worst_action = 1;
      worst_reward = data_y[num_rows+i];
    
      for( d = 0; d < num_cols_y; d++ )
      {
         if( data_y[d*num_rows+i] > best_reward )
         {
            best_action = d;
            best_reward = data_y[d*num_rows+i];
         }
         if( data_y[d*num_rows+i] < worst_reward )
         {
            worst_action = d;
            worst_reward = data_y[d*num_rows+i];
         }
      }
      (*best_actions)[i] = best_action;
      (*worst_actions)[i] = worst_action;
      assert(worst_action != best_action);
      assert(worst_reward <= best_reward);
   }
}


/** Find an optimal depth=1 (ie at most one single split) tree for a set of units 
 */
static
void level_one_learning(
   NODE*                 node,               /**< uninitialised tree to be populated with optimal tree */
   const SORTED_SET**    sorted_sets,        /**< sorted sets for the units, one for each covariate */
   const double*         data_x,             /**< covariates, data_x+(j*num_rows) points to values for covariate j */
   const double*         data_y,             /**< gammas, data_y+(d*num_rows) points to values for reward d */
   int                   num_rows,           /**< number of units in full dataset */
   int                   num_cols_x,         /**< number of covariates */
   int                   num_cols_y,         /**< number of rewards/actions */
   const int*            best_actions,       /**< best_actions[i] is the best action for unit i */
   const int*            worst_actions,      /**< worst_actions[i] is the worst action for unit i */
   WORKSPACE*            workspace,          /**< workspace */
   int*                  perfect             /**< *perfect=1 iff each unit assigned its best action, else *perfect=0 */
   )
{
   int p;
   double best_reward;
   int first_reward = 1;
   /* next 2 lines just get pointers to pre-allocated space */
   double* nosplit_rewards = get_rewards_space(workspace);
   SORTED_SET* right_sorted_set = get_unint_sorted_set(workspace);
   NODE* left_child;
   NODE* right_child;
   
   int* elts;
   int nelts;

   double best_left_reward;
   double best_right_reward;
   int best_left_action;
   int best_right_action;

#ifdef VERYVERBOSE
   double bestsplitval;
   int bestp;
#endif

   
   /* get reward for each action if no split were done */
   find_nosplit_rewards((const SORTED_SET**) sorted_sets, num_cols_y, data_y, num_rows, nosplit_rewards);

   get_children(node, &left_child, &right_child);
   
   /* consider each covariate for splitting */
   for( p = 0; !(*perfect) && p < num_cols_x; p++)
   {
      const double* data_xp = data_x+(p*num_rows);
      double splitval;
      
      /* initialise all left rewards to 0 */
      double* left_rewards = get_rewards_space_zeroed(workspace);

      /* make very shallow copy of sorted set for covariate p */
      very_shallow_copy(sorted_sets[p], right_sorted_set);
      
      /* consider each split (x[p] <= splitval, x[p] > splitval) of the data */
      while( next_shallow_split(right_sorted_set, data_xp, &splitval, &elts, &nelts) )
      {
         double this_reward;
         
         update_left_rewards(left_rewards, elts, nelts, data_y, num_rows, num_cols_y);
         update_best_left_reward(left_rewards, num_cols_y, &best_left_reward, &best_left_action);
         update_best_right_reward(left_rewards, nosplit_rewards, num_cols_y, &best_right_reward, &best_right_action);

         /* get reward for this split */
         this_reward = best_left_reward + best_right_reward;

         /* if best so far, update */
         if( first_reward || this_reward > best_reward ) 
         {
            best_reward = this_reward;
            record_split(node, p, splitval, best_reward);
            make_leaf(left_child, best_left_reward, best_left_action);
            make_leaf(right_child, best_right_reward, best_right_action);
#ifdef VERYVERBOSE
            bestsplitval = splitval;
            bestp = p;
#endif
            first_reward = 0;
         }
      }
   }
#ifdef VERYVERBOSE
   printf("Best split for depth=1 tree is split value %g for covariate %d with reward %g.\n",
      bestsplitval, bestp, best_reward);
#endif
}

/** on return, `node` will be the root of an optimal tree of depth `depth` for the data
 * represented by `sorted_sets`
 */
static
void find_best_split(
   NODE*                 node,               /**< uninitialised tree to be populated with optimal tree */
   int                   depth,              /**< depth of tree */
   const SORTED_SET**    sorted_sets,        /**< sorted sets for the units, one for each covariate */
   int                   min_node_size,      /**< smallest terminal node size */
   const double*         data_x,             /**< covariates, data_x+(j*num_rows) points to values for covariate j */
   const double*         data_y,             /**< gammas, data_y+(d*num_rows) points to values for reward d */
   int                   num_rows,           /**< number of units in full dataset */
   int                   num_cols_x,         /**< number of covariates */
   int                   num_cols_y,         /**< number of rewards/actions */
   const int*            best_actions,       /**< best_actions[i] is the best action for unit i */
   const int*            worst_actions,      /**< worst_actions[i] is the worst action for unit i */
   WORKSPACE*            workspace,          /**< working space */
   int*                  perfect             /**< *perfect=1 iff each unit assigned its best action, else *perfect=0 */
  )
{

   int p;
   int pure;
   double best_reward;
   int first_reward = 1;
   NODE* left_child;
   NODE* right_child;

   SORTED_SET** left_sorted_sets;
   SORTED_SET** right_sorted_sets;

   assert( node != NULL );
   assert( depth >= 0 );
   assert( sorted_sets != NULL);
   assert( sorted_sets[0] != NULL);
   assert( min_node_size >= 1 );
   assert( num_rows >= 1 );
   assert( num_cols_x >= 0 );
   assert( num_cols_y >= 1 );
   assert( num_cols_x == 0 || data_x != NULL );
   assert( data_y != NULL );
   assert( are_sorted_sets(sorted_sets, data_x, num_rows, num_cols_x) );

#ifdef VERBOSE
   printf("Looking for an optimal depth=%d tree for a dataset of size %d.\n", depth, get_size(sorted_sets));
#endif
   
   /* determine whether the dataset is pure, i.e. whether each unit has same best action */
   pure = is_pure(sorted_sets, best_actions);

   /* if this dataset is pure we can make a perfect leaf */
   if( pure )
      *perfect = 1;

   /* nothing further to do if depth limit reached or if too few datapoints for splitting or if dataset is pure */
   if( depth == 0 || get_size(sorted_sets) <= min_node_size || pure )
   {
      int best_action;

      /* find best action and its associated reward */
      find_best_reward(sorted_sets, data_y, num_rows, num_cols_y, workspace, &best_reward, &best_action);

      /* make node a leaf with found best action and associated reward */
      make_leaf(node, best_reward, best_action);

      /* all done */
      return;
   }

   /* use specialised routine if only one split allowed */
   if( depth == 1 )
   {
      level_one_learning(node, sorted_sets, data_x, data_y,
         num_rows, num_cols_x, num_cols_y, best_actions, worst_actions, workspace, perfect); 
      
      return;
   }

   /* get left and right subtrees */
   get_children(node, &left_child, &right_child);

   assert( left_child != NULL );
   assert( right_child != NULL );

   /* consider each covariate for splitting (stopping if we find a perfect split) */
   for( p = 0; !(*perfect) && p < num_cols_x; p++)
   {
      const double* data_xp = data_x+(p*num_rows);

      double this_reward;

      int left_perfect = 0;
      int right_perfect = 0;

      double splitval;

      /* initialise so that each left_sorted_set (for each covariate) is empty and each 
         right_sorted_set is a copy of the input sorted set (for that covariate) */
      initialise_sorted_sets(sorted_sets, depth, num_cols_x, workspace, &left_sorted_sets, &right_sorted_sets);

      assert( are_sorted_sets((const SORTED_SET**) left_sorted_sets, data_x, num_rows, num_cols_x) );
      assert( are_sorted_sets((const SORTED_SET**) right_sorted_sets, data_x, num_rows, num_cols_x) );

      /* consider each split (x[p] <= splitval, x[p] > splitval) of the data */
      while( !(*perfect) && next_split(left_sorted_sets, right_sorted_sets, p, data_xp, num_cols_x, workspace, 
            &splitval, NULL, NULL) )
      {
         
         assert( are_sorted_sets((const SORTED_SET**) left_sorted_sets, data_x, num_rows, num_cols_x) );
         assert( are_sorted_sets((const SORTED_SET**) right_sorted_sets, data_x, num_rows, num_cols_x) );

#ifdef VERYVERBOSE
         printf("Working on split value %g for covariate %d.\n", splitval, p);
#endif
         
         /* find optimal tree for left data set */
         find_best_split(left_child, depth-1, (const SORTED_SET**) left_sorted_sets, min_node_size, data_x, data_y,
            num_rows, num_cols_x, num_cols_y, best_actions, worst_actions, workspace, &left_perfect); 

         /* find optimal tree for right data set */
         find_best_split(right_child, depth-1, (const SORTED_SET**) right_sorted_sets, min_node_size, data_x, data_y,
            num_rows, num_cols_x, num_cols_y, best_actions, worst_actions, workspace, &right_perfect); 

         /* tree is perfect if and only if both left and right tree are perfect */
         *perfect = left_perfect && right_perfect;

         /* get reward for this split */
         this_reward = get_reward(left_child) + get_reward(right_child);

         /* if best so far, update */
         if( first_reward || this_reward > best_reward ) 
         {
            best_reward = this_reward;
            record_split(node, p, splitval, best_reward);
            record_best_tree(workspace, node, depth);

            first_reward = 0;
         }
      }
   }

   /* set node to best tree */
   retrieve_best_tree(workspace, node, depth);

#ifdef VERYVERBOSE
   printf("Best split for depth=%d tree for dataset of size %d is split value %g for covariate %d with reward %g.\n", depth, get_size(sorted_sets), get_value(node), get_index(node), get_reward(node));
#endif

}



/** Find an optimal policy tree of given maximal depth from a non-empty dataset 
 * @return An optimal policy tree
 */
NODE* tree_search_simple(
  int                    depth,              /**< (maximum) depth of returned tree */
  int                    min_node_size,      /**< smallest terminal node size */
  const double*          data_x,             /**< covariates, data_x+(j*num_rows) points to values for covariate j */
  const double*          data_y,             /**< gammas, data_y+(d*num_rows) points to values for reward d */
  int                    num_rows,           /**< number of units in full dataset */
  int                    num_cols_x,         /**< number of covariates */
  int                    num_cols_y          /**< number of rewards/actions */
  )
{
   NODE* tree = NULL;
   SORTED_SET** initial_sorted_sets = NULL;
   WORKSPACE* workspace = NULL;
   int* best_actions = NULL;
   int* worst_actions = NULL;
   int perfect;

   /* Permitted: depth 0 trees, no covariates */
   /* Not permitted: no data, min_node_size < 1, number of actions < 1 */ 
   
   assert( depth >= 0 );
   assert( min_node_size >= 1 );
   assert( num_rows >= 1 );
   assert( num_cols_x >= 0 );
   assert( num_cols_y >= 1 );
   assert( num_cols_x == 0 || data_x != NULL );
   assert( data_y != NULL );

   /* if no covariates are supplied then no splits are possible, so fix depth of tree to 0 */
   if( num_cols_x == 0 )
      depth = 0;
   
   /* make tree of right depth, values for its nodes will be overwritten with those defining an optimal tree */
   tree = make_tree(depth);

   /* for each of the num_cols_x covariates (if any), use the values of that covariate to sort the
    * (indices of the) datapoints and store the result   
    */
   initial_sorted_sets = make_initial_sorted_sets(data_x, num_rows, num_cols_x);
   assert( are_sorted_sets( (const SORTED_SET**) initial_sorted_sets, data_x, num_rows, num_cols_x) );
   
   /* create working spaces of various sorts (trees, sorted sets, arrays of rewards, etc) */
   workspace = make_workspace(depth, initial_sorted_sets, num_rows, num_cols_x, num_cols_y);

   /* compute and store best and worst actions for each unit */
   store_best_worst_actions(data_y, num_rows, num_cols_y, &best_actions, &worst_actions);

   assert( best_actions != NULL );
   assert( worst_actions != NULL );
   
   /* record that we have yet to find a 'perfect' tree for this data */
   perfect = 0;

   /* find the optimal tree */
   find_best_split(tree, depth, (const SORTED_SET**) initial_sorted_sets, min_node_size, data_x, data_y,
      num_rows, num_cols_x, num_cols_y, best_actions, worst_actions, workspace, &perfect); 

   /* free memory */
   free_workspace(workspace, depth, num_cols_x);
   free_sorted_sets(initial_sorted_sets, num_cols_x);

   /* remove any nodes below leaves, and merge leaves with the same action */
   /* fix_tree(tree); */

   /* return optimal tree */
   return tree;
}
