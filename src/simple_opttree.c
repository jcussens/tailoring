/** @file simple_opttree.c
 *  @brief Functions for building a tree
 *  @author James Cussens
 */

#include "simple_opttree.h"
#include "units.h"
#include "workspace.h"
#include <assert.h>
#include <stdlib.h>
/* #include <stdio.h> */

/* #define VERYVERBOSE  */

#ifdef VERYVERBOSE
#define VERBOSE
#endif

#ifdef VERBOSE
#include <stdio.h>
#endif

/** TODO */
static
double update_bestworstdiff(
   const ELEMENT*        elts,               /**< elements */
   int                   nelts,              /**< number of elements */
   int                   num_rows,           /**< number of units */
   const double*         data_y,             /**< gammas, data_y+(d*num_rows) points to values for reward d */
   const int*            best_actions,       /**< best actions */
   const int*            worst_actions       /**< worst actions */
   )
{
   int i;
   double res = 0.0;

   for( i = 0; i < nelts; i++ )
   {
      ELEMENT elt = elts[i];
      const double* tmp = data_y + elt;
      res += (*(tmp+best_actions[elt]*num_rows) - *(tmp+worst_actions[elt]*num_rows)); 
   }
   return res;
}



/**< update the total reward for each action for a set of 'left' units due to new units being moved into this 'left' set */
static
void update_left_rewards(
   double*               left_rewards,       /**< rewards for each action for a 'left' set of units */
   ELEMENT*              elts,               /**< units just moved into 'left' set of units */
   int                   nelts,              /**< number of units just moved into 'left' set of units */
   const double*         data_y,             /**< data_y[d*num_rows+elt] is the reward for action d for unit elt */
   int                   num_rows,           /**< number of units in full dataset */
   int                   num_cols_y          /**< number of actions */
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

/** find the best action and associated reward for a set of 'right' units */ 
static
void find_best_right_action
(
   const double*         left_rewards,       /**< rewards for each action for a 'left' set of units */
   const double*         nosplit_rewards,    /**< rewards for each action for the set of 'left' and 'right' units combined */
   int                   num_cols_y,         /**< number of actions */
   double*               best_right_reward,  /**< *best_right_reward will be the reward associated with the best action */
   int*                  best_right_action   /**< *best_right_action will be the best action */
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

/** find the best action and associated reward for a set of 'left' units */ 
static
void find_best_left_action(
   const double*         left_rewards,       /**< rewards for each action for a 'left' set of units */
   int                   num_cols_y,         /**< number of actions */
   double*               best_left_reward,   /**< *best_left_reward will be the reward associated with the best action */
   int*                  best_left_action    /**< *best_left_action will be the best action */
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

/** For each unit find and record the best and worst actions for that unit */
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
   CONST_UNITS           units,              /**< the units */
   const double*         data_x,             /**< covariates, data_x+(j*num_rows) points to values for covariate j */
   const double*         data_y,             /**< gammas, data_y+(d*num_rows) points to values for reward d */
   int                   num_rows,           /**< number of units in full dataset */
   int                   num_cols_x,         /**< number of covariates */
   int                   num_cols_y,         /**< number of actions */
   const int*            best_actions,       /**< best_actions[i] is the best action for unit i */
   const int*            worst_actions,      /**< worst_actions[i] is the worst action for unit i */
   WORKSPACE*            workspace           /**< workspace */
   )
{
   int p;
   double best_reward;
   int first_reward = 1;
   double* nosplit_rewards = get_rewards_space(workspace);
   NODE* left_child;
   NODE* right_child;
   
   ELEMENT* elts;
   int nelts;

   double best_left_reward;
   double best_right_reward;
   int best_left_action;
   int best_right_action;

#ifdef VERYVERBOSE
   double bestsplitval;
   int bestp;
#endif

   UNITS right_units;
   int optimal_tree_found;
   
   /* get best possible reward */
   double rewardub = reward_ub(units, data_y, num_rows, best_actions);
   
   /* get reward for each action if no split were done */
   find_nosplit_rewards(units, num_cols_y, data_y, num_rows, nosplit_rewards);

   get_children(node, &left_child, &right_child);
   
   /* consider each covariate for splitting */
   optimal_tree_found = 0;
   for( p = 0; !optimal_tree_found && p < num_cols_x; p++)
   {
      const double* data_xp = data_x+(p*num_rows);
      double splitval;
      
      /* initialise all left rewards to 0 */
      double* left_rewards = get_rewards_space_zeroed(workspace);

      /* initialise the index which will specify the current split */
      int idx = 0;

      /* initialise so that right_units is a copy of units 
         ready for splitting on covariate p */
      shallow_initialise_units(units, p, num_cols_x, workspace, &right_units);

      assert( units_ok((CONST_UNITS) right_units, p, data_x, num_rows, num_cols_x) );

      /* consider each split (x[p] <= splitval, x[p] > splitval) of the data 
         elts[0] ... elts[nelts-1] are the units moved from right to left */
      while( next_shallow_split( (CONST_UNITS) right_units, p, idx, data_xp, &splitval, &elts, &nelts) )
      {
         double this_reward;
         
         update_left_rewards(left_rewards, elts, nelts, data_y, num_rows, num_cols_y);
         find_best_left_action(left_rewards, num_cols_y, &best_left_reward, &best_left_action);
         find_best_right_action(left_rewards, nosplit_rewards, num_cols_y, &best_right_reward, &best_right_action);

         /* get reward for this split */
         this_reward = best_left_reward + best_right_reward;

         assert( this_reward <= rewardub );

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

            /* if this reward is best possible, just stop searching */
            if( this_reward >= rewardub )
            {
               optimal_tree_found = 1;
               break;
            }
         }

         /* update index to move split along */
         idx += nelts;
      }
   }
#ifdef VERYVERBOSE
   printf("Best split for depth=1 tree is split value %g for covariate %d with reward %g.\n",
      bestsplitval, bestp, best_reward);
#endif
}

/** on return, `node` will be the root of an optimal tree of depth `depth` for the data
 * represented by `units`
 */
static
void find_best_split(
   NODE*                 node,               /**< uninitialised tree to be populated with optimal tree */
   int                   depth,              /**< depth of tree */
   CONST_UNITS           units,              /**< the units */
   int                   min_node_size,      /**< smallest terminal node size */
   const double*         data_x,             /**< covariates, data_x+(j*num_rows) points to values for covariate j */
   const double*         data_y,             /**< gammas, data_y+(d*num_rows) points to values for reward d */
   int                   num_rows,           /**< number of units in full dataset */
   int                   num_cols_x,         /**< number of covariates */
   int                   num_cols_y,         /**< number of actions */
   const int*            best_actions,       /**< best_actions[i] is the best action for unit i */
   const int*            worst_actions,      /**< worst_actions[i] is the worst action for unit i */
   WORKSPACE*            workspace,          /**< working space */
   int                   reward_cutoff_set,  /**< whether a reward cutoff has been set */
   double                reward_cutoff,      /**< if reward_cutoff_set, only interested in trees with reward above this value. it's OK to abort search
                                              * if we know we will not find a tree above this value */
   int*                  tree_set            /**< *tree_set=1 if search for an optimal tree was not aborted */
  )
{

   int p;
   int pure;
   double best_reward;
   double best_possible_reward;
   int best_reward_set = 0;
   NODE* left_child;
   NODE* right_child;

   UNITS left_units;
   UNITS right_units;

   int optimal_tree_found;
   
   assert( node != NULL );
   assert( depth >= 0 );
   assert( units != NULL);
   assert( min_node_size >= 1 );
   assert( num_rows >= 1 );
   assert( num_cols_x >= 0 );
   assert( num_cols_y >= 1 );
   assert( num_cols_x == 0 || data_x != NULL );
   assert( data_y != NULL );
   assert( units_ok(units, -1, data_x, num_rows, num_cols_x) );

#ifdef VERBOSE
   printf("Looking for an optimal depth=%d tree for a dataset of size %d.\n", depth, get_size(units));
#endif

   /* get quick poor bound on best possible reward */
   best_possible_reward = reward_ub(units, data_y, num_rows, best_actions);

   /* if no chance of exceeding the cutoff, just abort */
   if( reward_cutoff_set && best_possible_reward <= reward_cutoff )
   {
      /* printf("size=%d, depth=%d, best_possible=%g, cutoff=%g\n", get_size(units), depth, best_possible_reward, reward_cutoff); */
      *tree_set = 0;
      return;
   }
   
   /* determine whether the dataset is pure, i.e. whether each unit has same best action */
   pure = is_pure(units, best_actions);

   /* nothing further to do if depth limit reached or if too few datapoints for splitting or if dataset is pure */
   if( depth == 0 || get_size(units) <= min_node_size || pure )
   {
      int best_action;

      /* find best action and its associated reward */
      find_best_action(units, data_y, num_rows, num_cols_y, workspace, &best_reward, &best_action);

      /* make node a leaf with found best action and associated reward */
      make_leaf(node, best_reward, best_action);

      if( reward_cutoff_set && get_reward(node) <= reward_cutoff )
         *tree_set = 0;
      else
         *tree_set = 1;

      /* all done */
      return;
   }

   /* use specialised routine if only one split allowed */
   if( depth == 1 )
   {
      level_one_learning(node, units, data_x, data_y,
         num_rows, num_cols_x, num_cols_y, best_actions, worst_actions, workspace); 

      if( reward_cutoff_set && get_reward(node) <= reward_cutoff )
         *tree_set = 0;
      else
         *tree_set = 1;
      
      return;
   }

   /* get left and right subtrees */
   get_children(node, &left_child, &right_child);

   assert( left_child != NULL );
   assert( right_child != NULL );

   /* consider each covariate for splitting (stopping if we find a provably optimal tree) */
   optimal_tree_found = 0;
   for( p = 0; !optimal_tree_found && p < num_cols_x; p++)
   {
      const double* data_xp = data_x+(p*num_rows);

      double this_reward;
      double last_reward;

      double splitval;

      int have_last_reward = 0;
      double bestworstdiff = 0.0;
      ELEMENT* elts;
      int nelts;
      
      /* initialise so that left_units is empty and right_units is a copy of units 
         ready for splitting on covariate p */
      initialise_units(units, p, depth, num_cols_x, workspace, &left_units, &right_units);

      assert( units_ok((CONST_UNITS) left_units, p, data_x, num_rows, num_cols_x) );
      assert( units_ok((CONST_UNITS) right_units, p, data_x, num_rows, num_cols_x) );

      /* consider each split (x[p] <= splitval, x[p] > splitval) of the data */
      while( !optimal_tree_found && next_split(left_units, right_units, p, data_xp, num_cols_x, &splitval, &elts, &nelts) )
      {

         int left_reward_cutoff_set = 0;
         int right_reward_cutoff_set = 0;
         int left_tree_set = 0;
         int right_tree_set = 0;

         /* set to dummy values */
         double left_reward_cutoff = 0.0;
         double right_reward_cutoff = 0.0;
         
         assert( units_ok((CONST_UNITS) left_units, p, data_x, num_rows, num_cols_x) );
         assert( units_ok((CONST_UNITS) right_units, p, data_x, num_rows, num_cols_x) );

         if( have_last_reward )
         {
            bestworstdiff += update_bestworstdiff(elts, nelts, num_rows, data_y, best_actions, worst_actions);
            if( last_reward + bestworstdiff <= best_reward )
            {
               continue;
            }
         }

         
#ifdef VERYVERBOSE
         printf("Working on split value %g for covariate %d.\n", splitval, p);
#endif
         
         /* find optimal tree for left data set */
         find_best_split(left_child, depth-1, (CONST_UNITS) left_units, min_node_size, data_x, data_y,
            num_rows, num_cols_x, num_cols_y, best_actions, worst_actions, workspace, left_reward_cutoff_set, left_reward_cutoff, &left_tree_set);

         /* no tree from this split */
         if( !left_tree_set )
         {
            /* printf("left tree aborted\n"); */
            continue;
         }
         
         if( best_reward_set )
         {
            right_reward_cutoff = best_reward - get_reward(left_child);
            right_reward_cutoff_set = 1;
         }
         else if( reward_cutoff_set )
         {
            right_reward_cutoff = reward_cutoff - get_reward(left_child);
            right_reward_cutoff_set = 1;
         }
         else
         {
            right_reward_cutoff_set = 0;
         }
         
         /* find optimal tree for right data set */
         find_best_split(right_child, depth-1, (CONST_UNITS) right_units, min_node_size, data_x, data_y,
            num_rows, num_cols_x, num_cols_y, best_actions, worst_actions, workspace, 
            right_reward_cutoff_set, right_reward_cutoff, &right_tree_set); 

         /* no tree from this split */
         if( !right_tree_set )
         {
            /* printf("right tree aborted\n"); */
            continue;
         }

         /* get reward for this split */
         this_reward = get_reward(left_child) + get_reward(right_child);

         /* record that reward for this split was found and reset bestworstdiff */
         last_reward = this_reward;
         have_last_reward = 1;
         bestworstdiff = 0.0;
         
         /* if best so far, update */
         if( !reward_cutoff_set || this_reward > reward_cutoff )
         {
            if( !best_reward_set || this_reward > best_reward ) 
            {
               best_reward = this_reward;
               record_split(node, p, splitval, best_reward);
               record_best_tree(workspace, node, depth);
               
               best_reward_set = 1;

               /* if there can be no better tree, stop looking for one!
                * use ">=" rather than"==" since should reduce numerical problems */
               if( best_reward >= best_possible_reward )
                  optimal_tree_found = 1;                 
            }
         }
      }
   }

   if( best_reward_set )
   {
      /* set node to best tree */
      retrieve_best_tree(workspace, node, depth);
      *tree_set = 1;
   }
   else
   {
      /* indicate that no acceptable tree found */
      *tree_set = 0;
   }

#ifdef VERYVERBOSE
   printf("Best split for depth=%d tree for dataset of size %d is split value %g for covariate %d with reward %g.\n",
      depth, get_size(units), get_value(node), get_index(node), get_reward(node));
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
  int                    num_cols_y,         /**< number of actions */
  double*                reward              /**< reward for optimal tree */
  )
{
   NODE* tree = NULL;
   UNITS units = NULL;
   WORKSPACE* workspace = NULL;
   int* best_actions = NULL;
   int* worst_actions = NULL;

   int reward_cutoff_set = 0;
   double reward_cutoff = 0.0; /* dummy value */
   int aborted;
   
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

   /* make initial set of units from covariate data */
   units = make_units(data_x, num_rows, num_cols_x);
   assert( units_ok( (CONST_UNITS) units, -1, data_x, num_rows, num_cols_x) );
   
   /* create working spaces of various sorts (trees, units, arrays of rewards, etc) */
   workspace = make_workspace(depth, (CONST_UNITS) units, num_rows, num_cols_x, num_cols_y);

   /* compute and store best and worst actions for each unit */
   store_best_worst_actions(data_y, num_rows, num_cols_y, &best_actions, &worst_actions);

   assert( best_actions != NULL );
   assert( worst_actions != NULL );
   
   /* find the optimal tree */
   find_best_split(tree, depth, (CONST_UNITS) units, min_node_size, data_x, data_y,
      num_rows, num_cols_x, num_cols_y, best_actions, worst_actions, workspace, 
      reward_cutoff_set, reward_cutoff, &aborted); 

   *reward = get_reward(tree);
   
   /* free memory */
   free(best_actions);
   free(worst_actions);
   free_workspace(workspace, depth, num_cols_x);
   free_units(units, num_cols_x);

   /* remove any nodes below leaves, and merge leaves with the same action */
   fix_tree(tree);

   /* return optimal tree */
   return tree;
}
