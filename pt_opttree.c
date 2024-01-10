#include <float.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include "pt_opttree.h"

#define INF DBL_MAX
#define LEFT 0
#define RIGHT 1
#define MIN(a,b) (((a)<(b))?(a):(b))
#define VERBOSE 0
#define VERY_VERBOSE 0
#define USE_PERFECT 1
#define USE_PREVIOUS 1
#define USE_BOUNDS 1
#define DEBUG 0

static
int check_perfect_pt(
   const NODE*           tree,               /**< allegedly perfect tree */
   const SORTED_SET*     sorted_set,         /**< data set for tree */
   const double*         data_x,             /**< covariates, data_x+(j*num_rows) points to values for covariate j */
   const int*            best_actions,       /**< best_actions[i] is the best action for unit i */
   int                   num_rows            /**< number of units */
   )
{
   return check_perfect(tree, get_size(sorted_set), get_elements(sorted_set), data_x, best_actions, num_rows);
}


/**
 * For each unit find and record the best and worst actions for that unit
 */
static
void store_best_worst_actions(
   const double*         data_y,             /**< rewards (column major) */
   int                   num_rows,           /**< number of units */
   int                   num_cols_y,         /**< number of rewards */
   int*                  best_actions,       /**< best actions */
   int*                  worst_actions       /**< worst actions */
  )
{
   int i;
   int best_action;
   int worst_action;
   int d;
   double best_reward;
   double worst_reward;

   assert(num_cols_y > 1);
   
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
      best_actions[i] = best_action;
      worst_actions[i] = worst_action;
      assert(worst_action != best_action);
      assert(worst_reward <= best_reward);
   }
}

/** find best action and associated reward for a set of units 
 * (i.e. a leaf)
 */
static
void find_best_reward(
   const int*            elements,           /**< set of units */
   int                   n,                  /**< number of units */
   const double*         data_y,             /**< rewards (column major) */
   const int             num_cols_y,         /**< number of actions/rewards */
   const int             num_rows,           /**< number of units in full dataset */
   double*               rewards,            /**< temp space to store num_cols_y rewards */
   double*               best_reward,        /**< on return *best_reward is the best reward */ 
   int*                  best_action         /**< if *best_action != -1 then 
                                              *best_action is known to be the best action, on return *best_action is the best action */ 
  )
{
  int d;
  const double* dyelt;
  int i;
  const double* offset;
  
  assert(elements != NULL);
  assert(data_y != NULL);
  assert(best_reward != NULL);
  assert(best_action != NULL);
  assert(num_cols_y >= 0);
  assert(rewards != NULL);
  assert(best_reward != NULL);
  assert(best_action != NULL);
  assert(*best_action >= -1);
  assert(*best_action < num_cols_y);
     
  /* special case when best action is already known */
  if( *best_action != -1 )
  {
     *best_reward = 0;
     offset = data_y + (*best_action)*num_rows;
     for( i = 0; i < n; i++ )
     {
        *best_reward += *(offset + elements[i]);
     }
     return;
  }
  
  for( d = 0; d < num_cols_y; d++ )
    rewards[d] = 0.0;

  /* for each element in set, update the
   * reward value for each possible action
   */
  for( i = 0; i < n; i++ )
  {
    dyelt = data_y + elements[i];
    for( d = 0; d < num_cols_y; d++ )
    {
        rewards[d] += *dyelt;
        dyelt += num_rows;
    }
  }

  *best_reward = rewards[0];
  *best_action = 0;
  for( d = 1; d < num_cols_y; d++ )
    if( rewards[d] > *best_reward )
    {
      *best_reward = rewards[d];
      *best_action = d;
    }
}


/** Find an optimal depth=1 (ie at most one single split) tree for a set of units 
 * If no split is optimal, then that is returned, ie a leaf is returned
 */
static
void level_one_learning(
   NODE*                 node,               /**< uninitialised tree to be populated with optimal tree */
   SORTED_SET**          sorted_sets,        /**< sorted sets for the units, one for each covariate */
   const int             split_step,         /**< consider splits every split_step'th possible split */
   const double*         data_x,             /**< covariates, data_x+(j*num_rows) points to values for covariate j */
   const double*         data_y,             /**< gammas, data_y+(d*num_rows) points to values for reward d */
   const int             num_rows,           /**< number of units in full dataset */
   const int             num_cols_x,         /**< number of covariates */
   const int             num_cols_y,         /**< number of rewards/actions */
   const int*            best_actions,       /**< best_actions[i] is the best action for unit i */
   const int*            worst_actions,      /**< worst_actions[i] is the worst action for unit i */
   double*               rewards,            /**< temporary storage for computing best rewards */
   double*               rewards2,           /**< temporary storage for computing best rewards */
   int*                  perfect             /**< *perfect=1 iff each unit assigned its best action, else *perfect=0 */
   )
{

  int d;
  double* nosplit_rewards = rewards2;
  double* left_rewards = rewards;

  double best_left_reward_for_split;
  int best_left_action_for_split = -1;
  double best_right_reward_for_split;
  int best_right_action_for_split = -1;

  double best_left_reward;
  int best_left_action = -1;
  double best_right_reward;
  int best_right_action = -1;
  double best_reward;
  int best_action;

  int elt;
  int p;

  int i;
  const double* dyelt;

  int best_split_var;
  double best_split_val;

  int n_left;
  int n_right;

  int left_perfect;
  int left_perfect_action;

  SORTED_SET* sorted_set0 = sorted_sets[0];

  assert(node != NULL);
  assert(has_bothchildren(node));
  assert(data_x != NULL);
  assert(data_y != NULL);
  assert(sorted_sets != NULL);

  find_nosplit_rewards(sorted_set0, data_y, num_rows, nosplit_rewards);
  
  /* find best reward if no split were done */
  /* we do not check whether sorted_sets[0] is pure, since this should
     have been done already */
  best_reward = nosplit_rewards[0];
  best_action = 0;
  for( d = 1; d < num_cols_y; d++ )
  {
    if( nosplit_rewards[d] > best_reward )
    {
      best_reward = nosplit_rewards[d];
      best_action = d;
    }
  }
  
  if( VERBOSE )
     printf("No split best reward for %d datapoints is %g, best action is %d\n", get_size(sorted_set0), best_reward, best_action); 

  if( USE_PERFECT )
  {
     /* search for a perfect split */
     for( p = 0; p < num_cols_x; p++)
     {
        const SORTED_SET* sorted_setp = sorted_sets[p];
        const double* data_xp = data_x+(p*num_rows);
        
        if( perfect_split(sorted_setp, best_actions, data_xp, data_y, num_rows, &best_split_val,
              &best_left_reward, &best_left_action, &best_right_reward, &best_right_action) );
        {
           /* we found a perfect split, set values for node */
           record_level_one_split(node, p, best_split_val, best_left_reward + best_right_reward,
              best_left_reward, best_left_action, best_right_reward, best_right_action);

           /* do a check, if debugging */
           assert(check_perfect_pt(node, sorted_setp, data_x, best_actions, num_rows));
           
           if( VERBOSE )
              printf("Found perfect split for depth-1 tree with %d datapoints, covariate=%d, split value=%g, reward=%g=%g+%g .\n",
                 get_size(sorted_setp), p, best_split_val,
                 best_left_reward+best_right_reward,best_left_reward,best_right_reward );

           *perfect = 1;
           return;
        }
     }
  }  

  *perfect = 0;
  /* search for best split */
  for( p = 0; p < num_cols_x; p++)
  {
    const SORTED_SET* sorted_set = sorted_sets[p];
    const double* data_xp = data_x+(p*num_rows);
    
    n_left = 0;
    n_right = sorted_set->n;

    /* initialise left rewards for this p */
    for( d = 0; d < num_cols_y; d++ )
      left_rewards[d] = 0.0;

    /* don't move last element from right to left since then we would have no split */
    for( i = 0; i < (sorted_set->n)-1; i++ )
    {
      elt = sorted_set->elements[i];

      /* update left rewards */
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

      n_left++;
      n_right--;

      
      if( VERY_VERBOSE )
         printf("covariate %d,moved=%d,n_left=%d,action 0 reward=%g\n",p,elt,n_left,left_rewards[0]); 

      assert( data_xp[elt] <= data_xp[sorted_set->elements[i+1]] );

      /* if a proper split see whether it's a best split */
      if( data_xp[elt] != data_xp[sorted_set->elements[i+1]] )
      {
         if( VERY_VERBOSE )
            printf("valid split at %d, since %d has different %d value\n",elt,sorted_set->elements[i+1],p); 

         /* if previous best left action is best action for elt then this remains the best_left_action_for_split */
         if( USE_PREVIOUS && best_left_action_for_split == best_actions[elt] )
         {
            assert(best_left_action_for_split != -1);
            best_left_reward_for_split = left_rewards[best_actions[elt]];
         }
         else
         {
            best_left_reward_for_split = left_rewards[0];
            best_left_action_for_split = 0;
            
            for( d = 1; d < num_cols_y; d++ )
            {
               if( left_rewards[d] > best_left_reward_for_split )
               {
                  best_left_reward_for_split = left_rewards[d];
                  best_left_action_for_split = d;
               }
            }
         }

         /* if previous best right action is worst action for elt then this remains the best_right_action_for_split */
         if( USE_PREVIOUS && best_right_action_for_split == worst_actions[elt] )
         {
            assert(best_right_action_for_split != -1);
            d = worst_actions[elt];
            best_right_reward_for_split = nosplit_rewards[d] - left_rewards[d];
         }
         else
         {
            best_right_reward_for_split = nosplit_rewards[0] - left_rewards[0];
            best_right_action_for_split = 0;
            for( d = 1; d < num_cols_y; d++ )
            {
               if( nosplit_rewards[d] - left_rewards[d] > best_right_reward_for_split )
               {
                  best_right_reward_for_split = nosplit_rewards[d] - left_rewards[d];
                  best_right_action_for_split = d;
               }
            }
         }

        /* update if new best split */
        if( best_left_reward_for_split + best_right_reward_for_split > best_reward)
        {
          best_reward = best_left_reward_for_split + best_right_reward_for_split;
          best_left_reward = best_left_reward_for_split;
          best_left_action = best_left_action_for_split;
          best_right_reward = best_right_reward_for_split;
          best_right_action = best_right_action_for_split;
          best_split_var = p;
          /* split value is last covariate value(=vp) on left, so split is xp <= vp */
          best_split_val = data_xp[elt];

          if( VERBOSE )
             printf("New best split for %d=%d+%d datapoints for depth 1 tree is: covariate=%d, split value=%g, reward=%g=%g+%g\n",
                sorted_set->n,n_left,n_right,p,best_split_val,best_reward,best_left_reward,best_right_reward);
        }
      }
      else
      {
         /* if no split possible then we don't have valid best actions to use for next split */
         best_left_action_for_split = -1;
         best_right_action_for_split = -1;
      }
    }
  }
  
  if( best_left_action != -1 )
  {
    /* we found a split which is better than not doing a split */
     record_level_one_split(node, best_split_var, best_split_val, best_reward,
        best_left_reward, best_left_action,
        best_right_reward, best_right_action);
  }
  else
  {
     make_leaf(node, best_reward, best_action);
  }
}

/** find a tree using a greedy approach */
static
void find_greedy_split(
   NODE*                 node,               /**< uninitialised tree to be populated with optimal tree */
   const int             depth,              /**< depth of tree */
   SORTED_SET**          sorted_sets,        /**< sorted sets for the units, one for each covariate */
   const int             split_step,         /**< consider splits every split_step'th possible split */
   const int             min_node_size,      /**< smallest terminal node size */
   const double*         data_x,             /**< covariates, data_x+(j*num_rows) points to values for covariate j */
   const double*         data_y,             /**< gammas, data_y+(d*num_rows) points to values for reward d */
   const int             num_rows,           /**< number of units in full dataset */
   const int             num_cols_x,         /**< number of covariates */
   const int             num_cols_y,         /**< number of rewards/actions */
   const int*            best_actions,       /**< best_actions[i] is the best action for unit i */
   const int*            worst_actions,      /**< worst_actions[i] is the worst action for unit i */
   NODE***               tmp_trees,          /**< trees of various depths for temporary storage of 'best so far' trees */
   SORTED_SET****        tmp_sorted_sets,    /**< e.g have tmp_sorted_sets[depth][LEFT][p] preallocated space */
   double*               rewards,            /**< temporary storage for computing best rewards */
   double*               rewards2,           /**< temporary storage for computing best rewards */
   int*                  perfect,            /**< *perfect=1 iff each unit assigned its best action, else *perfect=0 */
   int                   root                /**< is this the root node? */ 
  )
{

  int best_action = -1;

  int p;
  int i;
  int j;
  NODE* left_child;
  NODE* right_child;
  NODE* best_left_child;   /* will store the best left child 'so far' */
  NODE* best_right_child;  /* will store the best right child 'so far' */

  int breakpoint_i;
  int breakpoint;
  
  double reward;

  double best_reward = -INF;
  int best_split_var = -1;
  double best_split_val = 0;

  int pp;

  SORTED_SET* left_sorted_set;
  SORTED_SET* right_sorted_set;

  int left_perfect;
  int right_perfect;
  
  assert(node != NULL);
  assert(tmp_trees != NULL);
  assert(data_x != NULL);
  assert(data_y != NULL);
  assert(sorted_sets != NULL);
  assert(depth >= 0);
  assert(num_rows > 0);
  assert(num_cols_x > 0);
  assert(num_cols_y > 0);

  /* we have yet to find a perfect tree for this dataset */
  *perfect = 0;

  if( USE_PERFECT )
  {
     /* determine whether the dataset is pure, i.e. whether each unit has same best action */
     best_action = pure(sorted_sets[0],best_actions);

     /* if this dataset is pure and we can make a perfect leaf */
     if( best_action != -1 )
        *perfect = 1;
  }

  /* nothing further to do for a leaf or if too few datapoints for splitting or if dataset is pure */
  if( depth == 0 || sorted_sets[0]->n <= min_node_size || *perfect )
  {
     /* find best reward with no split 
      * ( if best_action != 1 then we already know the best action and
      * find_best_reward doesn't try to recompute it  
      */
     find_best_reward(sorted_sets[0]->elements,sorted_sets[0]->n,data_y,num_cols_y,num_rows,rewards,&best_reward,&best_action);

    if( VERBOSE && *perfect )
       printf("Node with %d datapoints for depth %d tree is pure with best action %d and reward=%g\n",
          sorted_sets[0]->n,depth,best_action,best_reward);

    /* populate node */
    make_leaf(node, best_reward, best_action);  

    assert(!(*perfect) || check_perfect_pt(node, sorted_sets[0], data_x, best_actions, num_rows) );
    
    return;
  }


  /* find best depth=1 tree for this dataset */
  /* NB this call will set *perfect incorrectly, but it is overwritten later, so no matter */
  level_one_learning(node, sorted_sets, split_step, data_x, data_y, num_rows, num_cols_x, num_cols_y, best_actions,
     worst_actions, rewards, rewards2, perfect);

  /* deal with the case where there is no depth=1 tree is better than a depth=0 tree */
  if( is_leaf(node) )
  {
     /* if the leaf were pure this would have been discovered earlier */
     assert(!(*perfect));

     /* don't choose an arbitrary split and keep going, instead just return the leaf */
     return;
  }
  
  /* initialise all sorted sets for splits of this dataset to be empty on left, and empty on right */
  for( pp = 0; pp < num_cols_x; pp++)
  {
     make_empty(tmp_sorted_sets[depth][LEFT][pp]);
     make_empty(tmp_sorted_sets[depth][RIGHT][pp]);
  }

  /* put each datapoint on either the left or right */

  assert(!is_leaf(node));
  for( i = 0; i < sorted_sets[0]->n; i++ )
  {
     int elt = sorted_sets[0]->elements[i];
     if( data_x[get_index(node)*num_rows+elt] <= get_value(node) )
        for( pp = 0; pp < num_cols_x; pp++)
           insert_element(tmp_sorted_sets[depth][LEFT][pp],elt);
     else
        for( pp = 0; pp < num_cols_x; pp++)
           insert_element(tmp_sorted_sets[depth][RIGHT][pp],elt);
  }

  get_children(node, &left_child, &right_child);
  
  find_greedy_split(left_child, depth-1, tmp_sorted_sets[depth][LEFT], split_step, min_node_size,
     data_x, data_y, num_rows, num_cols_x, num_cols_y, best_actions, worst_actions,
     tmp_trees, tmp_sorted_sets, rewards, rewards2, &left_perfect, 0);

  find_greedy_split(right_child, depth-1, tmp_sorted_sets[depth][RIGHT], split_step, min_node_size,
     data_x, data_y, num_rows, num_cols_x, num_cols_y, best_actions, worst_actions,
     tmp_trees, tmp_sorted_sets, rewards, rewards2, &right_perfect, 0);

  *perfect = left_perfect && right_perfect;

  /* at this point the reward recorded for node will be that for a single split since level_one_learning is used,
   * so now need to update with correct reward
   */
  set_reward(node,get_reward(left_child)+get_reward(right_child));

  assert(!(*perfect) || check_perfect_pt(node, sorted_sets[0], data_x, best_actions, num_rows) );
}


/** on return, `node` will be the root of an optimal tree of depth `depth` for the data
 * represented by `sorted_sets`
 * sorted_sets[0], sorted_sets[1], ... are the same set, they only differ in the order of the elements
 * sorted_sets[p] has its elements sorted according to their rank which is a strict total order
 * consistent with order given by the pth covariate
 * arguments `split_step` to `num_cols_y` are constant values needed to determine optimal tree
 * arguments tmp_trees to rewards2 are pre-allocated working space
 */
static
void find_best_split(
   NODE*                 node,               /**< uninitialised tree to be populated with optimal tree */
   const int             depth,              /**< depth of tree */
   SORTED_SET**          sorted_sets,        /**< sorted sets for the units, one for each covariate */
   const int             split_step,         /**< consider splits every split_step'th possible split */
   const int             min_node_size,      /**< smallest terminal node size */
   const double*         data_x,             /**< covariates, data_x+(j*num_rows) points to values for covariate j */
   const double*         data_y,             /**< gammas, data_y+(d*num_rows) points to values for reward d */
   const int             num_rows,           /**< number of units in full dataset */
   const int             num_cols_x,         /**< number of covariates */
   const int             num_cols_y,         /**< number of rewards/actions */
   const int*            best_actions,       /**< best_actions[i] is the best action for unit i */
   const int*            worst_actions,      /**< worst_actions[i] is the worst action for unit i */
   NODE***               tmp_trees,          /**< trees of various depths for temporary storage of 'best so far' trees */
   SORTED_SET****        tmp_sorted_sets,    /**< e.g have tmp_sorted_sets[depth][LEFT][p] preallocated space */
   double*               rewards,            /**< temporary storage for computing best rewards */
   double*               rewards2,           /**< temporary storage for computing best rewards */
   int*                  perfect,            /**< *perfect=1 iff each unit assigned its best action, else *perfect=0 */
   int                   root                /**< is this the root node? */ 
  )
{

  int best_action = -1;

  int p;
  int i;
  int j;
  NODE* left_child;
  NODE* right_child;
  NODE* best_left_child;   /* will store the best left child 'so far' */
  NODE* best_right_child;  /* will store the best right child 'so far' */

  int breakpoint_i;
  int breakpoint;
  
  double reward;

  double best_reward = -INF;
  int best_split_var = -1;
  double best_split_val = 0;

  int pp;

  SORTED_SET* left_sorted_set;
  SORTED_SET* right_sorted_set;

  int* ml;
  int* mr;

  int n_calls = 0;
  int n_potential_calls = 0;

  NODE* f0;
  
  assert(node != NULL);
  assert(tmp_trees != NULL);
  assert(data_x != NULL);
  assert(data_y != NULL);
  assert(sorted_sets != NULL);
  assert(depth >= 0);
  assert(num_rows > 0);
  assert(num_cols_x > 0);
  assert(num_cols_y > 0);

  /* we have yet to find a perfect tree for this dataset */
  *perfect = 0;

  if( USE_PERFECT )
  {
     /* determine whether the dataset is pure, i.e. whether each unit has same best action */
     best_action = pure(sorted_sets[0],best_actions);

     /* if this dataset is pure and we can make a perfect leaf */
     if( best_action != -1 )
        *perfect = 1;
  }

  /* nothing further to do for a leaf or if too few datapoints for splitting or if dataset is pure */
  if( depth == 0 || sorted_sets[0]->n <= min_node_size || *perfect )
  {
     /* find best reward with no split 
      * ( if best_action != 1 then we already know the best action and
      * find_best_reward doesn't try to recompute it  
      */
     find_best_reward(sorted_sets[0]->elements,sorted_sets[0]->n,data_y,num_cols_y,num_rows,rewards,&best_reward,&best_action);

    if( VERBOSE && *perfect )
       printf("Node with %d datapoints for depth %d tree is pure with best action %d and reward=%g\n",
          sorted_sets[0]->n,depth,best_action,best_reward);

    /* populate node */
    make_leaf(node, best_reward, best_action);

    assert(!(*perfect) || check_perfect_pt(node, sorted_sets[0], data_x, best_actions, num_rows) );
    
    return;
  }

  if( depth == 1 )
  {
     level_one_learning(node, sorted_sets, split_step, data_x, data_y, num_rows, num_cols_x, num_cols_y, best_actions,
        worst_actions, rewards, rewards2, perfect);
    return;
  }

  /* all of these trees are unitialised at this point */
  get_children(node, &left_child, &right_child);
  best_left_child = tmp_trees[depth-1][LEFT];
  best_right_child = tmp_trees[depth-1][RIGHT];
  
  find_greedy_split(node, depth, sorted_sets, split_step, min_node_size,
     data_x, data_y, num_rows, num_cols_x, num_cols_y, best_actions, worst_actions,
     tmp_trees, tmp_sorted_sets, rewards, rewards2, perfect, 0);

  /* this tree is best so far */
  best_reward = get_reward(node);
  best_split_var = get_index(node);
  best_split_val = get_value(node);
  best_action = get_action(node);
  tree_copy(left_child,best_left_child);
  tree_copy(right_child,best_right_child);

  
  if( DEBUG && root )
  {
     printf("Greedily-found tree has reward %g\n", get_reward(node));
  }


  /* find f_{0}, the right tree we get with a dummy split where all units are on the
     right and none on the left
  */

  /* will eventually pre-allocate space for f0, doing it here is temporary */
  
  f0 = make_tree(depth-1);
  find_greedy_split(f0, depth-1, sorted_sets, split_step, min_node_size,
     data_x, data_y, num_rows, num_cols_x, num_cols_y, best_actions, worst_actions,
     tmp_trees, tmp_sorted_sets, rewards, rewards2, perfect, 0);

  if( DEBUG && root )
  {
     printf("Dummy split has reward %g\n", get_reward(f0));
  }

  
  /* similarly need to store elts added/removed from left/right since most recent */
  ml = (int*) malloc (sorted_sets[0]->n * sizeof(int));
  mr = (int*) malloc (sorted_sets[0]->n * sizeof(int));

  /* consider each covariate for splitting */
  for( p = 0; !(*perfect) && p < num_cols_x; p++)
  {

     /* initially assume that any left or right trees are not perfect */
     int left_perfect = 0;
     int right_perfect = 0;

     /* when looking for an optimal tree for a left/right dataset we can sometimes
      * save time if we have the 'previous' left or right tree, ie the one for the previous split
      */
     int have_previous_left_child;
     int have_previous_right_child;

     /* if the previous left or right tree is available sometimes we know it is optimal for the current
        left or right dataset */
     int can_use_previous_left_child;
     int can_use_previous_right_child;
     
     int n_left = 0;
     int n_right = sorted_sets[0]->n;
    
     /* set sorted_set to be the dataset sorted according to covariate p */
     SORTED_SET* sorted_set = sorted_sets[p];
     const double* data_xp = data_x+(p*num_rows);

     /* upper bounds on reward for an optimal tree for current left and right sets */
     double ul;
     double ur;

     /* rewards for most recently found left and right trees */
     double rl;
     double rr;

     /* how many in each of these three sets */
     int nml = 0;
     int nmr = 0;

     int continue_flag;
     
     /* f0 is tree from dummy split, just before first real split 
      * Note that all trees for empty set (on left) are optimal, here we just choose f0 arbitrarily
      */
     tree_copy(f0,left_child);
     tree_copy(f0,right_child);

     /* correct the reward for (empty) left set */
     set_reward(left_child, 0.0);

     /* rl,rr are rewards for most recently found left and right trees, resp. */
     rl = get_reward(left_child);
     rr = get_reward(right_child);

     have_previous_left_child = 1;
     have_previous_right_child = 1;
     
     /* initialise all sorted sets for splits of this dataset to be empty on left, and full on right */
     for( pp = 0; pp < num_cols_x; pp++)
     {
        make_empty(tmp_sorted_sets[depth][LEFT][pp]);
        copy_sorted_set(sorted_sets[pp],tmp_sorted_sets[depth][RIGHT][pp]);
     }
    
     /* don't move last element since then we would have no split */
     for( i = 0; !(*perfect) && i < (sorted_set->n)-1; i++ )
     {
        int elt = sorted_set->elements[i];

        for( pp = 0; pp < num_cols_x; pp++)
        {
           insert_element(tmp_sorted_sets[depth][LEFT][pp],elt);
           remove_element(tmp_sorted_sets[depth][RIGHT][pp],elt);
        }
        
        n_left++;
        n_right--;

        ml[nml++] = elt;
        mr[nmr++] = elt;

        if(DEBUG && root)
           printf("X%d<=%g, i=%d, best=%g, worst=%g\n",p,data_xp[elt],i,data_y[best_actions[elt]*num_rows+elt], data_y[worst_actions[elt]*num_rows+elt]);

        /* check sorted set is indeed sorted for covariate p */
        assert( data_xp[elt] <= data_xp[sorted_set->elements[i+1]] );
        
        /* if a proper split see whether it's a best split */
        if( data_xp[elt] < data_xp[sorted_set->elements[i+1]] )
        {
           /* in the worst case have to make 2 recursive calls */
           n_potential_calls += 2;

           /* if have_previous_left_tree then left_child is an optimal tree for the previous split
            * if this tree assigns each elt in M (those moved since previous left child)
            * its best action then left_child is also optimal
            * for the left part of the current split
            */
           can_use_previous_left_child = 0;
           if( USE_PREVIOUS && have_previous_left_child )
           {
              can_use_previous_left_child = 1;
              for(j = 0; j < nml; j++)
              {
                 elt = ml[j];
                 if (assigned_action(left_child, data_x, num_rows, elt) != best_actions[elt] )
                 {
                    can_use_previous_left_child = 0;
                    break;
                 }
              }

              if( can_use_previous_left_child )
              {
                 for(j = 0; j < nml; j++)
                 {
                    elt = ml[j];
                    update_rewards(left_child, data_x, num_rows, elt, data_y[best_actions[elt]*num_rows+elt]);
                 }
                 /* nothing 'pending' */
                 nml = 0;
                 /* this is now the reward for most recent left tree */
                 set_reward(left_child,rl);
                 /* left child is optimal, so also an upper bound */
                 ul = rl;
              }
           }

           /* as above but for right side */
           can_use_previous_right_child = 0;
           if( USE_PREVIOUS && have_previous_right_child )
           {
              can_use_previous_right_child = 1;
              for(j = 0; j < nmr; j++)
              {
                 elt = mr[j];
                 if (assigned_action(right_child, data_x, num_rows, elt) != worst_actions[elt] )
                 {
                    can_use_previous_right_child = 0;
                    break;
                 }
              }

              if( can_use_previous_right_child )
              {
                 for(j = 0; j < nmr; j++)
                 {
                    elt = mr[j];
                    update_rewards(right_child, data_x, num_rows, elt, -data_y[worst_actions[elt]*num_rows+elt]);
                 }
                 nmr = 0;
                 rr = get_reward(right_child);
                 ur = rr;
              }
           }

           if( DEBUG && root )
              printf("Use previous: %d,%d|%d,%d. ",
                 have_previous_left_child, can_use_previous_left_child,
                 have_previous_right_child,can_use_previous_right_child);

           /* update have_previous_left_child, have_previous_right_child
            * for next loop. At this point we have not computed new trees so
            * have_previous_left_child iff can_use_previous_left_child (similarly for right)
            * values will be overwritten if a left/right tree is computed
            */

           have_previous_left_child = can_use_previous_left_child;
           have_previous_right_child = can_use_previous_right_child;

           if( USE_BOUNDS && !(can_use_previous_left_child  && can_use_previous_right_child) )
           {
              /* if we can't use both previous trees, then we might still be able to use
                 upper bounds to prevent recursive calls */

              if( !can_use_previous_left_child )
              {
                 /* don't already have upper bound for left, so compute one */
                 ul = rl;
                 for( j = 0; j < nml; j++ )
                 {
                    elt = ml[j];
                    ul += data_y[best_actions[elt]*num_rows+elt];
                 }
              }

              if( !can_use_previous_right_child )
              {
                 /* don't already have upper bound for right, so compute one */
                 ur = rr;
                 for( j = 0; j < nmr; j++ )
                 {
                    elt = mr[j];
                    ur -= data_y[worst_actions[elt]*num_rows+elt];
                 }
              }
           }

           
           if(DEBUG && root)
              printf("ub is %g=%g+%g.\n", ul+ur,ul,ur);

           /* just continue if this split cannot beat the best split */
           continue_flag = ( ul + ur <= best_reward );

           if( !continue_flag && !can_use_previous_left_child )
           {
              /* have to do recursive call */
              n_calls++;
              find_best_split(left_child, depth-1, tmp_sorted_sets[depth][LEFT], split_step, min_node_size,
                 data_x, data_y, num_rows, num_cols_x, num_cols_y, best_actions, worst_actions,
                 tmp_trees, tmp_sorted_sets, rewards, rewards2, &left_perfect, 0);
              have_previous_left_child = 1;
              rl = get_reward(left_child);
              nml = 0;
           }

           if( !continue_flag && !can_use_previous_right_child )
           {
              /* have to do recursive call */
              n_calls++;
              find_best_split(right_child, depth-1, tmp_sorted_sets[depth][RIGHT], split_step, min_node_size,
                 data_x, data_y, num_rows, num_cols_x, num_cols_y, best_actions, worst_actions,
                 tmp_trees, tmp_sorted_sets, rewards, rewards2, &right_perfect, 0);
              have_previous_right_child = 1;
              rr = get_reward(right_child);
              nmr = 0;
           }

           /* only check whether we have a new incumbent if we actually have both the left and right tree */
           /* NB have_previous_left_child indicates a left tree found on this iteration, variable has the name
            * to indicate this in next iteration, ditto for have_previous_right_child
            */
           if( !continue_flag )
           {
              reward = get_reward(left_child) + get_reward(right_child);
              if( DEBUG && root )
                 printf("Reward is %g=%g+%g\n", reward, get_reward(left_child), get_reward(right_child));
              
              *perfect = left_perfect && right_perfect;
           
              if( reward > best_reward )
              {
                 best_reward = reward;
                 best_split_var = p;
                 best_split_val = data_xp[elt];
                 /* if( VERBOSE || root ) */
                 /* { */
                 /*    printf("New best split for %d=%d+%d datapoints for depth %d tree is: covariate=%d, split value=%g, reward=%g=%g+%g\n", */
                 /*       sorted_set->n,n_left,n_right,depth,p,best_split_val,reward,left_child->reward,right_child->reward); */
                    
                 /*    if( *perfect ) */
                 /*       printf("Split was perfect, no need to consider any others\n"); */
                 /* } */
                 
                 /* save best left and right trees found so far */
                 /* tree_copy(source,target) */
                 tree_copy(left_child,best_left_child);
                 tree_copy(right_child,best_right_child);
              }
           }
           
        } /* scope is a particular split */
     } /* scope is an element moved from right to left */
  }  /* scope is a particular covariate */

  /* populate node */
  if( best_split_var == -1 )
  {
     make_leaf(node, best_reward, best_action);
  }
  else
  {
     record_split(node, best_split_var, best_split_val, best_reward);
     /* retrieve saved best left and right trees */
     /* tree_copy(source,target) */
     tree_copy(best_left_child,left_child);
     tree_copy(best_right_child,right_child);
  }

  assert(!(*perfect) ||  check_perfect_pt(node, sorted_sets[0], data_x, best_actions, num_rows) );

  free(ml);
  free(mr);
  tree_free(f0);

  /* if(root) */
  /*    printf("\n%d recursive calls executed out of a possible %d recursive calls.\n\n", n_calls, n_potential_calls); */
}

/*
 * Notes: A tree of depth 0 is a single node (with no children) where only the reward and action_id members are meaningful.
 */

/** Find an optimal policy tree of given maximal depth 
 * @return An optimal policy tree
 */
NODE* tree_search_jc_policytree(
  int                    depth,              /**< (maximum) depth of returned tree */
  int                    split_step,         /**< consider splits every split_step'th possible split */
  int                    min_node_size,      /**< smallest terminal node size */
  const double*          data_x,             /**< covariates, data_x+(j*num_rows) points to values for covariate j */
  const double*          data_y,             /**< gammas, data_y+(d*num_rows) points to values for reward d */
  int                    num_rows,           /**< number of units in full dataset */
  int                    num_cols_x,         /**< number of covariates */
  int                    num_cols_y          /**< number of rewards/actions */
  )
{

  NODE* tree;
  int* indices;
  int i;
  int p;
  NODE*** tmp_trees;
  int* tmp_indices;
  double* rewards;
  double* rewards2;

  SIMPLE_SORTED_SET** initial_sorted_sets = (SIMPLE_SORTED_SET**) malloc(num_cols_x*sizeof(SIMPLE_SORTED_SET*));
  
  /* have pre-allocated memory for each depth, each direction, each covariate */
  SIMPLE_SORTED_SET**** tmp_sorted_sets;

  int* best_actions;
  int* worst_actions;

  int perfect;
  double ub;
  
  /* data_x[col * num_rows + row] is the value of covariate 'col' in datapoint
   * 'row', so column major storage, all values for covariate 0 first
   *
   * data_y[col * num_rows + row] is the value of reward 'col' in datapoint
   * 'row', so column major storage, all values of reward 0 first.
   */
  
  /* make skeleton tree of correct depth */
  tree = make_tree(depth);
  tmp_indices = (int*) malloc(num_rows*sizeof(int));
  
  /* make initial sorted sets */
  for( p = 0; p < num_cols_x; p++)
    initial_sorted_sets[p] = make_sorted_set(num_rows,data_x+p*num_rows,tmp_indices); 
  free(tmp_indices);

  /* make temporary trees (for storing 'best so far' trees) */
  tmp_trees = (NODE***) malloc(depth*sizeof(NODE**));
  for(i = 0; i < depth; i++)
  {
    tmp_trees[i] = (NODE**) malloc(2*sizeof(NODE*));
    tmp_trees[i][LEFT] = make_tree(i);
    tmp_trees[i][RIGHT] = make_tree(i);
  }

  /* make temporary sorted sets */
  tmp_sorted_sets = make_sorted_set_array(initial_sorted_sets, depth, num_cols_x);
  
  rewards = (double*) malloc(num_cols_y*sizeof(double));
  rewards2 = (double*) malloc(num_cols_y*sizeof(double));

  best_actions = (int*) malloc(num_rows*sizeof(int));
  worst_actions = (int*) malloc(num_rows*sizeof(int));
  
  store_best_worst_actions(data_y, num_rows, num_cols_y, best_actions, worst_actions);

  ub = 0.0;
  for(i = 0; i < num_rows; i++ )
  {
     ub += data_y[best_actions[i]*num_rows+i];
  }
  /* printf("A perfect tree has reward %g.\n", ub); */
  
  perfect = 0;
  find_best_split(tree, depth, initial_sorted_sets, split_step, min_node_size, data_x, data_y, num_rows, num_cols_x, num_cols_y, best_actions, worst_actions,
     tmp_trees, tmp_sorted_sets, rewards, rewards2, &perfect, 1);  /* tmp_sorted_sets, rewards and rewards2 are temporary reusable storage */

  free(best_actions);
  free(worst_actions);
  
  free(rewards);
  free(rewards2);
    
  /* remove temporary sorted sets */
  free_sorted_set_array(tmp_sorted_sets, depth, num_cols_x);
    
  /* remove temporary trees */
  for(i = 0; i < depth; i++)
  {
    tree_free(tmp_trees[i][0]);
    tree_free(tmp_trees[i][1]);
    free(tmp_trees[i]);
  }
  free(tmp_trees);

  /* free sorted sets */
  for( p = 0; p < num_cols_x; p++)
     free_sorted_set(initial_sorted_sets[p]);
  free(initial_sorted_sets);

  fix_tree(tree);
  
  return tree;
}

