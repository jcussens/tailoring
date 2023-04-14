#include <float.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include "opttree.h"

#define INF DBL_MAX
#define LEFT 0
#define RIGHT 1
#define MIN(a,b) (((a)<(b))?(a):(b))
#define VERBOSE 0
#define VERY_VERBOSE 0


/** 
 * sorted set, policy tree style
*/
struct sorted_set
{
  int*                   elements;           /**< the sorted set */
  int                    n;                  /**< size of subset */
  int*                   rank;               /**< if i < j then rank[elements[i]] < rank[elements[j]]  */
};
typedef struct sorted_set SORTED_SET;

/** Determine whether a sorted set is 'pure'.
 * A pure sorted set is one where each unit has the same best action
 * @return The best action if the set is pure, otherwise -1
 */
static
int pure(
   const SORTED_SET*     sorted_set,         /**< sorted set */
   const int*            best_actions        /**< best_actions[i] is the best action for unit i */
  )
{
  int i;
  int best_action = 0;

  if( sorted_set->n > 0)
     best_action = best_actions[sorted_set->elements[0]];
  
  for( i = 1; i < sorted_set->n; i++)
     if( best_action != best_actions[sorted_set->elements[i]] )
        return -1;

  return best_action;
}

/**
 * For each unit find and record the best action for that unit
 * NB memory pointed to by returned pointer must be freed at some point
 * @return pointer to the array of best actions (one best action for each unit, in order)
 */
static
int* store_best_actions(
  const double*          data_y,             /**< rewards (column major) */
  int                    num_rows,           /**< number of units */
  int                    num_cols_y          /**< number of rewards */
  )
{
  int i;
  int* best_actions = (int*) malloc(num_rows*sizeof(int));
  int best_action;
  int d;
  double best_reward;
  
  for( i = 0; i < num_rows; i++ )
  {
    best_action = 0;
    best_reward = data_y[i];
    
    for( d = 1; d < num_cols_y; d++ )
    {
      if( data_y[d*num_rows+i] > best_reward )
      {
        best_action = d;
        best_reward = data_y[d*num_rows+i];
      }
    }
    best_actions[i] = best_action;
  }
  return best_actions;
}

/** print out a sorted set (for debugging only, not currently used) */
static
void print_set(
   const SORTED_SET*     sorted_set          /**< sorted set */
   )
{
   int i;
   int elt;
   
   for(i = 0; i < sorted_set->n; i++)
   {
      elt = sorted_set->elements[i];
      printf("%d(%d), ",elt,sorted_set->rank[elt]);
   }
   printf("\n");
}


/** remove an element from a sorted set, assuming it is already there 
 * @return 1 if the element is a member of the set, else -1
 */
static
int remove_element(
   SORTED_SET*           sorted_set,         /**< sorted set */
   int                   elt                 /**< element to remove */
  )
{

   int left;
   int right;
   int mid;
   int midsave;
   int midval;
   int eltval;
   int* rank;
   int n;
   int* set;  
   int retval = -1;
   
   assert(sorted_set != NULL);
   
   rank = sorted_set->rank;
   n = sorted_set->n;
   set = sorted_set->elements;  
   eltval = rank[elt];

   left = 0;
   right = n-1;

   while( left <= right )
   {
      mid = (left+right)/2;
      midval = rank[set[mid]];
      
      if( midval < eltval )
         left = mid+1;
      else if ( midval > eltval )
         right = mid - 1;
      else
      {
         retval = 1;
         break;
      }
   }

   if( retval == 1 )
   {
      /* unless the last element, have to shift elements */
      if( mid < n-1)
         memmove(set+mid,set+mid+1,(n-mid-1)*sizeof(int));
      sorted_set->n--;
   }
   
   return retval;
}

/** insert an element into a sorted set, assuming it is not already there */
static
void insert_element(
   SORTED_SET*           sorted_set,         /**< sorted set */
   int                   elt                 /**< element to insert */ 
  )
{
  int left;
  int right;
  int mid;
  int eltval;
  int* rank;
  int n;
  int* set;  
  
  assert(sorted_set != NULL);
  assert(elt >= 0);
  
  rank = sorted_set->rank;
  n = sorted_set->n;
  set = sorted_set->elements;  
  eltval = rank[elt];

  assert(rank != NULL);
  assert(set != NULL);
  assert(n >= 0);
  
  /* if set not empty look for an index 'mid' such that
   * rank[set[mid]] < eltval < rank[set[mid+1]]
   * can then insert elt by shifting set[mid+1], set[mid+2],,, and then
   * doing set[mid+1] = elt
   */

  if( n == 0)
  {
     set[0] = elt;
     sorted_set->n = 1;
     return;
  }
  else if( eltval < rank[set[0]] )
  {
    /* insert at beginning */
    mid = -1;
  }
  else if ( eltval > rank[set[n-1]] )
  {
    /* insert at end */
    mid = n-1;
  }
  else
  {
     /* since set is ordered there will be a value mid such that
      * rank[set[mid]] < eltval < rank[set[mid+1]]
      * where -1 < mid < n-1 (i.e. mid+1 is the index of an existing element )
      */

     assert(n > 1);
     left = 0;                /* smallest possible value for mid */
     right = n - 2;           /* largest possible value for mid */

     while( left <= right )
     {
        mid = (left + right)/2;  /* left <= mid <= right */
        
        if( rank[set[mid]] > eltval )
        {
           /* rank[set[mid]] too big, update biggest possible value for mid */
           right = mid - 1;
        }
        else if ( rank[set[mid+1]] < eltval )
        {
           /* rank[set[mid+1]] too small, update smallest possible value for mid*/
           left = mid + 1;
        }
        else
           break;
     }
  }
  
  /* move, if necessary */
  if( mid < n )
    memmove(set+mid+2,set+mid+1,(n-mid-1)*sizeof(int));
  
  /* insert at mid */
  set[mid+1] = elt;
  sorted_set->n++;
}


/** print a policy tree (debugging only) */
static
void print_tree(
   const NODE*           tree,               /**< root of tree to print */ 
   const char**          covnames            /**< covnames[i] is the name of covariate i */
  )
{
  printf("node = %p\n", (void*) tree);
  if( tree-> index != -1)
    printf("covariate = %s\n", covnames[tree->index]);
  printf("value = %g\n", tree->value);
  printf("reward = %g\n", tree->reward);
  printf("action_id = %d\n", tree->action_id);
  printf("left_child = %p\n", (void*) tree->left_child);
  printf("right_child = %p\n", (void*) tree->right_child);
  printf("\n");

  if(tree->left_child != NULL)
    print_tree(tree->left_child,covnames);

  if(tree->right_child != NULL)
    print_tree(tree->right_child,covnames);
}


/** merge a pair of ordered subsequences to get a single ordered subsequence
 * left run is indices[ileft:iright-1]
 * right run is indices[iright:iend-1]
 */
static
void bottomupmerge(
   const int*            indices,            /**< global array */
   int                   ileft,              /**< index of first element of left subsequence */
   int                   iright,             /**< index of first element of right subsequence */
   int                   iend,               /**< iend-1 is index of last element of right subsequence */
   int*                  tmp,                /**< on output tmp[ileft:iend-1] has the merged subsequence */
   const double*         data_x_p            /**< data_x_p[i] is ordering key for i */
  )
{
  int i = ileft;
  int j = iright;
  int k;

  for( k = ileft; k < iend; k++ )
  {
    if( i < iright && (j >= iend || data_x_p[indices[i]] <= data_x_p[indices[j]]) )
      tmp[k] = indices[i++];
    else
      tmp[k] = indices[j++];
  }
}

/** sort an array by bottom up merge sort */
static
void bottomupmergesort(
   int*                  indices,            /**< array to sort */
   int*                  tmp,                /**< on output, sorted array */              
   int                   n,                  /**< length of array */
   const double*         data_x_p            /**< data_x_p[i] is ordering key for i */
  )
{
  int width;
  int i;
  const size_t size = n*sizeof(int);
  
  for( width = 1; width < n; width *= 2)
  {
    for( i = 0; i < n; i += 2*width)
      bottomupmerge(indices, i, MIN(i+width,n), MIN(i+2*width,n), tmp, data_x_p);

    memcpy(indices,tmp,size);
  }
}

/** make a sorted set the empty set */
static
void make_empty(
   SORTED_SET*           sorted_set          /**< sorted set */
  )
{
   sorted_set->n = 0;
}

/** copy a sorted set, overwriting an existing target set */
static
void copy_sorted_set(
   const SORTED_SET*     source,             /**< source sorted set */
   SORTED_SET*           target              /**< target sorted set */
   )
{
   assert( source != NULL );
   assert( target != NULL );
   assert( source->elements != NULL );
   assert( target->elements != NULL );
   assert( source->rank != NULL );
   assert( source->rank == target->rank );
   
   memcpy(target->elements,source->elements,source->n*sizeof(int));
   target->n = source->n; 
}

/** create a new uninitialised subset from existing sorted subset */
static
void new_sorted_subset(
   const SORTED_SET*     source,             /**< source sorted set */
   SORTED_SET*           target              /**< target sorted set */
  )
{

   assert( source != NULL );
   assert( target != NULL );
   assert( source->elements != NULL );
   assert( source->rank != NULL );

   /* ensure (typically more than) enough space */
   target->elements = (int*) malloc((source->n)*sizeof(int));
   target->rank = source->rank;
   /* don't set target->n since size of target not determined yet */
}


/** make a sorted set from elements 0,...,num_indices-1 
 * this function allocates all the necessary memory for the sorted set
 * @return The sorted set
 */
static
SORTED_SET* make_sorted_set(
   const int             num_indices,            /**< size of set */
   const double*         data_xx,                /**< data_xx[i] is ordering key for element i */
   int*                  tmp                     /**< temporary storage for ordering */
  )
{
  SORTED_SET* sorted_set;
  int i;
  int* elements;
  int* rank;

  sorted_set = (SORTED_SET*) malloc(sizeof(SORTED_SET));
  elements = (int *) malloc(num_indices*sizeof(int));
  rank = (int *) malloc(num_indices*sizeof(int));
  for( i = 0; i < num_indices; i++)
    elements[i] = i;

  bottomupmergesort(elements,tmp,num_indices,data_xx);

  for( i = 0; i < num_indices; i++)
    rank[elements[i]] = i;
  
  sorted_set->elements = elements;
  sorted_set->n = num_indices;
  sorted_set->rank = rank;
  
  return sorted_set;
}

/**
 * Make a policy tree of required depth with dummy values,
 * which can later be overwritten
 * @return a policy tree of required depth
 * or NULL if there is insufficient memory
 */
static
NODE* make_tree(
   int                   depth               /**< required depth of tree */
  )
{
  NODE* node;

  assert(depth >= 0); 
  
  node = (NODE*) malloc(sizeof(NODE));

  if( node == NULL )
  {
    /* not enough memory! */
    return NULL;
  }

  /* explicitly set default values for all members to keep valgrind happy */

  node->index = -1; /* no splitting covariate so far, may never be one...*/
  node->value = 0;
  node->reward = 0;
  node->action_id = -1;
  
  if( depth > 0 )
  {
    node->left_child = make_tree(depth-1);
    node->right_child = make_tree(depth-1);
    if( node->left_child == NULL || node->right_child == NULL )
      /* not enough memory */
      return NULL;
  }
  else
  {
    /* for leaf nodes set children explicitly to NULL */
    node->left_child = NULL;
    node->right_child = NULL;
  }
  return node;
}


/**
 * copy data from source to target tree
 * assumes both trees have same depth
 */
static
void tree_copy(
   const NODE*           source,             /**< source tree */
   NODE*                 target              /**< target tree */
  )
{

  assert(source != NULL);
  assert(target != NULL);
  
  target->index = source->index;
  target->value = source->value;
  target->reward = source->reward;
  target->action_id = source->action_id;
  if( source->left_child != NULL)
    tree_copy(source->left_child,target->left_child);
  if( source->right_child != NULL)
    tree_copy(source->right_child,target->right_child);
}

/* /\* */
/*  * delete a tree (free the memory it occupied) */
/*  *\/ */
/* void tree_free( */
/*   NODE* node */
/*   ) */
/* { */
/*   assert(node != NULL); */
  
/*   if( node->left_child != NULL ) */
/*   { */
/*     tree_free(node->left_child); */
/*     tree_free(node->right_child); */
/*   } */
/*   free(node); */
/* } */


/** find best action and associated reward for a set of units 
 * (i.e. a leaf)
 */
static
void find_best_reward(
   const SORTED_SET*     sorted_set,         /**< set of units */
   const double*         data_y,             /**< rewards (column major) */
   const int             num_cols_y,         /**< number of actions/rewards */
   const int             num_rows,           /**< number of units in full dataset */
   double*               rewards,            /**< temp space to store num_cols_y rewards */
   double*               best_reward,        /**< on return *best_reward is the best reward */ 
   int*                  best_action         /**< on return *best_action is the best action */ 
  )
{
  int d;
  const double* dyelt;
  int i;
  
  assert(sorted_set != NULL);
  assert(data_y != NULL);
  assert(best_reward != NULL);
  assert(best_action != NULL);
  assert(num_cols_y >= 0);

  for( d = 0; d < num_cols_y; d++ )
    rewards[d] = 0.0;

  /* for each element in set, update the
   * reward value for each possible action
   */
  for( i = 0; i < sorted_set->n; i++ )
  {
    dyelt = data_y + sorted_set->elements[i];
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

/** Find an optimal depth=1 (ie a single split) tree for a set of units */
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
   double*               rewards,            /**< temporary storage for computing best rewards */
   double*               rewards2            /**< temporary storage for computing best rewards */
   )
{

  int d;
  double* nosplit_rewards = rewards2;
  double* left_rewards = rewards;

  double best_left_reward_for_split;
  int best_left_action_for_split;
  double best_right_reward_for_split;
  int best_right_action_for_split;

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


  assert(node != NULL);
  assert(node->left_child != NULL);
  assert(node->right_child != NULL);
  assert(data_x != NULL);
  assert(data_y != NULL);
  assert(sorted_sets != NULL);

  
  /* find reward for each action if no split were done */
  for( d = 0; d < num_cols_y; d++ )
      nosplit_rewards[d] = 0.0;
  /* use set ordered on first covariate, would get same result
   * if any other covariate were used
   */
  for( i = 0; i < (sorted_sets[0])->n; i++ )
  {
    dyelt = data_y + (sorted_sets[0])->elements[i];
    for( d = 0; d < num_cols_y; d++ )
    {
      nosplit_rewards[d] += *dyelt;
      dyelt += num_rows;
    }
  }
  
  /* find best reward if no split were done */
  best_reward = nosplit_rewards[0];
  best_action = 0;
  for( d = 1; d < num_cols_y; d++ )
    if( nosplit_rewards[d] > best_reward )
    {
      best_reward = nosplit_rewards[d];
      best_action = d;
    }

    if( VERBOSE )
       printf("No split best reward for %d datapoints is %g, best action is %d\n", sorted_sets[0]->n, best_reward, best_action); 
  
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

         
        /* find best left and right reward+action */
        best_left_reward_for_split = left_rewards[0];
        best_left_action_for_split = 0;
        best_right_reward_for_split = nosplit_rewards[0] - left_rewards[0];
        best_right_action_for_split = 0;
        for( d = 1; d < num_cols_y; d++ )
        {
          if( left_rewards[d] > best_left_reward_for_split )
          {
            best_left_reward_for_split = left_rewards[d];
            best_left_action_for_split = d;
          }
          if( nosplit_rewards[d] - left_rewards[d] > best_right_reward_for_split )
          {
            best_right_reward_for_split = nosplit_rewards[d] - left_rewards[d];
            best_right_action_for_split = d;
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
          best_split_val = data_xp[elt];

          if( VERBOSE )
             printf("New best split for %d=%d+%d datapoints for depth 1 tree is: covariate=%d, split value=%g, reward=%g=%g+%g\n",
                sorted_set->n,n_left,n_right,p,best_split_val,best_reward,best_left_reward,best_right_reward);

        }
      }
    }
  }
  
  if( best_left_action != -1 )
  {
    /* we found a split which is better than not doing a split */
    node->index = best_split_var;
    node->value = best_split_val;
    node->reward = best_reward;
    node->left_child->index = -1;
    node->left_child->reward = best_left_reward;
    node->left_child->action_id = best_left_action;
    node->right_child->reward = best_right_reward;
    node->right_child->action_id = best_right_action;
  }
  else
  {
    node->index = -1;
    node->reward = best_reward;
    node->action_id = best_action;
  }
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
   NODE***               tmp_trees,          /**< trees of various depths for temporary storage */
   SORTED_SET****        tmp_sorted_sets,    /**< e.g have tmp_sorted_sets[depth][LEFT][p] preallocated space */
   double*               rewards,            /**< temporary storage for computing best rewards */
   double*               rewards2            /**< temporary storage for computing best rewards */
  )
{

  int best_action;

  int p;
  int i;
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

  int best_reward_for_all;

  
  assert(node != NULL);
  assert(tmp_trees != NULL);
  assert(data_x != NULL);
  assert(data_y != NULL);
  assert(sorted_sets != NULL);
  assert(depth >= 0);
  assert(num_rows > 0);
  assert(num_cols_x > 0);
  assert(num_cols_y > 0);
  

  /* nothing further to do for a leaf or if too few datapoints for splitting or if dataset is pure */
  if( depth == 0 || sorted_sets[0]->n <= min_node_size || pure(sorted_sets[0],best_actions) != -1)
  {
    /* find best reward with no split */
    find_best_reward(sorted_sets[0],data_y,num_cols_y,num_rows,rewards,&best_reward,&best_action);

    if( VERBOSE && pure(sorted_sets[0],best_actions) != -1)
       printf("Node with %d datapoints for depth %d tree is pure with best action %d and reward=%g\n",
          sorted_sets[0]->n,depth,best_action,best_reward);

    
    /* populate node */
    node->index = -1;
    node->reward = best_reward;
    node->action_id = best_action;
    return;
  }

  if( depth == 1 )
  {
     level_one_learning(node, sorted_sets, split_step, data_x, data_y, num_rows, num_cols_x, num_cols_y, best_actions,
      rewards, rewards2);
    return;
  }

  left_child = node->left_child;
  right_child = node->right_child;
  best_left_child = tmp_trees[depth-1][LEFT];
  best_right_child = tmp_trees[depth-1][RIGHT];

  /* consider each covariate for splitting */
  for( p = 0; p < num_cols_x; p++)
  {
     
     int n_left = 0;
     int n_right = sorted_sets[0]->n;

     
    SORTED_SET* sorted_set = sorted_sets[p];
    const double* data_xp = data_x+(p*num_rows);

    /* initialise all sorted sets to be empty on left, and full on right */
    for( pp = 0; pp < num_cols_x; pp++)
    {
      make_empty(tmp_sorted_sets[depth][LEFT][pp]);
      copy_sorted_set(sorted_sets[pp],tmp_sorted_sets[depth][RIGHT][pp]);
    }
    
    /* don't move last one since then we would have no split */
    /* if splitting is not necessary then will just get same best tree on both sides */
    for( i = 0; i < (sorted_set->n)-1; i++ )
    {
      int elt = sorted_set->elements[i];

      for( pp = 0; pp < num_cols_x; pp++)
      {
        insert_element(tmp_sorted_sets[depth][LEFT][pp],elt);
        remove_element(tmp_sorted_sets[depth][RIGHT][pp],elt);
      }

      n_left++;
      n_right--;

      
      assert( data_xp[elt] <= data_xp[sorted_set->elements[i+1]] );

      /* if a proper split see whether it's a best split */
      if( data_xp[elt] < data_xp[sorted_set->elements[i+1]] )
      {
         find_best_split(left_child, depth-1, tmp_sorted_sets[depth][LEFT], split_step, min_node_size, data_x, data_y, num_rows, num_cols_x, num_cols_y, best_actions,
          tmp_trees, tmp_sorted_sets, rewards, rewards2);
         find_best_split(right_child, depth-1, tmp_sorted_sets[depth][RIGHT], split_step, min_node_size, data_x, data_y, num_rows, num_cols_x, num_cols_y, best_actions,
          tmp_trees, tmp_sorted_sets, rewards, rewards2);

        reward = left_child->reward + right_child->reward;
        
        if( reward > best_reward )
        {
          best_reward = reward;
          best_split_var = p;
          best_split_val = data_xp[elt];
          if( VERBOSE )
             printf("New best split for %d=%d+%d datapoints for depth %d tree is: covariate=%d, split value=%g, reward=%g=%g+%g\n",
                sorted_set->n,n_left,n_right,depth,p,best_split_val,reward,left_child->reward,right_child->reward);

          /* save best left and right trees found so far */
          /* tree_copy(source,target) */
          tree_copy(left_child,best_left_child);
          tree_copy(right_child,best_right_child);
        }
      }
    }
  }

  /* populate node */
  node->index = best_split_var;
  node->value = best_split_val;
  node->reward = best_reward;
  if( best_split_var != -1 )
  {
    /* retrieve saved best left and right trees */
    /* tree_copy(source,target) */
    tree_copy(best_left_child,left_child);
    tree_copy(best_right_child,right_child);
  }
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

  SORTED_SET** initial_sorted_sets = (SORTED_SET**) malloc(num_cols_x*sizeof(SORTED_SET*));
  
  /* have pre-allocated memory for each depth, each direction, each covariate */
  SORTED_SET**** tmp_sorted_sets;

  int* best_actions;
  
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
  tmp_sorted_sets = (SORTED_SET****) malloc((depth+1)*sizeof(SORTED_SET***));
  tmp_sorted_sets[0] = NULL; /* nothing for depth 0 */
  for(i = 1; i < depth+1; i++)
  {
    tmp_sorted_sets[i] = (SORTED_SET***) malloc(2*sizeof(SORTED_SET**));
    tmp_sorted_sets[i][LEFT] = (SORTED_SET**) malloc(num_cols_x*sizeof(SORTED_SET*));
    tmp_sorted_sets[i][RIGHT] = (SORTED_SET**) malloc(num_cols_x*sizeof(SORTED_SET*));

    for( p = 0; p < num_cols_x; p++)
    {
      tmp_sorted_sets[i][LEFT][p] = (SORTED_SET*) malloc(sizeof(SORTED_SET));
      new_sorted_subset(initial_sorted_sets[p],tmp_sorted_sets[i][LEFT][p]);
      tmp_sorted_sets[i][RIGHT][p] = (SORTED_SET*) malloc(sizeof(SORTED_SET));
      new_sorted_subset(initial_sorted_sets[p],tmp_sorted_sets[i][RIGHT][p]);
    }
  }
  
  rewards = (double*) malloc(num_cols_y*sizeof(double));
  rewards2 = (double*) malloc(num_cols_y*sizeof(double));

  best_actions = store_best_actions(data_y, num_rows, num_cols_y);

  find_best_split(tree, depth, initial_sorted_sets, split_step, min_node_size, data_x, data_y, num_rows, num_cols_x, num_cols_y, best_actions,
    tmp_trees, tmp_sorted_sets, rewards, rewards2);  /* these 3 temporary reusable storage */

  free(best_actions);
  
  free(rewards);
  free(rewards2);
    
  /* remove temporary sorted sets */
  for(i = 1; i < depth+1; i++)
  {
     for( p = 0; p < num_cols_x; p++)
     {
        /* don't need to free ->rank since these just shared the pointer */
        free(tmp_sorted_sets[i][LEFT][p]->elements);
        free(tmp_sorted_sets[i][LEFT][p]);
        free(tmp_sorted_sets[i][RIGHT][p]->elements);
        free(tmp_sorted_sets[i][RIGHT][p]);
     }
     free(tmp_sorted_sets[i][LEFT]);
     free(tmp_sorted_sets[i][RIGHT]);
     free(tmp_sorted_sets[i]);
  }
  free(tmp_sorted_sets);
    
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
  {
    free(initial_sorted_sets[p]->elements);
    /* have to free rank */
    free(initial_sorted_sets[p]->rank);
    free(initial_sorted_sets[p]);
  }
  free(initial_sorted_sets);

  return tree;
}

