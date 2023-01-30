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


/* 
 * a sorted set is a single subset of a universal set, we include the ordering of the universal set for each covariate
 * and also breakpoints for each covariate
*/
struct sorted_set
{
  int** sorted_universe;    /** sorted universe set for each covariate */
  int n_universe;           /** size of universe */
  double** data_xx;         /** data_x[p][elt] is the sort value for elt for covariate p */
  char* present;            /** present[elt] = 1 if i in set */
  int n;                    /** size of subset */
  int** breakpoints;        /** array of (ordered) breakpoints for each covariate, 
                                a breakpoint is an index i into su such that
                                su[i-1]<su[i] in underlying ordering */
  int* n_breakpoints;       /** number of breakpoints */
};
typedef struct sorted_set SORTED_SET;



/* debugging only */
static
void print_tree(
  NODE* tree,
  char** covnames
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


/*
 * left run is indices[ileft:iright-1]
 * right run is indices[iright:iend-1]
 */
static
void bottomupmerge(
  int* indices,
  int ileft,
  int iright,
  int iend,
  int* tmp,
  double* data_x_p
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

static
void bottomupmergesort(
  int* indices,
  int* tmp,
  int n,
  double* data_x_p
  )
{
  int width;
  int i;
  size_t size = n*sizeof(int);
  
  for( width = 1; width < n; width *= 2)
  {
    for( i = 0; i < n; i += 2*width)
      bottomupmerge(indices, i, MIN(i+width,n), MIN(i+2*width,n), tmp, data_x_p);

    memcpy(indices,tmp,size);
  }
}

static
double find_split_val(
  SORTED_SET* sorted_set,
  int p,
  int breakpoint
  )
{
  return (sorted_set->data_xx)[p][sorted_set->sorted_universe[p][breakpoint]];
}

/** insert elements into a sorted set, assuming they are not already there */
static
void insert_elements(
  SORTED_SET* sorted_set,
  int* elts,
  int n
  )
{
  int i;
  
  assert(elts != NULL);
  assert(sorted_set != NULL);
  
  for( i = 0; i < n; i++)
    sorted_set->present[elts[i]] = 1;

  sorted_set->n += n;
}

/** remove elements from a sorted set, assuming they are currently there */
static
void remove_elements(
  SORTED_SET* sorted_set,
  int* elts,
  int n
  )
{
  int i;
  
  assert(elts != NULL);
  assert(sorted_set != NULL);
  
  for( i = 0; i < n; i++)
    sorted_set->present[elts[i]] = 0;

  sorted_set->n -= n;
}

static
void next_breakpoint(
  SORTED_SET* sorted_set,
  int p,
  int* indices,
  int* n_indices,
  int current_bki,
  int current_bk,
  int* next_bki,
  int* next_bk
  )
{

  int i;
  int elt;

  assert( sorted_set != NULL );
  assert( indices != NULL );
  assert( n_indices != NULL );
  assert( next_bki != NULL );
  assert( next_bk != NULL );
  assert( p >= 0);
  assert( sorted_set->n_breakpoints[p] > 0);
  assert( current_bki >= -1 );
  assert( current_bk >= 0 );
  
  /* bki is the index of the breakpoint in sorted_set->breakpoints
   * bk is the breakpoint which is an index into sorted_set->sorted_universe
   */
  *n_indices = 0;
  *next_bki = -1;
  *next_bk = -1;
  
  for( *next_bki = current_bki+1; *next_bki < sorted_set->n_breakpoints[p]; (*next_bki)++)
  {
    /* get the next breakpoint */
    assert(*next_bki >= 0);
    *next_bk = sorted_set->breakpoints[p][*next_bki];
    assert(*next_bk >= 0);
    
    /* inspect elements in sorted_universe[current_bk:next_bk-1] */
    for( i = current_bk; i < *next_bk; i++)
    {
      elt = sorted_set->sorted_universe[p][i];
      if( sorted_set->present[elt] )
        indices[(*n_indices)++] = elt;
    }
    if( *n_indices > 0 )
      return;

    /* no elements present in this subset, keep looking .. */
    current_bk = *next_bk;
  }
}

/** make a sorted set the empty set */
static
void make_empty(
  SORTED_SET* sorted_set
  )
{
  memset(sorted_set->present,0,(sorted_set->n_universe)*sizeof(char));
  sorted_set->n = 0;
}

/** copy a sorted set, overwriting an existing target set */
static
void copy_sorted_set(
  SORTED_SET* source,
  SORTED_SET* target
  )
{

  memcpy(target->present,source->present,source->n_universe*sizeof(char));
  target->n = source->n; 
}

/** create a new uninitialised subset from existing sorted subset, universe set is shared */
static
void new_sorted_subset(
  SORTED_SET* source,
  SORTED_SET* target
  )
{
  target->sorted_universe = source->sorted_universe;
  target->n_universe = source->n_universe;
  target->data_xx = source->data_xx;
  /* present is only place where source and target can differ */
  /* don't bother to copy present values */
  target->present = (char* ) malloc(source->n_universe*sizeof(char));
  target->breakpoints = source->breakpoints;
  target->n_breakpoints = source->n_breakpoints;
}


/* make a (multi-)sorted set from elements 0,...,num_indices-1 */
static
SORTED_SET* make_sorted_set(
  int num_indices,            /* size of set */
  int num_cols_x,             /* number of covariates */
  double** data_xx,           /* ordering key values, one for each covariate */
  int* tmp                    /* temporary storage for ordering */
  )
{
  SORTED_SET* sorted_set;
  int n_breakpoints;
  int i;
  int* indices;
  int p;
  
  indices = (int *) malloc(num_indices*sizeof(int));
  sorted_set = (SORTED_SET*) malloc(sizeof(SORTED_SET));
  sorted_set->present = (char *) malloc(num_indices*sizeof(char));

  for( i = 0; i < num_indices; i++)
  {
    indices[i] = i;
    sorted_set->present[i] = 1;
  }
  
  sorted_set->n_universe = num_indices;
  sorted_set->data_xx = data_xx;
  sorted_set->n = num_indices;
  sorted_set->sorted_universe = (int**) malloc(num_cols_x * sizeof(int*));
  sorted_set->breakpoints = (int**) malloc(num_cols_x * sizeof(int*));
  sorted_set->n_breakpoints = (int*) malloc(num_cols_x * sizeof(int));
  for( p = 0; p < num_cols_x; p++ )
  {
    sorted_set->sorted_universe[p] = (int*) malloc(num_indices * sizeof(int));
    sorted_set->breakpoints[p] = (int*) malloc(num_indices * sizeof(int));
    memcpy(sorted_set->sorted_universe[p],indices,num_indices*sizeof(int));
    bottomupmergesort(sorted_set->sorted_universe[p],tmp,num_indices,data_xx[p]);

    n_breakpoints = 0;
    for(i = 1; i < num_indices; i++)
    {
      if( data_xx[p][sorted_set->sorted_universe[p][i-1]] < data_xx[p][sorted_set->sorted_universe[p][i]] )
        sorted_set->breakpoints[p][n_breakpoints++] = i;
    }
    sorted_set->breakpoints[p] = (int*) realloc(sorted_set->breakpoints[p],n_breakpoints*sizeof(int));
    sorted_set->n_breakpoints[p] = n_breakpoints;
  }
  free(indices);
  
  return sorted_set;
}

/*
 * Return a 'skeleton' policy tree of required depth
 * or NULL if there is insufficient memory
 */
static
NODE* make_tree(
  int depth
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
  node->action_id = 0;
  
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
  NODE* source,
  NODE* target
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

/*
 * delete a tree (free the memory it occupied)
 */
void tree_free(
  NODE* node
  )
{
  assert(node != NULL);
  
  if( node->left_child != NULL )
  {
    tree_free(node->left_child);
    tree_free(node->right_child);
  }
  free(node);
}


static
void find_best_reward(
  SORTED_SET* sorted_set,
  double* data_y,
  int num_cols_y,
  int num_rows,
  double* rewards,
  double* best_reward,
  int* best_action
  )
{
  int d;
  int elt;
  
  assert(sorted_set != NULL);
  assert(data_y != NULL);
  assert(best_reward != NULL);
  assert(best_action != NULL);
  assert(num_cols_y >= 0);

  for( d = 0; d < num_cols_y; d++ )
    rewards[d] = 0.0;

  for( elt = 0; elt < sorted_set->n_universe; elt++ )
    if( sorted_set->present[elt] )
    {
      double* dy = data_y;
      for( d = 0; d < num_cols_y; d++ )
      {
        rewards[d] += dy[elt];
        dy += num_rows;
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

/* find optimal depth=1 (ie a single split) tree */
static
void level_one_learning(
  NODE* node,               /** skeleton tree to be populated */
  SORTED_SET* sorted_set,   /** sorted set */
  int split_step,           /** consider splits every split_step'th possible split */
  int min_node_size,        /** smallest terminal node size */
  double* data_x,           /** covariates, data_x[j] are values for covariate j */
  double* data_y,           /** gammas, y[i] are the rewards for unit i */
  int num_rows,             /** number of units */
  int num_cols_x,           /** number of covariates */
  int num_cols_y,           /** number of rewards */
  int* best_actions,        /** best_actions[i] is the best action for unit i */
  int* tmp_indices,                 /** stores indices moved from right branch to left branch */
  double* rewards,                  /** temporary storage for computing best rewards */
  double* rewards2                  /** temporary storage for computing best rewards */
  )
{

  int d;
  double* nosplit_rewards = rewards2;
  double* left_rewards = rewards;

  double left_reward;
  int left_action;
  double best_left_reward;
  int best_left_action = -1;
  double right_reward;
  int right_action;
  double best_right_reward;
  int best_right_action = -1;

  double best_reward;
  int best_action;

  int elt;
  int p;

  int breakpoint_i;
  int breakpoint;
  int i;

  int best_split_var;
  double best_split_val;

  int n_tmp_indices;

  assert(node != NULL);
  assert(node->left_child != NULL);
  assert(node->right_child != NULL);
  assert(data_x != NULL);
  assert(data_y != NULL);
  assert(sorted_set != NULL);
  assert(tmp_indices != NULL);

  
  /* initialise right rewards */
  /* perhaps do one-off compression? */
  for( d = 0; d < num_cols_y; d++ )
      nosplit_rewards[d] = 0.0;
  for( elt = 0; elt < sorted_set->n_universe; elt++ )
    if( sorted_set->present[elt] )
    {
      double* dy = data_y;
      for( d = 0; d < num_cols_y; d++ )
      {
        nosplit_rewards[d] += dy[elt];
        dy += num_rows;
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
  
  /* search for best split */
  for( p = 0; p < num_cols_x; p++)
  {
    /* initialise left rewards for this p */
    for( d = 0; d < num_cols_y; d++ )
      left_rewards[d] = 0.0;
    
    /* find first breakpoint, if any, for covariate p s.t. at least one datapoint moves from left to right */
    next_breakpoint(sorted_set,p,tmp_indices,&n_tmp_indices,-1,0,&breakpoint_i,&breakpoint);
    
    while( n_tmp_indices > 0 )
    {
      /* update left rewards */
      double* dy = data_y;
      if( num_cols_y == 2 )
      {
         double tmp0 = 0.0;
         double tmp1 = 0.0;
         for( i = 0; i < n_tmp_indices; i++ )
         {
            int elt = tmp_indices[i]; 
            tmp0 += data_y[elt];
            tmp1 += data_y[num_rows+elt];
         }
         left_rewards[0] += tmp0;
         left_rewards[1] += tmp1;
      }
      else
      {
         for( d = 0; d < num_cols_y; d++ )
         {
            double tmp = 0.0;
            for( i = 0; i < n_tmp_indices; i++ )
               tmp += dy[tmp_indices[i]];
            left_rewards[d] += tmp;
            dy += num_rows;      
         }
      }
      /* find best left and right reward+action */
      left_reward = left_rewards[0];
      left_action = 0;
      right_reward = nosplit_rewards[0] - left_rewards[0];
      right_action = 0;
      for( d = 1; d < num_cols_y; d++ )
      {
        if( left_rewards[d] > left_reward )
        {
          left_reward = left_rewards[d];
          left_action = d;
        }
        if( nosplit_rewards[d] - left_rewards[d] > right_reward )
        {
          right_reward = nosplit_rewards[d] - left_rewards[d];
          right_action = d;
        }
      }
      
      /* update if new best split */
      if( left_reward + right_reward > best_reward)
      {
        best_reward = left_reward + right_reward;
        best_left_reward = left_reward;
        best_left_action = left_action;
        best_right_reward = right_reward;
        best_right_action = right_action;
        best_split_var = p;
        best_split_val = find_split_val(sorted_set,p,breakpoint);
      }
      
      /* find next breakpoint */
      next_breakpoint(sorted_set,p,tmp_indices,&n_tmp_indices,breakpoint_i,breakpoint,&breakpoint_i,&breakpoint);
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

static
int pure(
  SORTED_SET* sorted_set,
  int* best_actions
  )
{
  int i;
  int best_action = 0;
  int first = 1;
  int n = sorted_set->n;
  
  for( i = 0; i < sorted_set->n_universe; i++)
  {
    if( sorted_set->present[i] )
    {
      if( first )
      {
        best_action = best_actions[i];
        first = 0;
      }
      else
      {
        if( best_actions[i] != best_action )
          return -1;
      }
      n--;
      if( n == 0 )
        return best_action;
    }
  }
  return best_action;
}

/* on return, `node` will be the root of an optimal tree of depth `depth` for the data
 * represented by `sorted_set`
 * arguments `split_step` to `num_cols_y` are constant values needed to determine optimal tree
 * arguments tmp_trees to rewards2 are pre-allocated working space
 */
static
void find_best_split(
  NODE* node,               /** skeleton tree to be populated */
  int depth,                /** depth of tree */
  SORTED_SET* sorted_set,   /** sorted set */
  int split_step,           /** consider splits every split_step'th possible split */
  int min_node_size,        /** smallest terminal node size */
  double* data_x,           /** covariates, data_x[j] are values for covariate j */
  double* data_y,           /** gammas, y[i] are the rewards for unit i */
  int num_rows,             /** number of units */
  int num_cols_x,           /** number of covariates */
  int num_cols_y,           /** number of rewards */
  int* best_actions,        /** best_actions[i] is the best action for unit i */
  NODE*** tmp_trees,                /** trees of various depths for temporary storage */
  int* tmp_indices,                 /** stores indices moved from right branch to left branch */
  SORTED_SET*** tmp_sorted_sets,  /** e.g have tmp_sorted_sets[depth][LEFT] preallocated space */
  double* rewards,                   /** temporary storage for computing best rewards */
  double* rewards2                   /** temporary storage for computing best rewards */
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
  int n_tmp_indices;

  SORTED_SET* left_sorted_set;
  SORTED_SET* right_sorted_set;

  int best_reward_for_all;
  
  assert(node != NULL);
  assert(tmp_trees != NULL);
  assert(data_x != NULL);
  assert(data_y != NULL);
  assert(sorted_set != NULL);
  assert(depth >= 0);
  assert(num_rows > 0);
  assert(num_cols_x > 0);
  assert(num_cols_y > 0);
  

  /* nothing further to do for a leaf or if too few datapoints for splitting */
  if( depth == 0 || sorted_set->n <= min_node_size )
  {
    /* find best reward with no split */
    find_best_reward(sorted_set,data_y,num_cols_y,num_rows,rewards,&best_reward,&best_action);

    /* populate node */
    node->index = -1;
    node->reward = best_reward;
    node->action_id = best_action;
    return;
  }

  /* best_reward_for_all = -1; */
  best_reward_for_all = pure(sorted_set,best_actions);
  if( best_reward_for_all != -1 )
  {
    /* printf("Index set of size %d is pure with all units having %d as best action\n", sorted_set->n, best_reward_for_all); */

    /* find best reward with no split */
    find_best_reward(sorted_set,data_y,num_cols_y,num_rows,rewards,&best_reward,&best_action);

    assert(best_reward == best_reward_for_all);
    
    /* populate node */
    node->index = -1;
    node->reward = best_reward;
    node->action_id = best_action;
    return;
  }
  
  if( depth == 1 )
  {
    level_one_learning(node, sorted_set, split_step, min_node_size, data_x, data_y, num_rows, num_cols_x, num_cols_y,  best_actions, tmp_indices,
      rewards, rewards2);
    return;
  }

  left_child = node->left_child;
  right_child = node->right_child;
  best_left_child = tmp_trees[depth-1][LEFT];
  best_right_child = tmp_trees[depth-1][RIGHT];
  left_sorted_set = tmp_sorted_sets[depth][LEFT];
  right_sorted_set = tmp_sorted_sets[depth][RIGHT];

  /* consider each covariate for splitting */
  for( p = 0; p < num_cols_x; p++)
  {
    /* reinitialise left and right sorted sets */
    make_empty(left_sorted_set);
    copy_sorted_set(sorted_set,right_sorted_set);

    /* find first breakpoint, if any, for covariate p s.t. at least one datapoint moves from left to right */
    next_breakpoint(sorted_set,p,tmp_indices,&n_tmp_indices,-1,0,&breakpoint_i,&breakpoint);
    
    while( n_tmp_indices > 0 )
    {
      /* move points from right to left  */
      insert_elements(left_sorted_set,tmp_indices,n_tmp_indices);
      remove_elements(right_sorted_set,tmp_indices,n_tmp_indices);

      find_best_split(left_child, depth-1, left_sorted_set, split_step, min_node_size, data_x, data_y, num_rows, num_cols_x, num_cols_y, best_actions,
         tmp_trees, tmp_indices, tmp_sorted_sets, rewards, rewards2);
      find_best_split(right_child, depth-1, right_sorted_set, split_step, min_node_size, data_x, data_y, num_rows, num_cols_x, num_cols_y, best_actions,
         tmp_trees, tmp_indices, tmp_sorted_sets, rewards, rewards2);

      reward = left_child->reward + right_child->reward;
        
      if( reward > best_reward )
      {
        best_reward = reward;
        best_split_var = p;
        best_split_val = find_split_val(sorted_set,p,breakpoint);
        /* tree_copy(source,target) */
        tree_copy(left_child,best_left_child);
        tree_copy(right_child,best_right_child);
      }

      /* find next breakpoint */
      next_breakpoint(sorted_set,p,tmp_indices,&n_tmp_indices,breakpoint_i,breakpoint,&breakpoint_i,&breakpoint);
    }
  }

  /* populate node */
  node->index = best_split_var;
  node->value = best_split_val;
  node->reward = best_reward;
  if( best_split_var != -1 )
  {
    /* tree_copy(source,target) */
    tree_copy(best_left_child,left_child);
    tree_copy(best_right_child,right_child);
  }
}

/*
 * For each datapoint (index) find and record best action for that datapoint
 * Return pointer to the array of best actions
 */
static
int* store_best_actions(
  double* data_y,     /** gammas (column major) */
  int num_rows,       /** number of units */
  int num_cols_y      /** number of rewards */
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

/*
 * Notes: A tree of depth 0 is a single node (with no children) where only the reward and action_id members are meaningful.
 */

NODE* tree_search_jc_discretedata(
  int depth,            /** (maximum) depth of returned tree */
  int split_step,       /** consider splits every split_step'th possible split */
  int min_node_size,    /** smallest terminal node size */
  double* data_x, /** covariates (column major) */
  double* data_y, /** gammas (column major) */
  int num_rows,         /** number of units */
  int num_cols_x,       /** number of covariates */
  int num_cols_y        /** number of rewards */
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

  SORTED_SET* sorted_set;
  SORTED_SET*** tmp_sorted_sets;

  double** data_xx;

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
  /* a bit neater to store pointer to each column */
  data_xx = (double**) malloc(num_cols_x*sizeof(double*));
  for( p = 0; p < num_cols_x; p++)
    data_xx[p] = data_x+p*num_rows;

  /* make initial (multi-) sorted set */
  sorted_set = make_sorted_set(num_rows,num_cols_x,data_xx,tmp_indices);

  /* make temporary trees (for storing 'best so far' trees) */
  tmp_trees = (NODE***) malloc(depth*sizeof(NODE**));
  for(i = 0; i < depth; i++)
  {
    tmp_trees[i] = (NODE**) malloc(2*sizeof(NODE*));
    tmp_trees[i][LEFT] = make_tree(i);
    tmp_trees[i][RIGHT] = make_tree(i);
  }

  /* make temporary sorted sets */
  tmp_sorted_sets = (SORTED_SET***) malloc((depth+1)*sizeof(SORTED_SET**));
  tmp_sorted_sets[0] = NULL; /* nothing for depth 0 */
  for(i = 1; i < depth+1; i++)
  {
    tmp_sorted_sets[i] = (SORTED_SET**) malloc(2*sizeof(SORTED_SET*));
    tmp_sorted_sets[i][LEFT] = (SORTED_SET*) malloc(sizeof(SORTED_SET));
    tmp_sorted_sets[i][RIGHT] = (SORTED_SET*) malloc(sizeof(SORTED_SET));

    new_sorted_subset(sorted_set,tmp_sorted_sets[i][LEFT]);
    new_sorted_subset(sorted_set,tmp_sorted_sets[i][RIGHT]);
  }
  
  rewards = (double*) malloc(num_cols_y*sizeof(double));
  rewards2 = (double*) malloc(num_cols_y*sizeof(double));

  best_actions = store_best_actions(data_y, num_rows, num_cols_y);
  
  find_best_split(tree, depth, sorted_set, split_step, min_node_size, data_x, data_y, num_rows, num_cols_x, num_cols_y, best_actions,
    tmp_trees, tmp_indices, tmp_sorted_sets, rewards, rewards2);  /* these 3 temporary reusable storage */

  free(best_actions);
  
  free(tmp_indices);
  free(data_xx);
  free(rewards);
  free(rewards2);
    
  /* remove temporary sorted sets */
  for(i = 1; i < depth+1; i++)
  {
    free(tmp_sorted_sets[i][LEFT]->present);
    free(tmp_sorted_sets[i][RIGHT]->present);
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

  /* free sorted set */
  for( p = 0; p < num_cols_x; p++)
  {
    free(sorted_set->sorted_universe[p]);
    free(sorted_set->breakpoints[p]);
  }
  free(sorted_set->sorted_universe);
  free(sorted_set->breakpoints);
  free(sorted_set->n_breakpoints);
  free(sorted_set->present);
  free(sorted_set);

  return tree;
}

