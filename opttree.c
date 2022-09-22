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
 * sorted set, policy tree style
*/
struct sorted_set
{
  int* elements;    /** the sorted set */
  int n;            /** size of subset */
  int* rank;        /** rank[elt] is the (unique) rank value for elt  */
};
typedef struct sorted_set SORTED_SET;


static
int print_set(
  SORTED_SET* sorted_set
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


/** remove element into a sorted set, assuming it is already there 
 */
static
int remove_element(
  SORTED_SET* sorted_set,
  int elt
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

/** insert element into a sorted set, assuming it not already there */
static
void insert_element(
  SORTED_SET* sorted_set,
  int elt
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

/** make a sorted set the empty set */
static
void make_empty(
  SORTED_SET* sorted_set
  )
{
  sorted_set->n = 0;
}

/** copy a sorted set, overwriting an existing target set */
static
void copy_sorted_set(
  SORTED_SET* source,
  SORTED_SET* target
  )
{

  memcpy(target->elements,source->elements,source->n*sizeof(int));
  target->n = source->n; 
}

/** create a new uninitialised subset from existing sorted subset */
static
void new_sorted_subset(
  SORTED_SET* source,
  SORTED_SET* target
  )
{
  size_t size = (source->n)*sizeof(int);
  
  target->elements = (int*) malloc(size);
  target->rank = (int*) malloc(size);
  memcpy(target->rank,source->rank,size);
}


/* make a sorted set from elements 0,...,num_indices-1 */
static
SORTED_SET* make_sorted_set(
  int num_indices,            /* size of set */
  double* data_xx,            /* ordering values */
  int* tmp                    /* temporary storage for ordering */
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
  double* dyelt;
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

/* find optimal depth=1 (ie a single split) tree */
static
void level_one_learning(
  NODE* node,                /** skeleton tree to be populated */
  SORTED_SET** sorted_sets,  /** sorted sets */
  int split_step,            /** consider splits every split_step'th possible split */
  double* data_x,            /** covariates, data_x[j] are values for covariate j */
  double* data_y,            /** gammas, y[i] are the rewards for unit i */
  int num_rows,              /** number of units */
  int num_cols_x,            /** number of covariates */
  int num_cols_y,            /** number of rewards */
  double* rewards,           /** temporary storage for computing best rewards */
  double** rewards2          /** temporary storage for computing best rewards */
  )
{

  int d;
  double* nosplit_rewards = rewards;
  double** left_rewards = rewards2;

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
  double* dyelt;

  int best_split_var;
  double best_split_val;


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

  /* printf("No split best reward is %g, best action is %d\n", best_reward, best_action); */
  
  /* search for best split */
  for( p = 0; p < num_cols_x; p++)
  {
    SORTED_SET* sorted_set = sorted_sets[p];
    double* data_xp = data_x+(p*num_rows);
    int* elements = (sorted_sets[p])->elements;
    
    /* left_rewards[d][i] is reward for action d if all elements up to and including
     * element index i are on the left */

    /* initialise left_rewards */
    elt = elements[0];
    for( d = 0; d < num_cols_y; d++ )
      left_rewards[d][0] = data_y[num_rows*d+elt];

    /* don't move last element from right to left since then we would have no split */
    for( i = 1; i < (sorted_set->n)-1; i++ )
    {
      elt = elements[i];
      for( d = 0; d < num_cols_y; d++ )
        left_rewards[d][i] = left_rewards[d][i-1] + data_y[num_rows*d+elt];
    }

    /* now find best split for covariate p */
    for( i = 0; i < (sorted_set->n)-1; i++ )
    {
      assert( data_xp[elements[i]] <= data_xp[elements[i+1]] );

      /* if a proper split see whether it's a best split */
      if( data_xp[elements[i]] < data_xp[elements[i+1]] )
      {
        /* find best left and right reward+action */
        best_left_reward_for_split = left_rewards[0][i];
        best_left_action_for_split = 0;
        best_right_reward_for_split = nosplit_rewards[0] - left_rewards[0][i];
        best_right_action_for_split = 0;
        for( d = 1; d < num_cols_y; d++ )
        {
          if( left_rewards[d][i] > best_left_reward_for_split )
          {
            best_left_reward_for_split = left_rewards[d][i];
            best_left_action_for_split = d;
          }
          if( nosplit_rewards[d] - left_rewards[d][i] > best_right_reward_for_split )
          {
            best_right_reward_for_split = nosplit_rewards[d] - left_rewards[d][i];
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
          best_split_val = data_xp[elements[i]];

          /* printf("New best reward is %g=%g+%g, best actions are %d,%d split var is %d, split val is %g\n", best_reward, */
          /*   best_left_reward, best_right_reward, best_left_action, best_right_action, p, best_split_val); */
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

/* on return, `node` will be the root of an optimal tree of depth `depth` for the data
 * represented by `sorted_sets`
 * sorted_sets[0], sorted_sets[1], ... are the same set, they only differ in the order of the elements
 * sorted_sets[p] has its elements sorted according to their rank which is a strict total order
 * consistent with order given by the pth covariate
 * arguments `split_step` to `num_cols_y` are constant values needed to determine optimal tree
 * arguments tmp_trees to rewards2 are pre-allocated working space
 */
static
void find_best_split(
  NODE* node,               /** skeleton tree to be populated */
  int depth,                /** depth of tree */
  SORTED_SET** sorted_sets, /** sorted sets, one for each covariate */
  int split_step,           /** consider splits every split_step'th possible split */
  int min_node_size,        /** smallest terminal node size */
  double* data_x,           /** covariates, data_x[j] are values for covariate j */
  double* data_y,           /** gammas, y[i] are the rewards for unit i */
  int num_rows,             /** number of units */
  int num_cols_x,           /** number of covariates */
  int num_cols_y,           /** number of rewards */
  NODE*** tmp_trees,                /** trees of various depths for temporary storage */
  SORTED_SET**** tmp_sorted_sets,   /** e.g have tmp_sorted_sets[depth][LEFT][p] preallocated space */
  double* rewards,                  /** temporary storage for computing best rewards */
  double** rewards2                  /** temporary storage for computing best rewards */
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
  
  assert(node != NULL);
  assert(tmp_trees != NULL);
  assert(data_x != NULL);
  assert(data_y != NULL);
  assert(sorted_sets != NULL);
  assert(depth >= 0);
  assert(num_rows > 0);
  assert(num_cols_x > 0);
  assert(num_cols_y > 0);
  

  /* nothing further to do for a leaf or if too few datapoints for splitting */
  if( depth == 0 || sorted_sets[0]->n <= min_node_size )
  {
    /* find best reward with no split */
    find_best_reward(sorted_sets[0],data_y,num_cols_y,num_rows,rewards,&best_reward,&best_action);

    /* populate node */
    node->index = -1;
    node->reward = best_reward;
    node->action_id = best_action;
    return;
  }

  if( depth == 1 )
  {
    level_one_learning(node, sorted_sets, split_step, data_x, data_y, num_rows, num_cols_x, num_cols_y,
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
    SORTED_SET* sorted_set = sorted_sets[p];
    double* data_xp = data_x+(p*num_rows);

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
      
      assert( data_xp[elt] <= data_xp[sorted_set->elements[i+1]] );

      /* if a proper split see whether it's a best split */
      if( data_xp[elt] < data_xp[sorted_set->elements[i+1]] )
      {
        find_best_split(left_child, depth-1, tmp_sorted_sets[depth][LEFT], split_step, min_node_size, data_x, data_y, num_rows, num_cols_x, num_cols_y,
          tmp_trees, tmp_sorted_sets, rewards, rewards2);
        find_best_split(right_child, depth-1, tmp_sorted_sets[depth][RIGHT], split_step, min_node_size, data_x, data_y, num_rows, num_cols_x, num_cols_y,
          tmp_trees, tmp_sorted_sets, rewards, rewards2);

        reward = left_child->reward + right_child->reward;
        
        if( reward > best_reward )
        {
          best_reward = reward;
          best_split_var = p;
          best_split_val = data_xp[elt];
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

NODE* tree_search(
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
  int d;
  NODE*** tmp_trees;
  int* tmp_indices;
  double* rewards;
  double** rewards2;

  SORTED_SET** initial_sorted_sets;
  /* have pre-allocated memory for each depth, each direction, each covariate */
  SORTED_SET**** tmp_sorted_sets;

  /* data_x[col * num_rows + row] is the value of covariate 'col' in datapoint
   * 'row', so column major storage, all values for covariate 0 first
   *
   * data_y[col * num_rows + row] is the value of reward 'col' in datapoint
   * 'row', so column major storage, all values of reward 0 first.
   */
  
  /* make skeleton tree of correct depth */
  tree = make_tree(depth);

  tmp_indices = (int*) malloc(num_rows*sizeof(int));
  initial_sorted_sets = (SORTED_SET**) malloc(num_cols_x*sizeof(SORTED_SET*));
  
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
  rewards2 = (double**) malloc(num_cols_y*sizeof(double*));
  for( d = 0; d < num_cols_y; d++ )
    rewards2[d] = (double*) malloc(num_rows*sizeof(double));
  
  find_best_split(tree, depth, initial_sorted_sets, split_step, min_node_size, data_x, data_y, num_rows, num_cols_x, num_cols_y,
    tmp_trees, tmp_sorted_sets, rewards, rewards2);  /* these 3 temporary reusable storage */

  free(rewards);
  for( d = 0; d < num_cols_y; d++ )
    free(rewards2[d]);
  free(rewards2);
  
  /* remove temporary sorted sets */
  for(i = 1; i < depth+1; i++)
  {
    for( p = 0; p < num_cols_x; p++)
    {
      free(tmp_sorted_sets[i][LEFT][p]->elements);
      free(tmp_sorted_sets[i][LEFT][p]->rank);
      free(tmp_sorted_sets[i][LEFT][p]);
      free(tmp_sorted_sets[i][RIGHT][p]->elements);
      free(tmp_sorted_sets[i][RIGHT][p]->rank);
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
    free(initial_sorted_sets[p]->rank);
    free(initial_sorted_sets[p]);
  }
  free(initial_sorted_sets);

  return tree;
}

