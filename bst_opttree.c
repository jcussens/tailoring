#include <float.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include "bst_opttree.h"
#include "bst.h"

#define INF DBL_MAX
#define LEFT 0
#define RIGHT 1
#define MIN(a,b) (((a)<(b))?(a):(b))
#define VERBOSE 0
#define VERY_VERBOSE 0


static
int pure(
   BST* sorted_set,
   int* best_actions
  )
{
  int elt;
  int best_action = 0;
  BST* bstnode;
  
  bstnode = find_minimum(sorted_set,&elt);

  if( bstnode != NULL )
  {
     /* set best action to that for first element found */
     best_action = best_actions[elt];

     for( bstnode = find_next(bstnode,&elt); bstnode != NULL; bstnode = find_next(bstnode,&elt) )
     {
        if( best_action != best_actions[elt] )
           return -1;
     } 
  }

  return best_action;
}

/*
 * For each datapoint (index) find and record best action for that datapoint
 * Return pointer to the array of best actions
 */
static
int* store_best_actions(
   double* data_y,     /** gammas (column major) */
  int num_rows,             /** number of units */
  int num_cols_y            /** number of rewards */
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

/* make a sorted set from elements 0,...,num_indices-1 */
static
BST* make_sorted_set(
   int num_indices,        /* size of set */
   double* data_xx,        /* ordering values */
  int* rank,                    /* rank[elt] is rank of elt on exit */
  int* tmp                      /* temporary storage for ordering */
  )
{
  BST* sorted_set;
  int i;
  int* elements;

  elements = (int *) malloc(num_indices*sizeof(int));
  for( i = 0; i < num_indices; i++)
    elements[i] = i;

  bottomupmergesort(elements,tmp,num_indices,data_xx);

  for( i = 0; i < num_indices; i++)
     rank[elements[i]] = i;
  
  sorted_set = make_bst(elements,num_indices,0,NULL);

  free(elements);

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
   BST* sorted_set,
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
  BST* bstnode;
  int elt;
  
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

  for( bstnode = find_minimum(sorted_set,&elt); bstnode != NULL; bstnode = find_next(bstnode,&elt) )
  {
     dyelt = data_y + elt;
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
  NODE* node,                      /** skeleton tree to be populated */
  BST** sorted_sets,  /** sorted sets */
  int n,                            /** number of datapoints reaching this node */
   int split_step,                  /** consider splits every split_step'th possible split */
   double* data_x,                  /** covariates, data_x[j] are values for covariate j */
   double* data_y,                  /** gammas, y[i] are the rewards for unit i */
   int num_rows,                    /** number of units */
   int num_cols_x,                  /** number of covariates */
   int num_cols_y,                  /** number of rewards */
   int* best_actions,         /** best_actions[i] is the best action for unit i */
  double* rewards,                 /** temporary storage for computing best rewards */
  double* rewards2                 /** temporary storage for computing best rewards */
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
  int eltnxt;
  int p;

  int i;
   double* dyelt;

  int best_split_var;
  double best_split_val;

  BST* bstnode;

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
  for( bstnode = find_minimum(sorted_sets[0],&elt); bstnode != NULL; bstnode = find_next(bstnode,&elt) )
  {
    dyelt = data_y + elt;
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
     printf("No split best reward for %d datapoints is %g, best action is %d\n", n, best_reward, best_action); 
  
  /* search for best split */
  for( p = 0; p < num_cols_x; p++)
  {
     BST* sorted_set = sorted_sets[p];
     double* data_xp = data_x+(p*num_rows);

     n_left = 0;
     n_right = n;
     
    /* initialise left rewards for this p */
    for( d = 0; d < num_cols_y; d++ )
      left_rewards[d] = 0.0;

    if( VERY_VERBOSE )
    {
      in_order_tree_walk(sorted_set);
      printf("\n");
      in_order_tree_walk_fast(sorted_set);
      printf("\n");

    }
    
    /* find lowest ranked element for covariate p */
    bstnode = find_minimum(sorted_set,&elt);

    while( bstnode != NULL )
    {

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
      
      /* find next element to see whether we have a valid split */
      bstnode = find_next(bstnode,&eltnxt);

      /* if no next element or next element has same value for covariate p, don't consider split */
      if( bstnode != NULL && data_xp[elt] != data_xp[eltnxt] )
      {
        if( VERY_VERBOSE )
          printf("valid split at %d, since %d has different %d value\n",elt,eltnxt,p); 

        
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
                  n,n_left,n_right,p,best_split_val,best_reward,best_left_reward,best_right_reward);
        }
      }
      elt = eltnxt; 
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
  NODE* node,                     /** skeleton tree to be populated */
   int depth,                /** depth of tree */
  BST** sorted_sets,              /** sorted sets, one for each covariate */
  int** ranks,              /** ranks[p][elt] is the rank for element elt for covariate p */
  int n,                          /** number of datapoints reaching this node (ie common size of sorted sets) */ 
   int split_step,           /** consider splits every split_step'th possible split */
   int min_node_size,        /** smallest terminal node size */
   double* data_x,           /** covariates, data_x[j] are values for covariate j */
   double* data_y,           /** gammas, y[i] are the rewards for unit i */
   int num_rows,             /** number of units (in full dataset) */
   int num_cols_x,           /** number of covariates */
   int num_cols_y,           /** number of rewards */
   int* best_actions,        /** best_actions[i] is the best action for unit i */
  NODE*** tmp_trees,              /** trees of various depths for temporary storage */
  BST**** tmp_sorted_sets,        /** e.g have tmp_sorted_sets[depth][LEFT][p] preallocated space */
  double* rewards,                /** temporary storage for computing best rewards */
  double* rewards2                /** temporary storage for computing best rewards */
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

  BST* left_sorted_set;
  BST* right_sorted_set;

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
  if( depth == 0 || n <= min_node_size || pure(sorted_sets[0],best_actions) != -1)
  {
     
    /* find best reward with no split */
    find_best_reward(sorted_sets[0],data_y,num_cols_y,num_rows,rewards,&best_reward,&best_action);

    if( VERBOSE && pure(sorted_sets[0],best_actions) != -1)
        printf("Node with %d datapoints for depth %d tree is pure with best action %d and reward=%g\n",
           n,depth,best_action,best_reward);

    
    /* populate node */
    node->index = -1;
    node->reward = best_reward;
    node->action_id = best_action;
    return;
  }

  if( depth == 1 )
  {
     level_one_learning(node, sorted_sets, n, split_step, data_x, data_y, num_rows, num_cols_x, num_cols_y, best_actions,
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
     int n_right = n;
     
     BST* sorted_set = sorted_sets[p];
      double* data_xp = data_x+(p*num_rows);
     
     BST* bstnode;
     int elt;
     int eltnxt;
    
     /* initialise all sorted sets to be empty on left, and full on right */
     for( pp = 0; pp < num_cols_x; pp++)
     {
        empty_data_bst(sorted_sets[pp],tmp_sorted_sets[depth][LEFT][pp]);
        copy_data_bst(sorted_sets[pp],tmp_sorted_sets[depth][RIGHT][pp]);
     }

    /* find lowest ranked element for covariate p */
    bstnode = find_minimum(sorted_set,&elt);

    while( bstnode != NULL )
    {
       /* move element from right to left branch */
       for( pp = 0; pp < num_cols_x; pp++)
       {
         printf("foo %d %d\n",pp,elt);
          insert(tmp_sorted_sets[depth][LEFT][pp],elt,ranks[pp][elt]);
          delete(tmp_sorted_sets[depth][RIGHT][pp],elt,ranks[pp][elt]);
       }
       n_left++;
       n_right--;

       /* find next element to see whether we have a valid split */
       bstnode = find_next(bstnode,&eltnxt);

      /* if no next element or next element has same value for covariate p, don't consider split */
      if( bstnode != NULL && data_xp[elt] != data_xp[eltnxt] )
      {
       
         find_best_split(left_child, depth-1, tmp_sorted_sets[depth][LEFT], ranks, n_left, split_step, min_node_size,
            data_x, data_y, num_rows, num_cols_x, num_cols_y, best_actions,
            tmp_trees, tmp_sorted_sets, rewards, rewards2);
         find_best_split(right_child, depth-1, tmp_sorted_sets[depth][RIGHT], ranks, n_right, split_step, min_node_size,
            data_x, data_y, num_rows, num_cols_x, num_cols_y, best_actions,
            tmp_trees, tmp_sorted_sets, rewards, rewards2);
      
         reward = left_child->reward + right_child->reward;
        
         if( reward > best_reward )
         {
            best_reward = reward;
            best_split_var = p;
            best_split_val = data_xp[elt];
            if( VERBOSE )
               printf("New best split for %d=%d+%d datapoints for depth %d tree is: covariate=%d, split value=%g, reward=%g=%g+%g\n",
                  n,n_left,n_right,depth,p,best_split_val,reward,left_child->reward,right_child->reward);
            /* save best left and right trees found so far */
            /* tree_copy(source,target) */
            tree_copy(left_child,best_left_child);
            tree_copy(right_child,best_right_child);
         }
      }
      elt = eltnxt;
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
   double* data_x,       /** covariates (column major) */
   double* data_y,       /** gammas (column major) */
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

  BST** initial_sorted_sets = (BST**) malloc(num_cols_x*sizeof(BST*));
  
  /* have pre-allocated memory for each depth, each direction, each covariate */
  BST**** tmp_sorted_sets;

  int* best_actions;

  int** ranks;
  
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
  ranks = (int **) malloc(num_cols_x*sizeof(int*));
  for( p = 0; p < num_cols_x; p++)
  {
     ranks[p] = (int *) malloc(num_rows*sizeof(int));
     initial_sorted_sets[p] = make_sorted_set(num_rows,data_x+p*num_rows,ranks[p],tmp_indices);
  }
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
  tmp_sorted_sets = (BST****) malloc((depth+1)*sizeof(BST***));
  tmp_sorted_sets[0] = NULL; /* nothing for depth 0 */
  for(i = 1; i < depth+1; i++)
  {
    tmp_sorted_sets[i] = (BST***) malloc(2*sizeof(BST**));
    tmp_sorted_sets[i][LEFT] = (BST**) malloc(num_cols_x*sizeof(BST*));
    tmp_sorted_sets[i][RIGHT] = (BST**) malloc(num_cols_x*sizeof(BST*));

    for( p = 0; p < num_cols_x; p++)
    {
       tmp_sorted_sets[i][LEFT][p] = copy_bst(initial_sorted_sets[p],NULL);
       tmp_sorted_sets[i][RIGHT][p] = copy_bst(initial_sorted_sets[p],NULL);
    }
  }
  
  rewards = (double*) malloc(num_cols_y*sizeof(double));
  rewards2 = (double*) malloc(num_cols_y*sizeof(double));

  best_actions = store_best_actions(data_y, num_rows, num_cols_y);

  find_best_split(tree, depth, initial_sorted_sets, ranks, num_rows, split_step, min_node_size, data_x, data_y, num_rows, num_cols_x, num_cols_y, best_actions,
    tmp_trees, tmp_sorted_sets, rewards, rewards2);  /* these 3 temporary reusable storage */

  free(best_actions);
  
  free(rewards);
  free(rewards2);
    
  /* remove temporary sorted sets */
  for(i = 1; i < depth+1; i++)
  {
    for( p = 0; p < num_cols_x; p++)
    {
       delete_bst(tmp_sorted_sets[i][LEFT][p]);
       delete_bst(tmp_sorted_sets[i][RIGHT][p]);
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
     delete_bst(initial_sorted_sets[p]);
  free(initial_sorted_sets);

  for( p = 0; p < num_cols_x; p++)
     free(ranks[p]);
  free(ranks);
  
  return tree;
}

