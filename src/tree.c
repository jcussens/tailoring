/** @file tree.c
 *  @brief Functions for policy trees
 *  @author James Cussens
 */
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <stddef.h>
#include "tree.h"

/** node of a policytree (for a particular dataset of units)
 *  the associated dataset of units is not represented, only the 
 *  reward for the associated dataset for the tree with the node as root
 */
struct node
{
   int                   index;              /**< splitting covariate */      
   double                value;              /**< split value */
   double                reward;             /**< reward for tree with this node on an unrepresented associated dataset of units */
   int                   action_id;          /**< best action (only meaningful for leaves) */
   struct node*          left_child;         /**< left child (or NULL) */
   struct node*          right_child;        /**< right child (or NULL) */
};

/* LOCAL FUNCTIONS */

/** prints a policy tree in policytree style to standard output */
static
void print_tree_policytree_rec(
   NODE*                 tree,               /**< policy tree to print */
   char**                covnames,           /**< covnames[j] is the name of covariate j */
   int                   level               /**< level of the tree */
  )
{

   /* indent */
   int i;

   for(i = 0; i < level; i++)
      printf("  ");
   
   if( tree->index != -1)
   {
      printf("split_variable: %s  split_value: %g\n", covnames[tree->index], tree->value);
      print_tree_policytree_rec(tree->left_child, covnames, level+1);
      print_tree_policytree_rec(tree->right_child, covnames, level+1);
   }
   else
   {
      printf("* action: %d\n", tree->action_id + 1);
   }
}

/** does a node split data using some value of some variable? */
#ifndef NDEBUG
static
int is_varsplit(
   NODE*                 node                /**< node */
   )
{
   assert(node != NULL);
   
   return (node->index != -1);
}
#endif

/** as soon as we hit a leaf, remove any subtrees underneath it
*/
static
void prune_tree(
   NODE*                 node                /**< root node */
   )
{

   assert(node != NULL);

   if( is_leaf(node) )
   {
      if( node->left_child != NULL )
      {
         tree_free(node->left_child);
         node->left_child = NULL;
      }
      if( node->right_child != NULL )
      {
         tree_free(node->right_child);
         node->right_child = NULL;
      }
   }
   else
   {
      assert(is_varsplit(node));
      assert(node->left_child != NULL);
      assert(node->right_child != NULL);
      prune_tree(node->left_child);
      prune_tree(node->right_child);
   }
}

/** if a node has two child leaves with the same action remove those children
 * and make node a leaf with that action
 */
static
void merge_leaves(
   NODE*                 node                /**< node */
   )
{
   assert(node != NULL);
   
   if( node->left_child != NULL )
      merge_leaves(node->left_child);

   if( node->right_child != NULL )
      merge_leaves(node->right_child);

   if(
      node->left_child != NULL && is_leaf(node->left_child) &&
      node->right_child != NULL && is_leaf(node->right_child) &&
      node->left_child->action_id == node->right_child->action_id )
   {
      /* node should not currently be a leaf */
      assert( !is_leaf(node) );

      /* this function should be called after dummy nodes removed, so leaves
         should have no children, so we just free them after using their values */
      assert(node->left_child->left_child == NULL);
      assert(node->left_child->right_child == NULL);
      assert(node->right_child->left_child == NULL);
      assert(node->right_child->right_child == NULL);
      
      /* make a new leaf */
      node->index = -1;
      node->value = 0.0;
      node->reward = node->left_child->reward + node->right_child->reward;
      node->action_id = node->left_child->action_id;
      free(node->left_child);
      free(node->right_child);
      node->left_child = NULL;
      node->right_child = NULL;
   }
}



/* PUBLIC FUNCTIONS */

/** print a policy tree (debugging only) 
 * if covariate names are not supplied then the indices for covariates are used
 */
void print_tree(
   const NODE*           tree,               /**< root of tree to print */ 
   const char**          covnames            /**< if covnames != NULL, then covnames[i] is the name of covariate i */
  )
{
   assert(tree != NULL);
   assert(tree->index == -1 || tree->left_child != NULL );
   assert(tree->index == -1 || tree->right_child != NULL );

   printf("node = %p\n", (void*) tree);
   printf("reward = %g\n", tree->reward);
   if( tree->index != -1)
   {
      if( covnames != NULL )
         printf("covariate = %s\n", covnames[tree->index]);
      else
         printf("covariate = %d\n", tree->index);
      printf("value = %g\n", tree->value);
      printf("left_child = %p\n", (void*) tree->left_child);
      printf("right_child = %p\n", (void*) tree->right_child);
   }
   else
   {
      printf("action_id = %d\n", tree->action_id);
   }
   printf("\n");
   
   if( tree->index != -1)
   {
      print_tree(tree->left_child,covnames);
      print_tree(tree->right_child,covnames);
   }
}

/** prints a policy tree in policytree style to standard output */
void print_tree_policytree(
   NODE*                 tree,               /**< policy tree to print */
   char**                covnames,           /**< covnames[j] is the name of covariate j */
   int                   depth,              /**< depth of the tree */
   int                   nactions,           /**< number of actions */
   char**                actionnames         /**< actionnames[j] is the name of action j */
  )
{
   int i;
   
   printf("policy_tree object\n");
   printf("Tree depth:   %d\n", depth);
   printf("Actions: ");
   for(i = 0; i < nactions; i++)
      printf(" %d: %s", i+1, actionnames[i]);
   printf("\nVariable splits:\n");

   print_tree_policytree_rec(tree, covnames, 0);   
}


/** is a node a leaf node? 
 * @return 1 if node is a leaf, else 0
*/
int is_leaf(
   NODE*                 node                /**< node */
   )
{
   assert(node != NULL);
   
   return (node->index == -1 && node->action_id != -1);
}


/** prints a policy tree to standard output */
void print_tree_raw(
   NODE*                 tree               /**< policy tree to print */
  )
{
   printf("<<<<<<<\n");
   printf("node = %p\n", (void*) tree);
   printf("index = %d\n", tree->index);
   printf("value = %g\n", tree->value);
   printf("reward = %g\n", tree->reward);
   printf("action_id = %d\n", tree->action_id);
   printf("left_child = %p\n", (void*) tree->left_child);
   printf("right_child = %p\n", (void*) tree->right_child);
   printf(">>>>>>>\n\n");
   
   if( tree->left_child != NULL )
      print_tree_raw(tree->left_child);
   
   if( tree->right_child != NULL )
      print_tree_raw(tree->right_child);
}

/**
 * Return a 'skeleton' policy tree of required depth
 * or NULL if there is insufficient memory
 */
NODE* make_tree(
   int                    depth               /**< depth of required tree */
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

  node->index = -1;     /* no splitting covariate so far, may never be one... */
  node->value = 0;
  node->reward = 0;
  node->action_id = -1; /* no associated action so far, may never be one... */
  node->left_child = NULL;
  node->right_child = NULL;
  
  if( depth > 0 )
  {
    node->left_child = make_tree(depth-1);
    node->right_child = make_tree(depth-1);
    if( node->left_child == NULL || node->right_child == NULL )
      /* not enough memory */
      return NULL;
  }

  return node;
}

/**
 * copy data from source tree to existing target tree
 * assumes both trees have same shape.
 * NB. No new memory is allocated
 */
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

/** add a given reward to those nodes in a tree which a particular element visits */
void update_rewards(
   NODE*                 tree,               /**< tree */
   const double*         data_x,             /**< covariates, data_x+(j*num_rows) points to values for covariate j */
   int                   num_rows,           /**< number of units */
   int                   elt,                /**< element */
   double                reward              /**< reward */
   )
{
   tree->reward += reward;

   if( tree->index == -1 )
      return;
   
   if( data_x[(tree->index)*num_rows+elt] <= tree->value )
      update_rewards(tree->left_child, data_x, num_rows, elt, reward);
   else
      update_rewards(tree->right_child, data_x, num_rows, elt, reward);
}

/** find the action assigned to a unit by a tree 
 * @return the assigned action
 */
int assigned_action(
   const NODE*           tree,               /**< tree */
   const double*         data_x,             /**< covariates, data_x+(j*num_rows) points to values for covariate j */
   int                   num_rows,           /**< number of units */
   int                   elt                 /**< element to find assigned action for */
   )
{
   assert( tree != NULL );
   assert( data_x != NULL );
   /* if node not a leaf, then both children are present */
   assert( tree->index == -1 || tree->left_child != NULL );
   assert( tree->index == -1 || tree->right_child != NULL );
   
   if( tree->index == -1 )
      return tree->action_id;
   else if( data_x[(tree->index)*num_rows+elt] <= tree->value )
      return assigned_action(tree->left_child, data_x, num_rows, elt);
   else
      return assigned_action(tree->right_child, data_x, num_rows, elt);
}

/** check whether a tree is 'perfect' for a set of units
 * if tree is not perfect a reason is printed to standard output
 * @return 1 if tree is perfect else 0
 */
int check_perfect(
   const NODE*           tree,               /**< allegedly perfect tree */
   int                   nunits,             /**< number of units for tree */
   const int*            units,              /**< units for tree */
   const double*         data_x,             /**< covariates, data_x+(j*num_rows) points to values for covariate j */
   const int*            best_actions,       /**< best_actions[i] is the best action for unit i */
   int                   num_rows            /**< number of units in entire dataset */
   )
{
   int i;

   for(i = 0; i < nunits; i++)
   {
      int elt = units[i];
      int assigned = assigned_action(tree, data_x, num_rows, elt);
      if( assigned != best_actions[elt] )
      {
         printf("Tree for %d elements is not, in fact, perfect!\n", nunits);
         printf("elt %d is assigned %d but best action is %d in following tree.\n", elt, assigned, best_actions[elt]);
         print_tree(tree, NULL);
         return 0;
      }
   }
   return 1;
}

/** delete a tree (free the memory it occupied) */
void tree_free(
   NODE*                 node                /**< root node of tree to free */
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




/** fix tree, remove 'dummy' nodes and then merge leaves with the same action
*/
void fix_tree(
   NODE*                 root                /**< root node */
   )
{
   prune_tree(root);
   merge_leaves(root);
}

/** return whether a node has both children */
int has_bothchildren(
   NODE*                 node                /**< node */
   )
{
   return node->left_child != NULL && node->right_child != NULL;
}


/** get the reward associated with a node */
double get_reward(
   NODE*                 node                /**< node */
   )
{
   return node->reward;
}

/** get the action associated with a node */
int get_action(
   NODE*                 node                /**< node */
   )
{
   return node->action_id;
}


/** set the reward associated with a node */
void set_reward(
   NODE*                 node,               /**< node */
   double                reward              /**< reward to associate with node */
   )
{
   node->reward = reward;
}


/** get the variable index associated with a node */
int get_index(
   NODE*                 node                /**< node */
   )
{
   return node->index;
}

/** get the split value associated with a node */
double get_value(
   NODE*                 node                /**< node */
   )
{
   return node->value;
}


/** get the children of a node */
void get_children(
   NODE*                 node,               /**< node */
   NODE**                left_child,         /**< pointer to left child */
   NODE**                right_child         /**< pointer to right child */
   )
{
   assert(node != NULL);
   
   *left_child = node->left_child;
   *right_child = node->right_child;
}
   
/** make a node a leaf and record reward and action for it */
void make_leaf(
   NODE*                 node,               /**< node */
   double                reward,             /**< total reward for units reaching this node */ 
   int                   action              /**< action for leaf */
   )
{
   assert(node != NULL); 
   
   node->index = -1;
   node->reward = reward;
   node->action_id = action;
}

/** record a variable split at a node */
void record_split(
   NODE*                 node,               /**< node */
   int                   split_var,          /**< (index of) var to split on */
   double                split_val,          /**< value of var to split on */
   double                reward              /**< total reward for units reaching this node */ 
   )
{
   assert(node != NULL); 
   assert(split_var > -1);

   node->index = split_var;
   node->value = split_val;
   node->reward = reward;
   node->action_id = -1;
}

/** record a variable split at a node and make its children leaves */
void record_level_one_split(
   NODE*                 node,               /**< node */
   int                   split_var,          /**< (index of) var to split on */
   double                split_val,          /**< value of var to split on */
   double                reward,             /**< total reward for units reaching this node */ 
   double                left_reward,        /**< reward for left leaf */
   int                   left_action,        /**< action for left leaf */
   double                right_reward,       /**< reward for right leaf */
   int                   right_action        /**< action for right leaf */
   )
{
   assert(node != NULL); 
   assert(node->left_child != NULL);
   assert(node->right_child != NULL);
   assert(split_var > -1);
   assert(left_action > -1);
   assert(right_action > -1);

   record_split(node, split_var, split_val, reward);

   node->left_child->index = -1;
   node->left_child->reward = left_reward;
   node->left_child->action_id = left_action;
   node->right_child->index = -1;
   node->right_child->reward = right_reward;
   node->right_child->action_id = right_action;
}

/** put a node and then all its descendants in the right place in an array of nodes */
static
void depth_first_nodes(
   NODE*                 node,               /**< current node */
   int                   depth,              /**< depth of current node */
   const int*            offset,             /**< offset for each depth */
   int*                  nd,                 /**< nd[d] is number of nodes at depth d already visited */ 
   NODE**                nodes               /**< array of nodes to populate */
   )
{
   /* put current node at correct position */
   nodes[offset[depth] + nd[depth]++] = node;

   if( node->left_child != NULL )
      depth_first_nodes(node->left_child, depth+1, offset, nd, nodes);

   if( node->right_child != NULL )
      depth_first_nodes(node->right_child, depth+1, offset, nd, nodes);
}

/** return nodes in a tree ordered breadth-first, left-right 
 * @return nodes in a tree ordered breadth-first, left-right 
 */
NODE** breadth_first_nodes(
   NODE*                 root,               /**< root node */
   int                   depth,              /**< depth of tree */
   int*                  num_nodes           /**< length of return array */
   )
{
   int i;
   int* offset = (int*) malloc((depth+1)*sizeof(int));
   int* nd;
   NODE** nodes;

   /* offset[i] is the number of nodes of depth less than i */
   *num_nodes = 1;
   for(i = 0; i <= depth; i++)
   {
      offset[i] = *num_nodes - 1;
      *num_nodes *= 2;
   }
   *num_nodes *= 2;
   (*num_nodes)--;

   nodes = (NODE**) calloc(*num_nodes,sizeof(NODE*));
   nd = (int*) calloc(depth+1,sizeof(int));

   depth_first_nodes(root, 0, (const int*) offset, nd, nodes);

   free(nd);
   free(offset);

   return nodes;
}

   
/** return the first (ie leftmost) node at a given depth, or NULL if there is none
 * @return the first (ie leftmost) node at a given depth, or NULL if there is none*/ 
NODE* first_at_depth(
   NODE*                 root,               /**< root node */
   int                   depth               /**< depth */
   )
{

   assert( root != NULL );
   assert( depth >= 0 );
   
   if( depth == 0 )
      return root;

   if( root->left_child == NULL )
      return NULL;

   return first_at_depth(root->left_child, depth-1);
}

/** return the next node at a given depth, or NULL if there is none
 * @return the next node at a given depth, or NULL if there is none*/ 
NODE* next_at_depth(
   NODE*                 root,               /**< root node */
   NODE*                 node,               /**< current node */
   int                   depth               /**< depth */
   )
{

   assert( root != NULL );
   assert( depth >= 0 );
   
   if( depth == 0 )
      return root;

   if( root->left_child == NULL )
      return NULL;

   return first_at_depth(root->left_child, depth-1);
}
