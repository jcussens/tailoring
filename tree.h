#ifndef __TREE_H__
#define __TREE_H__

#include "type_all.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Return a 'skeleton' policy tree of required depth
 * or NULL if there is insufficient memory
 */
NODE* make_tree(
   int                    depth               /**< depth of required tree */
   );

/** delete a tree (free the memory it occupied) */
void tree_free(
   NODE*                 node                /**< root node of tree to free */
   );

/**
 * copy data from source tree to existing target tree
 * assumes both trees have same shape.
 * NB. No new memory is allocated
 */
void tree_copy(
   const NODE*           source,             /**< source tree */
   NODE*                 target              /**< target tree */
   );

/** add a given reward to those nodes in a tree which a particular element visits */
void update_rewards(
   NODE*                 tree,               /**< tree */
   const double*         data_x,             /**< covariates, data_x+(j*num_rows) points to values for covariate j */
   int                   num_rows,           /**< number of units */
   int                   elt,                /**< element */
   double                reward              /**< reward */
   );


/** find the action assigned to a unit by a tree 
 * @return the assigned action
 */
int assigned_action(
   const NODE*           tree,               /**< tree */
   const double*         data_x,             /**< covariates, data_x+(j*num_rows) points to values for covariate j */
   int                   num_rows,           /**< number of units */
   int                   elt                 /**< element to find assigned action for */
   );

int check_perfect(
   const NODE*           tree,               /**< allegedly perfect tree */
   int                   nunits,             /**< number of units for tree */
   const int*            units,              /**< units for tree */
   const double*         data_x,             /**< covariates, data_x+(j*num_rows) points to values for covariate j */
   const int*            best_actions,       /**< best_actions[i] is the best action for unit i */
   int                   num_rows            /**< number of units in entire dataset */
   );


/** return whether a node has both children */
int has_bothchildren(
   NODE*                 node                /**< node */
   );

/* get the reward associated with a node */
double get_reward(
   NODE*                 node                /**< node */
   );

/* set the reward associated with a node */
void set_reward(
   NODE*                 node,               /**< node */
   double                reward              /**< reward to associate with node */
   );


/* get the variable index associated with a node */
int get_index(
   NODE*                 node                /**< node */
   );

/* get the split value associated with a node */
double get_value(
   NODE*                 node                /**< node */
   );

   /* get the action associated with a node */
int get_action(
   NODE*                 node                /**< node */
   );

/* get the children of a node */
void get_children(
   NODE*                 node,               /**< node */
   NODE**                left_child,         /**< pointer to left child */
   NODE**                right_child         /**< pointer to right child */
   );

/** make a node a leaf and record reward and action for it */
void make_leaf(
   NODE*                 node,               /**< node */
   double                reward,             /**< total reward for units reaching this node */ 
   int                   action              /**< action for leaf */
   );

/** record a variable split at a node */
void record_split(
   NODE*                 node,               /**< node */
   int                   split_var,          /**< (index of) var to split on */
   double                split_val,          /**< value of var to split on */
   double                reward              /**< total reward for units reaching this node */ 
   );

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
   );

/* is a node a leaf node? */
int is_leaf(
   NODE*                 node                /**< node */
   );

/** prints a policy tree in policytree style to standard output */
void print_tree_policytree(
   NODE*                 tree,               /**< policy tree to print */
   char**                covnames,           /**< covnames[j] is the name of covariate j */
   int                   depth,              /**< depth of the tree */
   int                   nactions,           /**< number of actions */
   char**                actionnames         /**< covnames[j] is the name of action j */
   );

/** print a policy tree (debugging only) 
 * if covariate names are not supplied then the indices for covariates are used
 */
void print_tree(
   const NODE*           tree,               /**< root of tree to print */ 
   const char**          covnames            /**< if covnames != NULL, then covnames[i] is the name of covariate i */
   );

/** prints a policy tree to standard output */
void print_tree_raw(
   NODE*                 tree               /**< policy tree to print */
   );

/** fix tree, remove 'dummy' nodes and then merge leaves with the same action
*/
void fix_tree(
   NODE*                 root                /**< root node */
   );


#ifdef __cplusplus
}
#endif

#endif
