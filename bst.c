/*
 * 
 *
 */

#include "bst.h"
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>

#define EMPTY_TREE '\000'
#define LEFT_NONEMPTY '\004'
#define RIGHT_NONEMPTY '\002'
#define ELEMENT_PRESENT '\001'

struct bst
{
   int elt;              /**< element stored at this node */
   int rank;             /**< rank of the element */
   char flags;           /**< bitset indicating absence/presence */
   struct bst* left;     /**< left subtree (or NULL if absent) */
   struct bst* right;    /**< right subtree (or NULL if absent) */
   struct bst* mother;   /**< mother node (or NULL if absent) */
};
typedef struct bst BST;  

void print_node(
  BST* bst
  )
{
  if( bst != NULL )
  {
    printf("node=%p,elt=%d,present=%d,rank=%d,flags=%d\nleft_nonempty=%d,right_nonempty=%d,left_child=%p,right_child=%p,mother=%p\n\n",
      (void*)bst,bst->elt,bst->flags & ELEMENT_PRESENT,bst->rank,bst->flags,
      bst->flags & LEFT_NONEMPTY,bst->flags & RIGHT_NONEMPTY,bst->left,bst->right,bst->mother);
  }
}

void print_tree(
  BST* bst
  )
{

  print_node(bst);
  
  if( bst != NULL )
  {
    print_tree(bst->left);
    print_tree(bst->right);
  }
}


/* return 1 if node represents an empty set, else 0 */
static
int empty_set(
   BST* node
   )
{
   assert(node != NULL);
   return node->flags == EMPTY_TREE;
}

/* return 4 if left branch of a non-null node is non-empty, else 0 */
static
int left_branch_nonempty(
   BST* node
   )
{
   assert(node != NULL);
   return node->flags & LEFT_NONEMPTY;
}

/* return 2 if right branch of a non-null node is non-empty, else 0 */
static
int right_branch_nonempty(
   BST* node
   )
{
   assert(node != NULL);
   return node->flags & RIGHT_NONEMPTY;
}


/* return 1 if element for non-null node is present in the set, else 0 */
static
int element_present(
   BST* node
   )
{
   assert(node != NULL);
   return node->flags & ELEMENT_PRESENT;
}

/* return 1 if node is the root node else 0 */
static
int is_root(
   BST* node
   )
{
   assert(node != NULL);
   return node->mother == NULL;
}

/* return 1 if node has a mother and is the left child of its mother else 0 */
static
int is_left_child(
   BST* node
   )
{
   assert(node != NULL);
   return node->mother == NULL || node->mother->left == node;
}

static
void mark_element_present(
   BST* node
   )
{
   assert(node != NULL);
   node->flags |= ELEMENT_PRESENT;
}

static
void mark_element_absent(
   BST* node
   )
{
   assert(node != NULL);
   node->flags &= ~ELEMENT_PRESENT;
}


static
void mark_left_branch_empty(
   BST* node
   )
{
   assert(node != NULL);
   node->flags &= ~LEFT_NONEMPTY;   
}

static
void mark_right_branch_empty(
   BST* node
   )
{
   assert(node != NULL);
   node->flags &= ~RIGHT_NONEMPTY;   
}

static
void mark_left_branch_nonempty(
   BST* node
   )
{
   assert(node != NULL);
   node->flags |= LEFT_NONEMPTY;   
}

static
void mark_right_branch_nonempty(
   BST* node
   )
{
   assert(node != NULL);
   node->flags |= RIGHT_NONEMPTY;   
}

static
int left_branch_marked_empty(
   BST* node
   )
{
   assert(node != NULL);
   return !(node->flags & LEFT_NONEMPTY);
}

static
int right_branch_marked_empty(
   BST* node
   )
{
   assert(node != NULL);
   return !(node->flags & RIGHT_NONEMPTY);
}

/** delete a binary search tree */
void delete_bst(
  BST* bst        /**< tree to be deleted */
  )
{
  if( bst->left != NULL)
    delete_bst(bst->left);
  if( bst->right != NULL)
    delete_bst(bst->right);
  free(bst);
}

/** make a binary search tree from an array of elements given in rank order */
BST* make_bst(
  int* elts,     /**< array of elements */
  int n,         /**< number of elements */
  int firstrank, /**< rank of first element */
  BST* mother    /**< mother node (or NULL) */
  )
{

  BST* bst;
  int mididx;
  
  if( n == 0)
    return NULL;

  /* if e.g. elts = {4,5,7,3,6}, so n = 5, we want 7 at index 5/2 = 2 to be
     the node element*/

  bst = (BST*) malloc(sizeof(BST));

  mididx = n/2;
  bst->elt = elts[mididx];
  bst->rank = firstrank+mididx;
  mark_element_present(bst);
  bst->left = NULL;
  bst->right = NULL;
  bst->mother = mother;
  if( mididx > 0 )
  {
    /* elts[0,...,mididx-1] goes left, mididx-1 elements */
    
    mark_left_branch_nonempty(bst);
    bst->left = make_bst(elts,mididx,firstrank,bst);
  }
  if( n-mididx-1 > 0 )
  {
    /* elts[mididx+1,...,n-1] goes right, n-(mididx+1) elements  */
    
    mark_right_branch_nonempty(bst);
    bst->right = make_bst(elts+mididx+1,n-mididx-1,firstrank+mididx+1,bst);
  }
  return bst;
}



/** if node is an empty tree, update ancestors, if any */
static
void empty_update_ancestors(
   BST* node
   )
{
   
   /* if tree is empty and not the root node */
   if( empty_set(node) && !is_root(node) )
   {
      /* update mother flags */ 
   
      if( is_left_child(node) )
         mark_left_branch_empty(node->mother);
      else
         mark_right_branch_empty(node->mother);

      empty_update_ancestors(node->mother);
   }
}

/** if node is an nonempty tree, update ancestors, if any */
static
void nonempty_update_ancestors(
   BST* node
   )
{
   
   /* if tree is non-empty and not the root node */
   if( !empty_set(node) && !is_root(node))
   {
      /* update mother flags, if necessary */
      /* only recurse if updating was done */
   
      if( is_left_child(node) && left_branch_marked_empty(node->mother))
         mark_left_branch_nonempty(node->mother);
      else if (right_branch_marked_empty(node->mother))
         mark_right_branch_nonempty(node->mother);
      else
         return;

      nonempty_update_ancestors(node->mother);
   }
}



/** delete an element, returning 1 if successful else 0 */
int delete(
  BST* node,  /**< tree containing element to delete */
  int elt,    /**< element to delete */
  int rank    /**< rank of element to delete */
  )
{
  if( node->elt == elt )
  {
    mark_element_absent(node);
    empty_update_ancestors(node);    
    return 1;
  }
  else if( rank < node->rank && node->left != NULL )
    return delete(node->left,elt,rank);
  else if(node->right != NULL )
    return delete(node->right,elt,rank);
  else
    return 0;
}

/** insert an element, returning 1 if successful else 0 */
int insert(
   BST* node, /**< tree containing element to insert */
   int elt,   /**< element to insert */
   int rank   /**< rank of element to insert */
   )
{
   if( node->elt == elt )
   {
      mark_element_present(node);
      nonempty_update_ancestors(node);     
      return 1;
   }
   else if( rank < node->rank && node->left != NULL )
      return insert(node->left,elt,rank);
   else if(node->right != NULL )
      return insert(node->right,elt,rank);
   else
      return 0;
}

/** find minimum element in a set, or NULL if set is empty */
BST* find_minimum(
  BST* node,  /**< node representing set */
  int* elt    /**< pointer to minimum element */
  )
{
  assert(node != NULL);
  assert(!left_branch_nonempty(node) || node->left != NULL);
  assert(!right_branch_nonempty(node) || node->right != NULL);
  
  if( left_branch_nonempty(node) )
    return find_minimum(node->left,elt);
  else if ( element_present(node) )
  {
    *elt = node->elt;
    return node;
  }
  else if( right_branch_nonempty(node) )
    return find_minimum(node->right,elt);
  else
    return NULL;
}


/** find next element, or NULL if there is not one */
BST* find_next(
  BST* node,  /**< node for current element */
  int* elt    /**< pointer to next element */
   )
{

   assert(node != NULL);
   assert(!left_branch_nonempty(node) || node->left != NULL);
   assert(!right_branch_nonempty(node) || node->right != NULL);

   if( right_branch_nonempty(node) )
     return find_minimum(node->right,elt);
   else
   {
     while( node != NULL && (is_root(node) || !is_left_child(node)) )
       node = node->mother;

     if( node ==  NULL)
       return NULL;

     node = node->mother;
     
     if( element_present(node) )
     {
       *elt = node->elt;
       return node;
     }
     else
       return find_next(node,elt);
   }
}

