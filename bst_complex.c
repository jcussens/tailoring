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
   assert(node->flags >= 0);
   node->flags |= ELEMENT_PRESENT;
   assert(node->flags >= 0);
}

static
void mark_element_absent(
   BST* node
   )
{
   assert(node != NULL);
   assert(node->flags >= 0);
   node->flags &= ~ELEMENT_PRESENT;
   assert(node->flags >= 0);
}


static
void mark_left_branch_empty(
   BST* node
   )
{
   assert(node != NULL);
   assert(node->flags >= 0);
   node->flags &= ~LEFT_NONEMPTY;
   assert(node->flags >= 0);
}

static
void mark_right_branch_empty(
   BST* node
   )
{
   assert(node != NULL);
   assert(node->flags >= 0);
   node->flags &= ~RIGHT_NONEMPTY;
   assert(node->flags >= 0);
}

static
void mark_left_branch_nonempty(
   BST* node
   )
{
   assert(node != NULL);
   assert(node->flags >= 0);
   node->flags |= LEFT_NONEMPTY;
   assert(node->flags >= 0);
}

static
void mark_right_branch_nonempty(
   BST* node
   )
{
   assert(node != NULL);
   assert(node->flags >= 0);
   node->flags |= RIGHT_NONEMPTY;
   assert(node->flags >= 0);
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


static
int is_empty(
   BST* bst
   )
{
   assert(bst != NULL);
   
   if( element_present(bst) )
      return 0;

   if( bst->left != NULL && !is_empty(bst->left) )
      return 0;

   if( bst->right != NULL && !is_empty(bst->right) )
      return 0;

   return 1;
   
}

static
int tree_correct(
   BST* bst
   )
{
   assert(bst != NULL);
   
   if( bst->left == NULL )
   {
      if( !left_branch_marked_empty(bst) )
      {
         printf("NULL left branch of %p is not marked empty!", (void*)bst);
         return 0;
      }
   }
   else
   {
      if( left_branch_marked_empty(bst) )
      {
         if( !is_empty(bst->left) )
         {
            printf("left branch of %p is marked empty but it (%p) is not !", (void*)bst,(void*)bst->left);
            return 0;
         }
      }
      else
      {
         if( is_empty(bst->left) )
         {
            printf("left branch of %p is not marked empty but it (%p) is empty !", (void*)bst,(void*)bst->left);
            return 0;
         }
      }
      if( !tree_correct(bst->left) )
      {
         printf("left branch of %p (%p) is not correct somehow", (void*)bst,(void*)bst->left);
         return 0;
      }
   }

   if( bst->right == NULL )
   {
      if( !right_branch_marked_empty(bst) )
      {
         printf("NULL right branch of %p is not marked empty!", (void*)bst);
         return 0;
      }
   }
   else
   {
      if( right_branch_marked_empty(bst) )
      {
         if( !is_empty(bst->right) )
         {
            printf("right branch of %p is marked empty but it (%p) is not !", (void*)bst,(void*)bst->right);
            return 0;
         }
      }
      else
      {
         if( is_empty(bst->right) )
         {
            printf("right branch of %p is not marked empty but it (%p) is empty !", (void*)bst,(void*)bst->right);
            return 0;
         }
      }
      if( !tree_correct(bst->right) )
      {
         printf("right branch of %p (%p) is not correct somehow", (void*)bst,(void*)bst->right);
         return 0;
      }
   }
   return 1;
}

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

/** print out set elements in order without exploiting information about empty subtrees */
void in_order_tree_walk(
   BST* bst           /**< BST */
   )
{
   if( bst != NULL )
   {
      in_order_tree_walk(bst->left);
      if( element_present(bst) )
         printf("rank=%d,elt=%d ",bst->rank,bst->elt);
      in_order_tree_walk(bst->right);
   }
}


/** print out set elements in order using information about empty subtrees */
void in_order_tree_walk_fast(
   BST* bst           /**< BST */
   )
{
   if( bst != NULL )
   {
      if( left_branch_nonempty(bst) )
         in_order_tree_walk(bst->left);
      if( element_present(bst) )
         printf("rank=%d,elt=%d ",bst->rank,bst->elt);
      if( right_branch_nonempty(bst) )
         in_order_tree_walk(bst->right);
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


/** copy presence/absence data from source BST to target BST
 * (source and target must have same shape, element and rank info)
 */ 
void copy_data_bst(
   BST* source,   /**< source */
  BST* target          /**< target */
  )
{
  if( source != NULL )
  {
    assert(target != NULL);
    
    target->flags = source->flags;
    copy_data_bst(source->left,target->left);
    copy_data_bst(source->right,target->right);
    assert(tree_correct(target));
  }
}

/** make target BST represent empty set with same shape etc as source
 * (source and target must be the same shape)
 */ 
void empty_data_bst(
   BST* source,   /**< source */
   BST* target          /**< target */
  )
{
  if( source != NULL )
  {
    assert(target != NULL);
    
    target->flags = EMPTY_TREE;
    empty_data_bst(source->left,target->left);
    empty_data_bst(source->right,target->right);
    assert(tree_correct(target));
  }
}



/** return a copy of a binary search tree */
BST* copy_bst(
   BST* source,  /**< tree to be copied */
  BST* mother         /**< mother of copy */
  )
{

  BST* target = NULL;
  
  if( source != NULL )
  {
    target = (BST*) malloc(sizeof(BST));
    target->elt = source->elt;
    target->rank = source->rank;
    target->flags = source->flags;
    target->left = NULL;
    target->right = NULL;
    target->mother = mother;
    if( source->left != NULL )
      target->left = copy_bst(source->left,target);
    if( source->right != NULL )
      target->right = copy_bst(source->right,target);
  }
  return target; 
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
  bst->flags = 0;
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

  assert(tree_correct(bst));
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
   assert(node != NULL);
   assert(tree_correct(node));
   
  if( node->elt == elt )
  {
    mark_element_absent(node);
    empty_update_ancestors(node);
    if( !tree_correct(node) )
       print_tree(node);
    assert(tree_correct(node));
    return 1;
  }
  else if( rank < node->rank && node->left != NULL )
    return delete(node->left,elt,rank);
  else if(node->right != NULL )
    return delete(node->right,elt,rank);
  else
    return 0;
}


/** return root node from any other node */
BST* root_node(
   BST* node
   )
{
   assert(node != NULL);
   while(node->mother != NULL)
      node = node->mother;
   return node;
}

/** insert an element, returning 1 if successful else 0 */
int insert(
   BST* node, /**< tree containing element to insert */
   int elt,   /**< element to insert */
   int rank   /**< rank of element to insert */
   )
{
   assert(node != NULL);
   /* if( !tree_correct(node) ) */
   /*    print_tree(node); */
   assert(tree_correct(root_node(node)));
   
   if( node->elt == elt )
   {
      mark_element_present(node);
      nonempty_update_ancestors(node);
      assert(tree_correct(root_node(node)));
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
static
BST* find_next_alt(
  BST* x,  /**< node for current element */
  int* elt    /**< pointer to next element */
   )
{
   BST* y;
   
   assert(x != NULL);
   assert(!left_branch_nonempty(x) || x->left != NULL);
   assert(!right_branch_nonempty(x) || x->right != NULL);

   if( right_branch_nonempty(x) )
     return find_minimum(x->right,elt);

   y = x->mother;
   while( y != NULL && x == y->right )
   {
      x = y;
      y = y->mother;
   }

   if( y ==  NULL )
      return NULL;

   if( element_present(y) )
   {
      *elt = y->elt;
      return y;
   }
   else
      return find_next(y,elt);
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

   return find_next_alt(node,elt);
   
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

