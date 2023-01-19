/* red-black trees all code taken directly from Cormen et al (elt is the single additional field)
 */

#include <assert.h>
#include <stdlib.h>
#include "rbt.h"

#define RED '\000'
#define BLACK '\001'

struct rbt
{
  int key;              /**< node key */
  char colour;          /**< colour (red or black ) */
  struct rbt* left;     /**< left subtree (nil if absent) */
  struct rbt* right;    /**< right subtree (nil if absent) */
  struct rbt* p;        /**< parent node (or nil if absent) */
  int elt;              /**< element stored at this node */
};


/** Free memory used by a red-black tree
 */
void free_rbt(
   RBT* t,       /**< tree to be freed */
   RBT* nil      /**< sentinel node */
   )
{
   if( t != nil )
   {
      free_rbt(t->left,nil);
      free_rbt(t->right,nil);
      free(t);
   }
}

/** Free memory used by sentinel node
 */
void free_sentinel(
   RBT* nil      /**< sentinel node */
   )
{
   free(nil);
}


/** Return a sentinel node
 * @return a sentinel node
 */
RBT* sentinel(
  )
{
  RBT* nil = (RBT*) malloc(sizeof(RBT));
  nil->colour = BLACK;
  return nil;
}

/** Rotate a node left, assuming that the node has a right child
 * @return the root of the tree after rotation
 */
static
RBT* left_rotate(
  RBT* t,         /**< root of the tree */
  RBT* x,         /**< node to rotate left */
  RBT* nil        /**< sentinel node */
  )
{
  RBT* y;

  assert(t != NULL);
  assert(x != NULL);
  assert(nil != NULL);

  assert(t->p == nil);
  assert(x->right != nil);

  y = x->right;
  x->right = y->left;
  if( y->left != nil )
    y->left->p = x;
  y->p = x->p;
  
  if( x->p == nil )
    t = y;
  else if( x == x->p->left )
    x->p->left = y;
  else
    x->p->right = y;
  
  y->left = x;
  x->p = y;

  return t;
}

/** Rotate a node right, assuming that the node has a left child
 * @return the root of the tree after rotation
 */
static
RBT* right_rotate(
  RBT* t,         /**< root of the tree */
  RBT* x,         /**< node to rotate right */
  RBT* nil        /**< sentinel node */
  )
{
  RBT* y;

  assert(t != NULL);
  assert(x != NULL);
  assert(nil != NULL);

  assert(t->p == nil);
  assert(x->left != nil);

  y = x->left;
  x->left = y->right;
  if( y->right != nil )
    y->right->p = x;
  y->p = x->p;
  
  if( x->p == nil )
    t = y;
  else if( x == x->p->right )
    x->p->right = y;
  else
    x->p->left = y;
  
  y->right = x;
  x->p = y;
  
  return t;
}

/* /\** Get the key and element stored at a node  */
/*  *\/ */
/* void tree_minimum( */
/*    RBT* node,   /\**< node *\/ */
/*    int* key,    /\**< on return *key will be key for node *\/ */
/*    int* elt     /\**< on returh *elt will be elt for key for node *\/ */
/*   ) */
/* { */
/*    assert(node != NULL); */
/*    assert(key != NULL); */
/*    assert(elt != NULL); */
   
/*    *key = node->key; */
/*    *elt = node->elt; */
/* } */

/** return the minimum element in a subtree rooted at a given node,
 * or nil if the tree is empty
 * @return the minimum element or nil if none
 */
RBT* tree_minimum(
  RBT* x,
  RBT* nil        /**< sentinel node */
  )
{
  assert(x != NULL);
  assert(nil != NULL);
  
  while( x->left != nil)
    x = x->left;

  return x;
}

/** return the successor of a given node, or nil if there is none 
 * @return succesor node or nil if there is none
*/
RBT* tree_successor(
  RBT* x,         /**< node whose successor is sought */
  RBT* nil        /**< sentinel node */
  )
{
  RBT* y;
  
  assert(x != NULL);
  assert(nil != NULL);

  assert(x != nil);
  
  if( x->right != nil )
     return tree_minimum(x->right,nil);

  y = x->p;
  while( y != nil && x == y->right )
  {
    x = y;
    y = y->p;
  }
  return y;
}

/** insert a node into a possibly empty tree and return root node
 * t=nil represent an empty tree 
 * Assumes node to be inserted has nil left and right children
 * @return the root of the tree after insertion
 */
static
RBT* tree_insert(
  RBT* t,         /**< root node of tree to insert into */
  RBT* z,         /**< node to insert */
  RBT* nil        /**< sentinel node */
  )
{
  RBT* y;
  RBT* x;

  assert(t != NULL);
  assert(z != NULL);
  assert(nil != NULL);

  assert(t->p == nil);
  assert(z != nil);
  assert(z->left == nil);
  assert(z->right == nil);

  y = nil;
  x = t;
  
  while( x != nil )
  {
    y = x;
    if( z->key < x->key )
      x = x->left;
    else
      x = x->right;
  }

  z->p = y;

  if( y == nil)
    return z;
  else
  {
    if( z->key < y->key )
      y->left = z;
    else
      y->right = z;
    return t;
  }
}

/** find the node with given key or nil if there isn't one 
 * @return node with given key, or nil
 */
RBT* iterative_tree_search(
  RBT* x,         /**< root node of tree */
  int key,        /**< key of sought node */
  RBT* nil        /**< sentinel node */
  )
{
  while( x != nil && x->key != key )
    if( key < x->key )
      x = x->left;
    else
      x = x->right;

  return x;
}

/** insert a node into a possibly empty tree and return root node,
 * and preserve red-black properties.
 * t=nil represent an empty tree
 * only the field values of the inserted node matter
 * @return the root of the tree after insertion
 */
RBT* rb_insert(
  RBT* t,         /**< root node of tree to insert into */
  RBT* x,         /**< node to insert */
  RBT* nil        /**< sentinel node */
  )
{

   RBT* y;

  assert(t != NULL);
  assert(x != NULL);
  assert(nil != NULL);

  assert(t->p == nil);
  assert(x != nil);

  x->left = nil;
  x->right = nil;

  t = tree_insert(t, x, nil);
  
  x->colour = RED;
  while( x->p != nil && x->p->colour == RED )
  {
    if( x->p == x->p->p->left )
    {
      y = x->p->p->right;
      if( y->colour == RED )
      {
        x->p->colour = BLACK;
        y->colour = BLACK;
        x->p->p->colour = RED;
        x = x->p->p;
      }
      else
      {
        if( x == x->p->right )
        {
          x = x->p;
          t = left_rotate(t, x, nil);
        }
        x->p->colour = BLACK;
        x->p->p->colour = RED;
        t = right_rotate(t, x->p->p, nil);
      }
    }
    else
    {
      y = x->p->p->left;
      if( y->colour == RED )
      {
        x->p->colour = BLACK;
        y->colour = BLACK;
        x->p->p->colour = RED;
        x = x->p->p;
      }
      else
      {
        if( x == x->p->left )
        {
          x = x->p;
          t = right_rotate(t, x, nil);
        }
        x->p->colour = BLACK;
        x->p->p->colour = RED;
        t = left_rotate(t, x->p->p, nil);
      }
    }
  }
  t->colour = BLACK;
  return t;
}

/** Fix a tree after deletion to restore red-black properties
 * @return the root of the tree after rotation
 */
static
RBT* rb_delete_fixup(
   RBT* t,         /**< root node of tree to correct */
   RBT* x,         /**< node to correct (could be the sentinel) */
   RBT* nil        /**< sentinel node */
   )
{

   RBT* w;
   
   assert(t != NULL);
   assert(x != NULL);
   assert(nil != NULL);
   
   assert(t->p == nil);

   while( x != t && x->colour == BLACK )
   {
      if( x->p->left == x )
      {
         w = x->p->right;

         if( w->colour == RED )
         {
            w->colour = BLACK;
            x->p->colour = RED;
            t = left_rotate(t, x->p, nil);
            w = x->p->right;
         }

         if( w->left->colour == BLACK && w->right->colour == BLACK )
         {
            w->colour = RED;
            x = x->p;
         }
         else
         {
            if( w->right->colour == BLACK )
            {
               w->left->colour = BLACK;
               w->colour = RED;
               t = right_rotate(t, w, nil);
               w = x->p->right;
            }

            w->colour = x->p->colour;
            x->p->colour = BLACK;
            w->right->colour = BLACK;
            t = left_rotate(t, x->p, nil);
            x = t;
         }
      }
      else
      {
         w = x->p->left;

         if( w->colour == RED )
         {
            w->colour = BLACK;
            x->p->colour = RED;
            t = right_rotate(t, x->p, nil);
            w = x->p->left;
         }

         if( w->right->colour == BLACK && w->left->colour == BLACK )
         {
            w->colour = RED;
            x = x->p;
         }
         else
         {
            if( w->left->colour == BLACK )
            {
               w->right->colour = BLACK;
               w->colour = RED;
               t = left_rotate(t, w, nil);
               w = x->p->left;
            }

            w->colour = x->p->colour;
            x->p->colour = BLACK;
            w->left->colour = BLACK;
            t = right_rotate(t, x->p, nil);
            x = t;
         }
      }
   }
   x->colour = BLACK;
   return t;
}

/* delete a node from a tree and return root node
 * the 'spliced out' node is also 'returned' via a pointer
 * If the node to delete z has two children then its successor y is spliced out
 * and the contents of y then overwrite the contents of z, otherwise z itself
 * is spliced out. The 'spliced out' node has same fields as z in any event.
 * @return the root of the tree after deletion
 */
RBT* rb_delete(
   RBT* t,         /**< root node of tree to insert into */
   RBT* z,         /**< node to delete */
   RBT** r,        /**< (pointer to) 'spliced out' node */
   RBT* nil        /**< sentinel node */
  )
{

   RBT* x;
   RBT* y;

   assert(t != NULL);
   assert(z != NULL);
   assert(r != NULL);
   assert(*r != NULL);
   assert(nil != NULL);
   
  assert(t->p == nil);
  assert(z != nil);
   
   if( z->left == nil || z->right == nil )
      y = z;
  else
     y = tree_successor(z, nil);
   
   if( y->left != nil )
      x = y->left;
   else
      x = y->right;
   
   x->p = y->p;
   
   if( y->p == nil )
      t = x;
   else if( y == y->p->left )
      y->p->left = x;
   else
      y->p->right = x;
   
   if( y != z )
   {
     /* if y != z then z is 'deleted' by having its
        fields overwritten by those of y
        but we also want y (the actually deleted, returned node)
        to have the fields that z has, so we swap fields
     */
     int tmp;

     tmp = z->key;
     z->key = y->key;
     y->key = tmp;

     tmp = z->elt;
     z->elt = y->elt;
     y->elt = tmp;
   }

   if( y->colour == BLACK )
      t = rb_delete_fixup(t, x, nil );

   *r = y;
   return t;
}


/** make a red-black tree from an array of integers
 * each integer is stored at a node, its key is its index in the input array
 * @return the red-black tree
 */
static
RBT* make_rbt0(
   int* elements,    /**< elements to store in tree, index in this array (+offset) is the key */
   int n,            /**< number of elements */
   RBT* p,           /**< parent of root node */
   int offset,       /**< offset to get key values correct */
   int depth,        /**< depth of tree within global tree */
   int maxdepth,     /**< depth of deepest non-nil nodes */
   RBT* nil          /**< sentinel node */
   )
{
   RBT* node;
   int mididx;

   if( n == 0)
      return nil;
   
   /* if n is odd,  mididx indexes the middle element, split is (n-1)/2,1,(n-1)/2 
    * if n is even, mididx indexes the greater of the two middle elements, split is n/2,1,n/2 - 1 
   */
   mididx = n/2;

   /* make node node */
   node = (RBT*) malloc(sizeof(RBT));
   node->key = mididx+offset;
   node->colour = ( depth == maxdepth ) ? RED : BLACK;
   node->elt = elements[mididx];
   node->p = p;
   node->left = make_rbt0( elements,          mididx,     node, offset,          depth+1, maxdepth, nil);
   node->right = make_rbt0(elements+mididx+1, n-mididx-1, node, offset+mididx+1, depth+1, maxdepth, nil);

   return node;
}

/** make a red-black tree from an array of integers
 * each integer is stored at a node, its key is its index in the input array
 * @return the red-black tree
 */
RBT* make_rbt(
   int* elements,    /**< elements to store in tree, index in this array (+offset) is the key */
   int n,            /**< number of elements */
   RBT* nil          /**< sentinel node */
   )
{
   int maxdepth = 0;
   int nleft = n;

   while( nleft > 1 )
   {
      nleft /= 2;
      maxdepth++;
   }

   return make_rbt0(elements, n, nil, 0, 0, maxdepth, nil);
}

/** get the element for a given node
 * @return the element for the given node
 */
int node_elt(
   RBT* node   /**< given node */
   )
{
   return node->elt;
}
