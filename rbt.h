#ifdef __cplusplus
extern "C" {
#endif

struct rbt;
typedef struct rbt RBT;

/** print the key-elements of a tree in key order */ 
void print_elts
(
   RBT* x,       /**< red-black tree */
   RBT* nil      /**< sentinel node */
   );


void print_rbt(
   RBT* x,       /**< red-black tree to print */
   RBT* nil      /**< sentinel node */
   );

/** make a red-black tree from an array of integers
 * each integer is stored at a node, its key is its index in the input array
 * @return the red-black tree
 */
RBT* make_rbt(
   int* elements,    /**< elements to store in tree, index in this array (+offset) is the key */
   int n,            /**< number of elements */
   RBT* nil          /**< sentinel node */
   );

/** Free memory used by a red-black tree
 */
void free_rbt(
   RBT* t,       /**< tree to be freed */
   RBT* nil      /**< sentinel node */
   );

/** Free memory used by sentinel node
 */
void free_sentinel(
   RBT* nil      /**< sentinel node */
   );


/** Return a sentinel node
 * @return a sentinel node
 */
RBT* sentinel(
   );

/** return the minimum element in a subtree rooted at a given node,
 * or nil if the tree is empty
 * @return the minimum element or nil if none
 */
RBT* tree_minimum(
   RBT* x,
   RBT* nil        /**< sentinel node */
   );


/** return the successor of a given node, or nil if there is none 
 * @return succesor node or nil if there is none
*/
RBT* tree_successor(
  RBT* x,         /**< node whose successor is sought */
  RBT* nil        /**< sentinel node */
   );

/** find the node with given key or nil if there isn't one 
 * @return node with given key, or nil
 */
RBT* iterative_tree_search(
  RBT* x,         /**< root node of tree */
  int key,        /**< key of sought node */
  RBT* nil        /**< sentinel node */
   );


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
   );


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
   );


/** get the element for a given node
 * @return the element for the given node
 */
int node_elt(
   RBT* node   /**< given node */
   );

#ifdef __cplusplus
}
#endif
