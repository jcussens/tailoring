#ifdef __cplusplus
extern "C" {
#endif

struct bst;
typedef struct bst BST;

/** copy presence/absence data from source BST to target BST
 * (source and target must have same shape, element and rank info)
 */ 
void copy_data_bst(
   BST* source,   /**< source */
  BST* target          /**< target */
   );

/** make target BST represent empty set with same shape etc as source
 * (source and target must be the same shape)
 */ 
void empty_data_bst(
   BST* source,   /**< source */
  BST* target          /**< target */
   );



/** return a copy of a binary search tree */
BST* copy_bst(
    BST* source,  /**< tree to be copied */
   BST* mother         /**< mother of copy */
   );


/** make a binary search tree from an array of elements given in rank order */
BST* make_bst(
    int* elts,   /**< array of elements */
   int n,             /**< number of elements */
   int firstrank,     /**< rank of first element */
   BST* mother        /**< mother node (or NULL) */
   );

/** delete an element, returning 1 if successful else 0 */
int delete(
   BST* node,  /**< tree containing element to delete */
   int elt,    /**< element to delete */
   int rank    /**< rank of element to delete */
   );

/** insert an element, returning 1 if successful else 0 */
int insert(
   BST* node, /**< tree containing element to insert */
   int elt,   /**< element to insert */
   int rank   /**< rank of element to insert */
   );

/** find minimum element in a set, or NULL if set is empty */
BST* find_minimum(
   BST* node,  /**< node representing set */
   int* elt    /**< pointer to minimum element */
   );

/** find next element, or NULL if there is not one */
BST* find_next(
  BST* node,  /**< node for current element */
  int* elt    /**< pointer to next element */
   );

/** delete a binary search tree */
void delete_bst(
  BST* bst        /**< tree to be deleted */
   );

#ifdef __cplusplus
}
#endif
