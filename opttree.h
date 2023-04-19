#ifdef __cplusplus
extern "C" {
#endif

/** node of a policytree (for a particular dataset of units)
 *  the associated dataset of units is not represented, only the 
 *  reward for the associated dataset for the tree with the node as root
 */
struct node
{
  int                    index;              /** splitting covariate */      
  double                 value;              /** split value */
  double                 reward;             /** reward for tree with this node on an unrepresented associated dataset of units */
  int                    action_id;          /** best action (only meaningful for leaves) */
  struct node*           left_child;         /** left child (or NULL) */
  struct node*           right_child;        /** right child (or NULL) */
};
typedef struct node NODE;


/** Find an optimal policy tree of given maximal depth 
 * @return An optimal policy tree
 */
NODE* tree_search_jc_policytree(
  int                    depth,              /**< (maximum) depth of returned tree */
  int                    split_step,         /**< consider splits every split_step'th possible split */
  int                    min_node_size,      /**< smallest terminal node size */
  const double*          data_x,             /**< covariates, data_x+(j*num_rows) points to values for covariate j */
  const double*          data_y,             /**< gammas, data_y+(d*num_rows) points to values for reward d */
  int                    num_rows,           /**< number of units in full dataset */
  int                    num_cols_x,         /**< number of covariates */
  int                    num_cols_y          /**< number of rewards/actions */
   );

/** Find an optimal policy tree of given maximal depth 
 * @return An optimal policy tree
 */
NODE* tree_search_jc_discretedata(
  int                    depth,              /**< (maximum) depth of returned tree */
  int                    split_step,         /**< consider splits every split_step'th possible split */
  int                    min_node_size,      /**< smallest terminal node size */
  const double*          data_x,             /**< covariates, data_x+(j*num_rows) points to values for covariate j */
  const double*          data_y,             /**< gammas, data_y+(d*num_rows) points to values for reward d */
  int                    num_rows,           /**< number of units in full dataset */
  int                    num_cols_x,         /**< number of covariates */
  int                    num_cols_y          /**< number of rewards/actions */
   );

void tree_free(
   NODE*
   );

#ifdef __cplusplus
}
#endif
