#ifdef __cplusplus
extern "C" {
#endif

struct node
{
  int index;                 /** splitting covariate */      
  double value;              /** split value */
  double reward;             /** reward (see above) */
  int action_id;             /** best action (only meaningful for leaves) */
  struct node* left_child;   /** left child (or NULL) */
  struct node* right_child;  /** right child (or NULL) */
};
typedef struct node NODE;

NODE* tree_search(
  const int depth,            /** (maximum) depth of returned tree */
  const int split_step,       /** consider splits every split_step'th possible split */
  const int min_node_size,    /** smallest terminal node size */
  const double* data_x,       /** covariates (column major) */
  const double* data_y,       /** gammas (column major) */
  const int num_rows,         /** number of units */
  const int num_cols_x,       /** number of covariates */
  const int num_cols_y        /** number of rewards */
   );

void tree_free(
   NODE*
   );

#ifdef __cplusplus
}
#endif
