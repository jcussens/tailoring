void tree_search(
  int depth,            /** (maximum) depth of returned tree */
  int split_step,       /** consider splits every split_step'th possible split */
  int min_node_size,    /** smallest terminal node size */
  double* data_x,       /** covariates (column major) */
  double* data_y,       /** gammas (column major) */
  int num_rows,         /** number of units */
  int num_cols_x,       /** number of covariates */
  int num_cols_y        /** number of rewards */
   );
