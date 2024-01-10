struct workspace
{
   double*               double_array;       /**< size is number of actions */
   SORTED_SET***         lefts;              /**< a sorted set for each depth and each covariate, each with space = number of units */
   SORTED_SET***         rights;             /**< a sorted set for each depth and each covariate, each with space = number of units */
};

/** make workspace to provide pre-allocated space for various functions */
WORKSPACE* make_workspace(
   int                   depth,              /**< (maximum) depth of returned tree */
   SORTED_SET**          initial_sorted_sets, /**< initial sorted sets */
   int                   num_rows,           /**< number of units in full dataset */
   int                   num_cols_x,         /**< number of covariates */
   int                   num_cols_y          /**< number of rewards/actions */
   )
{
   WORKSPACE* workspace;

   workspace->double_array = (double*) malloc(num_cols_y*sizeof(double));

   return workspace;
}

/** return array of doubles, one double for each action */
double* get_double_array(
   WORKSPACE*            workspace           /**< workspace */
   )
{
   return workspace->double_array:
}

SORTED_SET** get_left_sorted_sets(
   WORKSPACE*            workspace,          /**< workspace */
   int                   depth
   )
{
   return workspace->lefts[depth];
}

SORTED_SET** get_right_sorted_sets(
   WORKSPACE*            workspace,          /**< workspace */
   int                   depth
   )
{
   return workspace->rights[depth];
}
