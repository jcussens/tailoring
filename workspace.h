#ifndef __WORKSPACE_H__
#define __WORKSPACE_H__

#ifdef __cplusplus
extern "C" {
#endif

struct workspace;
typedef struct workspace WORKSPACE;  /**< Work space structure */

/** make workspace to provide pre-allocated space for various functions */
WORKSPACE* make_workspace(
   int                   depth,              /**< (maximum) depth of returned tree */
   SORTED_SET**          initial_sorted_sets, /**< initial sorted sets */
   int                   num_rows,           /**< number of units in full dataset */
   int                   num_cols_x,         /**< number of covariates */
   int                   num_cols_y          /**< number of rewards/actions */
   );

/** return array of doubles, one double for each reward */
double* get_double_array(
   WORKSPACE*            workspace           /**< workspace */
   );


#ifdef __cplusplus
}
#endif

#endif
