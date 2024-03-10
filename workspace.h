#ifndef __WORKSPACE_H__
#define __WORKSPACE_H__

#include "type_all.h"

#ifdef __cplusplus
extern "C" {
#endif

/** make workspace to provide pre-allocated space for various functions */
WORKSPACE* make_workspace(
   int                   depth,              /**< (maximum) depth of returned tree */
   SORTED_SET**          initial_sorted_sets, /**< initial sorted sets */
   int                   num_rows,           /**< number of units in full dataset */
   int                   num_cols_x,         /**< number of covariates */
   int                   num_cols_y          /**< number of rewards/actions */
   );

/** free the workspace */
void free_workspace(
   WORKSPACE*            workspace,          /**< workspace */
   int                   depth,              /**< (maximum) depth of returned tree */
   int                   num_cols_x          /**< number of covariates */
   );

/** return array of doubles, one double for each reward */
double* get_rewards(
   WORKSPACE*            workspace           /**< workspace */
   );

/** return array of zeroes, one zero for each reward */
double* get_rewards2(
   WORKSPACE*            workspace           /**< workspace */
   );


/** get left sorted sets associated with a given depth */
SORTED_SET** get_left_sorted_sets(
   WORKSPACE*            workspace,          /**< workspace */
   int                   depth               /**< depth */
   );

/** get right sorted sets associated with a given depth */
SORTED_SET** get_right_sorted_sets(
   WORKSPACE*            workspace,          /**< workspace */
   int                   depth               /**< depth */
   );


/** record a tree of given depth in workspace as the best tree */
void record_best_tree(
   WORKSPACE*            workspace,          /**< workspace */
   NODE*                 tree,               /**< tree */
   int                   depth               /**< depth of tree */
   );

/** retrieve the best tree of given depth from workspace */
void retrieve_best_tree(
   WORKSPACE*            workspace,          /**< workspace */
   NODE*                 tree,               /**< tree */
   int                   depth               /**< depth of tree */
   );


/** return an uninitialised sorted set */
SORTED_SET* get_sorted_set(
   WORKSPACE*            workspace           /**< workspace */
   );


#ifdef __cplusplus
}
#endif

#endif
