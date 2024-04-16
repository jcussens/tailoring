/** @file simple_opttree.h
 *  @brief Functions for building a tree
 *  @author James Cussens
 */

#ifndef __SIMPLE_OPTTREE_H__
#define __SIMPLE_OPTTREE_H__

#include "tree.h"

#ifdef __cplusplus
extern "C" {
#endif

/** Find an optimal policy tree of given maximal depth 
 * @return An optimal policy tree
 */
NODE* tree_search_simple(
  int                    depth,              /**< (maximum) depth of returned tree */
  int                    min_node_size,      /**< smallest terminal node size */
  const double*          data_x,             /**< covariates, data_x+(j*num_rows) points to values for covariate j */
  const double*          data_y,             /**< 'gammas', data_y+(d*num_rows) points to values for reward d */
  int                    num_rows,           /**< number of units in full dataset */
  int                    num_cols_x,         /**< number of covariates */
  int                    num_cols_y,         /**< number of actions */
  double*                reward              /**< reward for optimal tree */
   );

#ifdef __cplusplus
}
#endif

#endif
