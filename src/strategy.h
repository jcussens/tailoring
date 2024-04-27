/** @file strategy.h
 *  @brief Solving strategy
 *  @author James Cussens
 */

#ifndef __STRATEGY_H__
#define __STRATEGY_H__

#ifdef __cplusplus
extern "C" {
#endif

#include "type_all.h"

/** return an uninitialised strategy */
STRATEGY* get_unint_strategy(
   void
   );

/** are we using sorted sets for datasets? */
int using_sorted_sets(
   const STRATEGY*       strategy            /**< solving strategy */
   );

/** decide to use sorted sets */
void use_sorted_sets(
   STRATEGY*             strategy            /**< solving strategy */
   );

/** decide to use simple sets */
void use_simple_sets(
   STRATEGY*             strategy            /**< solving strategy */
   );

void decide_datatype(
   STRATEGY*             strategy,           /**< solving strategy */
   const double*         data_x,             /**< covariates, data_x+(j*num_rows) points to values for covariate j */
   int                   num_rows,           /**< number of units in full dataset */
   int                   num_cols_x          /**< number of covariates */
   );

#ifdef __cplusplus
}
#endif

#endif
