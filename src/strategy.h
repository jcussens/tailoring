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

/** return whether we wish to compute an upper bound on reward for each (sub-) dataset */
int find_reward_ub(
   const STRATEGY*       strategy            /**< solving strategy */
   );


/** set whether we wish to compute an upper bound on reward for each (sub-) dataset */
void set_find_reward_ub(
   STRATEGY*             strategy,           /**< solving strategy */
   int                   val                 /**< 0 for no, 1 for yes */
   );

/** return whether we wish to find dummy split rewards */
int find_dummy_split_reward(
   const STRATEGY*       strategy            /**< solving strategy */
   );

/** set whether we wish to find dummy split rewards */
void set_find_dummy_split_reward(
   STRATEGY*             strategy,           /**< solving strategy */
   int                   val                 /**< 0 for no, 1 for yes */
   );

/** return whether we wish to use last rewards rewards */
int use_last_rewards(
   const STRATEGY*       strategy            /**< solving strategy */
   );

/** set whether we wish to use last rewards */
void set_use_last_rewards(
   STRATEGY*             strategy,           /**< solving strategy */
   int                   val                 /**< 0 for no, 1 for yes */
   );

/** return whether we wish to use cutoffs */
int use_cutoffs(
   const STRATEGY*       strategy            /**< solving strategy */
   );

/** set whether we wish to use cutoffs */
void set_use_cutoffs(
   STRATEGY*             strategy,           /**< solving strategy */
   int                   val                 /**< 0 for no, 1 for yes */
   );


#ifdef __cplusplus
}
#endif

#endif
