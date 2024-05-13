/** @file units.h
 *  @brief Dispatches to appropriate data type for units
 *  @author James Cussens
 */
#ifndef __UNITS_H__
#define __UNITS_H__

#ifdef __cplusplus
extern "C" {
#endif

#include "type_all.h"

/** make unit set from covariate data */
UNITS make_units(
   const STRATEGY*       strategy,           /**< tree-building strategy */
   const double*         data_x,             /**< covariates, data_x+(j*num_rows) points to values for covariate j */
   int                   num_rows,           /**< number of units in full dataset */
   int                   num_cols_x          /**< number of covariates */
   );

/** free (space for) units */
void free_units(
   const STRATEGY*       strategy,           /**< tree-building strategy */
   UNITS                 units,              /**< units */
   int                   num_cols_x          /**< number of covariates */
   );

/** make a 'shallow' copy of units, these copies should be freed by calling shallow_free_units */
UNITS shallow_copy_units(
   const STRATEGY*       strategy,           /**< tree-building strategy */
   CONST_UNITS           units,              /**< source units */
   int                   num_cols_x          /**< number of covariates */
   );

/** free (space for) units previously created as shallow copies */
void shallow_free_units(
   const STRATEGY*       strategy,           /**< tree-building strategy */
   UNITS                 units,              /**< units */
   int                   num_cols_x          /**< number of covariates */
   );

/** for debugging only: check that a non-empty collection of sorted sets represent the same underlying set and that each is appropriately sorted */
int units_ok(
   const STRATEGY*       strategy,           /**< tree-building strategy */
   CONST_UNITS           units,              /**< units */
   int                   p,                  /**< if !=-1 then units must be ready for splitting on covariate p */
   const double*         data_x,             /**< covariates, data_x+(j*num_rows) points to values for covariate j */
   int                   num_rows,           /**< number of units in full dataset */
   int                   num_cols_x          /**< number of covariates */
   );

/** get the reward for a set if all units in the set were assigned their best action *
 * @return the reward for a set if all units in the set were assigned their best action *
 */
double get_reward_ub(
   const STRATEGY*       strategy,           /**< tree-building strategy */
   CONST_UNITS           units,              /**< set */
   const double*         data_y,             /**< gammas, data_y+(d*num_rows) points to values for reward d */
   int                   num_rows,           /**< number of units in full dataset */
   const int*            best_actions        /**< best_actions[i] is the best action for unit i */
   );


/** Determine whether a set is 'pure'.
 * A pure set is one where each unit has the same best action
 * @return 1 if the set is pure, else 0
 */
int is_pure(
   const STRATEGY*       strategy,           /**< tree-building strategy */
   CONST_UNITS           units,              /**< units */
   const int*            best_actions        /**< best_actions[i] is the best action for unit i */
   );

/** for each action find (total) reward if that action applied to each unit in a set of units */
void find_nosplit_rewards(
   const STRATEGY*       strategy,           /**< tree-building strategy */
   CONST_UNITS           units,              /**< units */
   int                   num_cols_y,         /**< number of actions */
   const double*         data_y,             /**< gammas, data_y+(d*num_rows) points to values for reward d */
   int                   num_rows,           /**< number of rows in the data */
   double*               nosplit_rewards     /**< space for computed no split rewards */
   );

/** find best action and its associated reward for a set of units */
void find_best_action(
   const STRATEGY*       strategy,           /**< tree-building strategy */
   CONST_UNITS           units,              /**< units */
   const double*         data_y,             /**< gammas, data_y+(d*num_rows) points to values for reward d */
   int                   num_rows,           /**< number of units in full dataset */
   int                   num_cols_y,         /**< number of rewards/actions */
   WORKSPACE*            workspace,          /**< workspace */
   double*               best_reward,        /**< (pointer to) best reward */
   int*                  best_action         /**< (pointer to) best action */
   );


/** get number of units in a set of units */
int get_size(
   const STRATEGY*       strategy,           /**< tree-building strategy */
   CONST_UNITS           units               /**< units */
   );

/** create a full 'right' set to prepare for looking for depth=1 splits using a given covariate */
void shallow_initialise_units(
   const STRATEGY*       strategy,           /**< tree-building strategy */
   CONST_UNITS           units,              /**< input units */
   int                   p,                  /**< splitting covariate */
   int                   num_cols_x,         /**< number of covariates */
   WORKSPACE*            workspace,          /**< workspace */
   UNITS*                right_units         /**< pointer to output right units */
   );


/** create an empty 'left' set of units and full 'right' set to prepare for looking for splits using a given covariate */
void initialise_units(
   const STRATEGY*       strategy,           /**< tree-building strategy */
   CONST_UNITS           units,              /**< input units */
   int                   p,                  /**< splitting covariate */
   int                   depth,              /**< depth of associated node */
   int                   num_cols_x,         /**< number of covariates */
   WORKSPACE*            workspace,          /**< workspace */
   UNITS*                left_units,         /**< pointer to output left units */
   UNITS*                right_units         /**< pointer to output right units */
   );

/** given that left_set and right_set are sorted according to covariate p, find next split (if any)  and associated split value
 * if there is a split, then both left_units and right_units are updated accordingly.
 * @return 1 if there is a next split, otherwise 0
 */
int next_split(
   const STRATEGY*       strategy,           /**< tree-building strategy */
   CONST_UNITS           units,              /**< input set */
   UNITS                 left_units,         /**< units on the left */
   UNITS                 right_units,        /**< units on the right */ 
   int                   p,                  /**< covariate to split on */
   const double*         data_xp,            /**< values for covariate to split on */
   int                   num_cols_x,         /**< number of covariates */
   double*               splitval,           /**< (pointer to) found value to split on */
   ELEMENT**             elts,               /**< (pointer to) elements moved */
   int*                  nelts,              /**< (pointer to) number of elements moved */
   int                   splitcount          /**< number of previous splits */
   );

/** find units with same covariate value starting from a given index */
int next_shallow_split(
   const STRATEGY*       strategy,           /**< tree-building strategy */
   CONST_UNITS           units,              /**< units */
   int                   p,                  /**< covariate to split on */
   int                   start,              /**< starting index */
   const double*         data_xp,            /**< values for covariate to split on */
   double*               splitval,           /**< (pointer to) found value to split on */
   ELEMENT**             elts,               /**< (pointer to) the elements moved */
   int*                  nelts,              /**< (pointer to) number of elements moved */
   int                   splitcount          /**< number of previous splits */
   );

/** return array of elements and size of array for a set of units */
void elements(
   const STRATEGY*       strategy,           /**< tree-building strategy */
   CONST_UNITS           units,              /**< units */
   ELEMENT**             elts,               /**< (pointer to) elements */
   int*                  nelts               /**< (pointer to) number of elements */
   );

/** is a covariate binary ? */
int is_binary(
   const STRATEGY*       strategy,           /**< tree-building strategy */
   CONST_UNITS           units,              /**< units */
   int                   p                   /**< covariate */
   );

/** update left rewards using single split on given covariate, return whether there is,
 * in fact, a split
 */
int update_left_rewards_from_full(
   const STRATEGY*       strategy,           /**< tree-building strategy */
   CONST_UNITS           units,              /**< units */
   int                   p,                  /**< covariate to split on */
   double*               left_rewards,       /**< rewards for each action for a 'left' set of units */
   const double*         data_xp,            /**< values for covariate to split on */
   const double*         data_y,             /**< data_y[d*num_rows+elt] is the reward for action d for unit elt */
   int                   num_rows,           /**< number of units in full dataset */
   int                   num_cols_y,         /**< number of actions */
   double*               splitval            /**< (pointer to) found value to split on */
   );


#ifdef __cplusplus
}
#endif

#endif
