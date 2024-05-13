/** @file sorted_set.h
 *  @brief Implements a set of units with an array for each covariate
 *  @author James Cussens
 */
#ifdef __cplusplus
extern "C" {
#endif

#include "type_all.h"

struct sorted_set;
typedef struct sorted_set SORTED_SET;  /**< Sorted sets */

/** make unit set from covariate data */
SORTED_SET** sorted_set_make_units(
   const double*         data_x,             /**< covariates, data_x+(j*num_rows) points to values for covariate j */
   int                   num_rows,           /**< number of units in full dataset */
   int                   num_cols_x          /**< number of covariates */
   );

/** free (space for) units */
void sorted_set_free_units(
   SORTED_SET**          units,              /**< units */
   int                   num_cols_x          /**< number of covariates */
   );

/** make a 'shallow' copy of units, these copies should be freed by calling shallow_free_units */
SORTED_SET** sorted_set_shallow_copy_units(
   const SORTED_SET**    units,              /**< source units */
   int                   num_cols_x          /**< number of covariates */
   );

/** free (space for) units previously created as shallow copies */
void sorted_set_shallow_free_units(
   SORTED_SET**          sorted_sets,        /**< units */
   int                   num_cols_x          /**< number of covariates */
   );

/** for debugging only: check that a non-empty collection of sorted sets represent the same underlying set and that each is appropriately sorted */
int sorted_set_units_ok(
   const SORTED_SET**    units,              /**< units */
   int                   p,                  /**< if !=-1 then units must be ready for splitting on covariate p */
   const double*         data_x,             /**< covariates, data_x+(j*num_rows) points to values for covariate j */
   int                   num_rows,           /**< number of units in full dataset */
   int                   num_cols_x          /**< number of covariates */
   );

/** get the reward for a set if all units in the set were assigned their best action *
 * @return the reward for a set if all units in the set were assigned their best action *
 */
double sorted_set_get_reward_ub(
   const SORTED_SET**    units,              /**< set */
   const double*         data_y,             /**< gammas, data_y+(d*num_rows) points to values for reward d */
   int                   num_rows,           /**< number of units in full dataset */
   const int*            best_actions        /**< best_actions[i] is the best action for unit i */
   );


/** Determine whether a set is 'pure'.
 * A pure set is one where each unit has the same best action
 * @return 1 if the set is pure, else 0
 */
int sorted_set_is_pure(
   const SORTED_SET**    units,              /**< units */
   const int*            best_actions        /**< best_actions[i] is the best action for unit i */
   );

/** for each action find (total) reward if that action applied to each unit in a set of units */
void sorted_set_find_nosplit_rewards(
   const SORTED_SET**    units,              /**< units */
   int                   num_cols_y,         /**< number of actions */
   const double*         data_y,             /**< gammas, data_y+(d*num_rows) points to values for reward d */
   int                   num_rows,           /**< number of rows in the data */
   double*               nosplit_rewards     /**< space for computed no split rewards */
   );

/** find best action and its associated reward for a set of units */
void sorted_set_find_best_action(
   const SORTED_SET**    units,              /**< units */
   const double*         data_y,             /**< gammas, data_y+(d*num_rows) points to values for reward d */
   int                   num_rows,           /**< number of units in full dataset */
   int                   num_cols_y,         /**< number of rewards/actions */
   WORKSPACE*            workspace,          /**< workspace */
   double*               best_reward,        /**< (pointer to) best reward */
   int*                  best_action         /**< (pointer to) best action */
   );


/** get number of units in a set of units */
int sorted_set_get_size(
   const SORTED_SET**    units               /**< units */
   );

/** create a full 'right' set to prepare for looking for depth=1 splits using a given covariate */
void sorted_set_shallow_initialise_units(
   const SORTED_SET**    units,              /**< input units */
   int                   p,                  /**< splitting covariate */
   int                   num_cols_x,         /**< number of covariates */
   WORKSPACE*            workspace,          /**< workspace */
   SORTED_SET***         right_units         /**< pointer to output right units */
   );


/** create an empty 'left' set of units and full 'right' set to prepare for looking for splits using a given covariate */
void sorted_set_initialise_units(
   const SORTED_SET**    units,              /**< input units */
   int                   p,                  /**< splitting covariate */
   int                   depth,              /**< depth of associated node */
   int                   num_cols_x,         /**< number of covariates */
   WORKSPACE*            workspace,          /**< workspace */
   SORTED_SET***         left_units,         /**< pointer to output left units */
   SORTED_SET***         right_units         /**< pointer to output right units */
   );

/** given that left_set and right_set are sorted according to covariate p, find next split (if any)  and associated split value
 * if there is a split, then both left_units and right_units are updated accordingly.
 * @return 1 if there is a next split, otherwise 0
 */
int sorted_set_next_split(
   const STRATEGY*       strategy,           /**< tree-building strategy */
   SORTED_SET**          left_units,         /**< units on the left */
   SORTED_SET**          right_units,        /**< units on the right */ 
   int                   p,                  /**< covariate to split on */
   const double*         data_xp,            /**< values for covariate to split on */
   int                   num_cols_x,         /**< number of covariates */
   double*               splitval,           /**< (pointer to) found value to split on */
   ELEMENT**             elts,               /**< (pointer to) elements moved */
   int*                  nelts               /**< (pointer to) number of elements moved */
   );

/** find units with same covariate value starting from a given index */
int sorted_set_next_shallow_split(
   const SORTED_SET**    right_units,        /**< units */
   int                   p,                  /**< covariate to split on */
   int                   start,              /**< starting index */
   const double*         data_xp,            /**< values for covariate to split on */
   double*               splitval,           /**< (pointer to) found value to split on */
   ELEMENT**             elts,               /**< (pointer to) the elements moved */
   int*                  nelts               /**< (pointer to) number of elements moved */
   );

/** return array of elements and size of array for a set of units */
void sorted_set_elements(
   const SORTED_SET**    units,              /**< units */
   ELEMENT**             elts,               /**< (pointer to) elements */
   int*                  nelts               /**< (pointer to) number of elements */
   );


#ifdef __cplusplus
}
#endif
