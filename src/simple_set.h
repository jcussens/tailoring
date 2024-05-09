/** @file simple_set.h
 *  @brief Prototypes for functions manipulating sets of units, using 'simple_set' approach
 *  @author James Cussens
 */
#ifdef __cplusplus
extern "C" {
#endif

#include "type_all.h"

struct simple_set;
typedef struct simple_set SIMPLE_SET;  /**< Simple sets */

/** make unit set from covariate data */
SIMPLE_SET* simple_set_make_units(
   const double*         data_x,             /**< covariates, data_x+(j*num_rows) points to values for covariate j */
   int                   num_rows,           /**< number of units in full dataset */
   int                   num_cols_x          /**< number of covariates */
   );

/** free (space for) units */
void simple_set_free_units(
   SIMPLE_SET*           units,              /**< units */
   int                   num_cols_x          /**< number of covariates */
   );

/** make a 'shallow' copy of units, these copies should be freed by calling shallow_free_units */
SIMPLE_SET* simple_set_shallow_copy_units(
   const SIMPLE_SET*     units,              /**< source units */
   int                   num_cols_x          /**< number of covariates */
   );

/** free (space for) units previously created as shallow copies */
void simple_set_shallow_free_units(
   SIMPLE_SET*           units,              /**< units */
   int                   num_cols_x          /**< number of covariates */
   );

/** for debugging only: check that there are no problems */
int simple_set_units_ok(
   const SIMPLE_SET*     units,              /**< units */
   int                   p,                  /**< if !=-1 then units must be ready for splitting on covariate p */
   const double*         data_x,             /**< covariates, data_x+(j*num_rows) points to values for covariate j */
   int                   num_rows,           /**< number of units in full dataset */
   int                   num_cols_x          /**< number of covariates */
   );

/** get the reward for a set if all units in the set were assigned their best action *
 * @return the reward for a set if all units in the set were assigned their best action *
 */
double simple_set_get_reward_ub(
   const SIMPLE_SET*     units,              /**< set */
   const double*         data_y,             /**< gammas, data_y+(d*num_rows) points to values for reward d */
   int                   num_rows,           /**< number of units in full dataset */
   const int*            best_actions        /**< best_actions[i] is the best action for unit i */
   );


/** Determine whether a set is 'pure'.
 * A pure set is one where each unit has the same best action
 * @return 1 if the set is pure, else 0
 */
int simple_set_is_pure(
   const SIMPLE_SET*     units,              /**< units */
   const int*            best_actions        /**< best_actions[i] is the best action for unit i */
   );

/** for each action find (total) reward if that action applied to each unit in a set of units */
void simple_set_find_nosplit_rewards(
   const SIMPLE_SET*     units,              /**< units */
   int                   num_cols_y,         /**< number of actions */
   const double*         data_y,             /**< gammas, data_y+(d*num_rows) points to values for reward d */
   int                   num_rows,           /**< number of rows in the data */
   double*               nosplit_rewards     /**< space for computed no split rewards */
   );

/** find best action and its associated reward for a set of units */
void simple_set_find_best_action(
   const SIMPLE_SET*     units,              /**< units */
   const double*         data_y,             /**< gammas, data_y+(d*num_rows) points to values for reward d */
   int                   num_rows,           /**< number of units in full dataset */
   int                   num_cols_y,         /**< number of rewards/actions */
   WORKSPACE*            workspace,          /**< workspace */
   double*               best_reward,        /**< (pointer to) best reward */
   int*                  best_action         /**< (pointer to) best action */
   );


/** get number of units in a set of units */
int simple_set_get_size(
   const SIMPLE_SET*     units               /**< units */
   );

/** create a full 'right' set to prepare for looking for depth=1 splits using a given covariate */
void simple_set_shallow_initialise_units(
   const STRATEGY*       strategy,           /**< tree-building strategy */
   const SIMPLE_SET*     units,              /**< input units */
   int                   p,                  /**< splitting covariate */
   int                   num_cols_x,         /**< number of covariates */
   WORKSPACE*            workspace,          /**< workspace */
   SIMPLE_SET**          right_units         /**< pointer to output right units */
   );


/** create an empty 'left' set of units and full 'right' set to prepare for looking for splits using a given covariate */
void simple_set_initialise_units(
   const STRATEGY*       strategy,           /**< tree-building strategy */
   const SIMPLE_SET*     units,              /**< input units */
   int                   p,                  /**< splitting covariate */
   int                   depth,              /**< depth of associated node */
   int                   num_cols_x,         /**< number of covariates */
   WORKSPACE*            workspace,          /**< workspace */
   SIMPLE_SET**          left_units,         /**< pointer to output left units */
   SIMPLE_SET**          right_units         /**< pointer to output right units */
   );

/** given that left_set and right_set are sorted according to covariate p, find next split (if any)  and associated split value
 * if there is a split, then both left_units and right_units are updated accordingly.
 * If splitting on a binary variable then a different approach is taken
 * @return 1 if there is a next split, otherwise 0
 */
int simple_set_next_split(
   const STRATEGY*       strategy,           /**< tree-building strategy */
   const SIMPLE_SET*     simple_set,         /**< input set */
   SIMPLE_SET*           left_units,         /**< units on the left */
   SIMPLE_SET*           right_units,        /**< units on the right */ 
   int                   p,                  /**< covariate to split on */
   const double*         data_xp,            /**< values for covariate to split on */
   int                   num_cols_x,         /**< number of covariates */
   double*               splitval,           /**< (pointer to) found value to split on */
   ELEMENT**             elts,               /**< (pointer to) elements moved */
   int*                  nelts,              /**< (pointer to) number of elements moved */
   int                   splitcount          /**< number of previous splits */
   );

/** find units with same covariate value starting from a given index */
int simple_set_next_shallow_split(
   const STRATEGY*       strategy,           /**< tree-building strategy */
   const SIMPLE_SET*     right_units,        /**< units */
   int                   p,                  /**< covariate to split on */
   int                   start,              /**< starting index */
   const double*         data_xp,            /**< values for covariate to split on */
   double*               splitval,           /**< (pointer to) found value to split on */
   ELEMENT**             elts,               /**< (pointer to) the elements moved */
   int*                  nelts,               /**< (pointer to) number of elements moved */
   int                   splitcount          /**< number of previous splits */
   );

/** return array of elements and size of array for a set of units */
void simple_set_elements(
   const SIMPLE_SET*     units,              /**< units */
   ELEMENT**             elts,               /**< (pointer to) elements */
   int*                  nelts               /**< (pointer to) number of elements */
   );

/** get number of distinct values of a covariate */
int nkeyvals(
   const SIMPLE_SET*     simple_set,         /**< simple set */
   int                   i                   /**< covariate */
   );

int simple_set_update_left_rewards_from_full(
   const SIMPLE_SET*     simple_set,         /**< units */
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
