/** @file units.c
 *  @brief Dispatches to appropriate data type for units
 *  @author James Cussens
 */
#ifdef __cplusplus
extern "C" {
#endif

#include "type_all.h"
#include "units.h"
#include "sorted_set.h"
#include "simple_set.h"
#include "strategy.h"

#define IFUSESORTED if( using_sorted_sets(strategy) )  
#define DISPATCH_NOUNITS(F,...) IFUSESORTED  \
   { return sorted_set_ ## F(__VA_ARGS__); } \
   else {return simple_set_ ## F(__VA_ARGS__); }
#define DISPATCH_UNITS(F,...) IFUSESORTED  \
   { return sorted_set_ ## F((SORTED_SET**) units, __VA_ARGS__); } \
   else {return simple_set_ ## F((SIMPLE_SET*) units, __VA_ARGS__); }
#define DISPATCH_UNITS_ONLY(F,...) IFUSESORTED  \
   { return sorted_set_ ## F((SORTED_SET**) units); } \
   else {return simple_set_ ## F((SIMPLE_SET*) units); }
#define DISPATCH_CONSTUNITS(F,...) IFUSESORTED  \
   { return sorted_set_ ## F((const SORTED_SET**) units, __VA_ARGS__); } \
   else {return simple_set_ ## F((const SIMPLE_SET*) units, __VA_ARGS__); }
#define DISPATCH_CONSTUNITS_ONLY(F,...) IFUSESORTED  \
   { return sorted_set_ ## F((const SORTED_SET**) units); } \
   else {return simple_set_ ## F((const SIMPLE_SET*) units); }
#define DISPATCH(F,A,B) IFUSESORTED \
   { return sorted_set_ ## F A; } \
   else {return simple_set_ ## F B; } 
#define DISPATCH_NOUNITS_CAST(R,F,...) IFUSESORTED  \
   { return (R) sorted_set_ ## F(__VA_ARGS__); } \
   else {return (R) simple_set_ ## F(__VA_ARGS__); }
#define DISPATCH_UNITS_CAST(R,F,...) IFUSESORTED                   \
   { return (R) sorted_set_ ## F((SORTED_SET**) units, __VA_ARGS__); }   \
   else {return (R) simple_set_ ## F((SIMPLE_SET*) units, __VA_ARGS__); }
#define DISPATCH_CONSTUNITS_CAST(R,F,...) IFUSESORTED                   \
   { return (R) sorted_set_ ## F((const SORTED_SET**) units, __VA_ARGS__); }   \
   else {return (R) simple_set_ ## F((const SIMPLE_SET*) units, __VA_ARGS__); } 
#define DISPATCH_CAST(R,F,A,B) IFUSESORTED \
   { return (R) sorted_set_ ## F A; } \
   else {return (R) simple_set_ ## F B; } 


/** make unit set from covariate data */
UNITS make_units(
   const STRATEGY*       strategy,           /**< tree-building strategy */
   const double*         data_x,             /**< covariates, data_x+(j*num_rows) points to values for covariate j */
   int                   num_rows,           /**< number of units in full dataset */
   int                   num_cols_x          /**< number of covariates */
   )
{
   DISPATCH_NOUNITS_CAST(UNITS, make_units, data_x, num_rows, num_cols_x)
}

/** free (space for) units */
void free_units(
   const STRATEGY*       strategy,           /**< tree-building strategy */
   UNITS                 units,              /**< units */
   int                   num_cols_x          /**< number of covariates */
   )
{
   DISPATCH_UNITS(free_units, num_cols_x)
}

/** make a 'shallow' copy of units, these copies should be freed by calling shallow_free_units */
UNITS shallow_copy_units(
   const STRATEGY*       strategy,           /**< tree-building strategy */
   CONST_UNITS           units,              /**< source units */
   int                   num_cols_x          /**< number of covariates */
   )
{
   DISPATCH_CONSTUNITS_CAST(UNITS, shallow_copy_units, num_cols_x)
}

/** free (space for) units previously created as shallow copies */
void shallow_free_units(
   const STRATEGY*       strategy,           /**< tree-building strategy */
   UNITS                 units,              /**< units */
   int                   num_cols_x          /**< number of covariates */
   )
{
   DISPATCH_UNITS(shallow_free_units, num_cols_x)
}

/** for debugging only: check that a non-empty collection of sorted sets represent the same underlying set and that each is appropriately sorted 
 * this is very slow, so by default no check is done
 */
int units_ok(
   const STRATEGY*       strategy,           /**< tree-building strategy */
   CONST_UNITS           units,              /**< units */
   int                   p,                  /**< if !=-1 then units must be ready for splitting on covariate p */
   const double*         data_x,             /**< covariates, data_x+(j*num_rows) points to values for covariate j */
   int                   num_rows,           /**< number of units in full dataset */
   int                   num_cols_x          /**< number of covariates */
   )
{
#ifdef CHECK_UNITSOK
   DISPATCH_CONSTUNITS(units_ok, p, data_x, num_rows, num_cols_x)
#else
   return 1;
#endif
}

/** get the reward for a set if all units in the set were assigned their best action *
 * @return the reward for a set if all units in the set were assigned their best action *
 */
double get_reward_ub(
   const STRATEGY*       strategy,           /**< tree-building strategy */
   CONST_UNITS           units,              /**< set */
   const double*         data_y,             /**< gammas, data_y+(d*num_rows) points to values for reward d */
   int                   num_rows,           /**< number of units in full dataset */
   const int*            best_actions        /**< best_actions[i] is the best action for unit i */
   )
{
   DISPATCH_CONSTUNITS(get_reward_ub, data_y, num_rows, best_actions)
}

/** Determine whether a set is 'pure'.
 * A pure set is one where each unit has the same best action
 * @return 1 if the set is pure, else 0
 */
int is_pure(
   const STRATEGY*       strategy,           /**< tree-building strategy */
   CONST_UNITS           units,              /**< units */
   const int*            best_actions        /**< best_actions[i] is the best action for unit i */
   )
{
   DISPATCH_CONSTUNITS(is_pure, best_actions)
}

/** for each action find (total) reward if that action applied to each unit in a set of units */
void find_nosplit_rewards(
   const STRATEGY*       strategy,           /**< tree-building strategy */
   CONST_UNITS           units,              /**< units */
   int                   num_cols_y,         /**< number of actions */
   const double*         data_y,             /**< gammas, data_y+(d*num_rows) points to values for reward d */
   int                   num_rows,           /**< number of rows in the data */
   double*               nosplit_rewards     /**< space for computed no split rewards */
   )
{
   DISPATCH_CONSTUNITS(find_nosplit_rewards, num_cols_y, data_y, num_rows, nosplit_rewards)
}

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
   )
{
   DISPATCH_CONSTUNITS(find_best_action, data_y, num_rows, num_cols_y, workspace, best_reward, best_action)
}

/** get number of units in a set of units */
int get_size(
   const STRATEGY*       strategy,           /**< tree-building strategy */
   CONST_UNITS           units               /**< units */
   )
{
   DISPATCH_CONSTUNITS_ONLY(get_size)
}

/** create a full 'right' set to prepare for looking for depth=1 splits using a given covariate */
void shallow_initialise_units(
   const STRATEGY*       strategy,           /**< tree-building strategy */
   CONST_UNITS           units,              /**< input units */
   int                   p,                  /**< splitting covariate */
   int                   num_cols_x,         /**< number of covariates */
   WORKSPACE*            workspace,          /**< workspace */
   UNITS*                right_units         /**< pointer to output right units */
   )
{
   DISPATCH(shallow_initialise_units, ((const SORTED_SET**) units, p, num_cols_x, workspace, (SORTED_SET***) right_units),
      ((const SIMPLE_SET*) units, p, num_cols_x, workspace, (SIMPLE_SET**) right_units))
}


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
   )
{
   DISPATCH(initialise_units, ((const SORTED_SET**) units, p, depth, num_cols_x, workspace, (SORTED_SET***) left_units, (SORTED_SET***) right_units),
      ((const SIMPLE_SET*) units, p, depth, num_cols_x, workspace, (SIMPLE_SET**) left_units, (SIMPLE_SET**) right_units))
}

/** given that left_set and right_set are sorted according to covariate p, find next split (if any)  and associated split value
 * if there is a split, then both left_units and right_units are updated accordingly.
 * @return 1 if there is a next split, otherwise 0
 */
int next_split(
   const STRATEGY*       strategy,           /**< tree-building strategy */
   UNITS                 left_units,         /**< units on the left */
   UNITS                 right_units,        /**< units on the right */ 
   int                   p,                  /**< covariate to split on */
   const double*         data_xp,            /**< values for covariate to split on */
   int                   num_cols_x,         /**< number of covariates */
   double*               splitval,           /**< (pointer to) found value to split on */
   ELEMENT**             elts,               /**< (pointer to) elements moved */
   int*                  nelts               /**< (pointer to) number of elements moved */
   )
{
   DISPATCH(next_split, ((SORTED_SET**) left_units, (SORTED_SET**) right_units, p, data_xp, num_cols_x, splitval, elts, nelts),
      ((SIMPLE_SET*) left_units, (SIMPLE_SET*) right_units, p, data_xp, num_cols_x, splitval, elts, nelts))
}

/** find units with same covariate value starting from a given index */
int next_shallow_split(
   const STRATEGY*       strategy,           /**< tree-building strategy */
   CONST_UNITS           units,              /**< units */
   int                   p,                  /**< covariate to split on */
   int                   start,              /**< starting index */
   const double*         data_xp,            /**< values for covariate to split on */
   double*               splitval,           /**< (pointer to) found value to split on */
   ELEMENT**             elts,               /**< (pointer to) the elements moved */
   int*                  nelts               /**< (pointer to) number of elements moved */
   )
{
   DISPATCH_CONSTUNITS(next_shallow_split, p, start, data_xp, splitval, elts, nelts)
}

/** return array of elements and size of array for a set of units */
void elements(
   const STRATEGY*       strategy,           /**< tree-building strategy */
   CONST_UNITS           units,              /**< units */
   ELEMENT**             elts,               /**< (pointer to) elements */
   int*                  nelts               /**< (pointer to) number of elements */
   )
{
   DISPATCH_CONSTUNITS(elements, elts, nelts)
}


#ifdef __cplusplus
}
#endif
