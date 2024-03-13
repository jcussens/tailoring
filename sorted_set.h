#ifdef __cplusplus
extern "C" {
#endif

#include "type_all.h"

/** for debugging only: check that a non-empty collection of sorted sets represent the same underlying set and that each is appropriately sorted */
int are_sorted_sets(
   CONST_UNITS           units,              /**< units */
   const double*         data_x,             /**< covariates, data_x+(j*num_rows) points to values for covariate j */
   int                   num_rows,           /**< number of units in full dataset */
   int                   num_cols_x          /**< number of covariates */
   );

UNITS make_initial_sorted_sets(
   const double*         data_x,             /**< covariates, data_x+(j*num_rows) points to values for covariate j */
   int                   num_rows,           /**< number of units in full dataset */
   int                   num_cols_x          /**< number of covariates */
   );


/** Determine whether a set is 'pure'.
 * A pure set is one where each unit has the same best action
 * @return 1 if the set is pure, else 0
 */
int is_pure(
   CONST_UNITS           units,              /**< units */
   const int*            best_actions        /**< best_actions[i] is the best action for unit i */
   );

/** get common size of sorted sets */
int get_size(
   CONST_UNITS           units               /**< units */
   );

/** find next splitting value ('splitval') for covariate p (if any) and move units from right to left so
 * that x[p] <= splitval for all units on left and x[p] > splitval on right. 
 * Return 1 if a split found, else 0
*/
int next_split(
   UNITS                 left_units,         /**< units on the left */
   UNITS                 right_units,        /**< units on the right */ 
   int                   p,                  /**< covariate to split on */
   const double*         data_xp,            /**< values for covariate to split on */
   int                   num_cols_x,         /**< number of covariates */
   WORKSPACE*            workspace,          /**< workspace */
   double*               splitval,           /**< (pointer to) found value to split on */
   int**                 elts,               /**< (pointer to) the elements moved */
   int*                  nelts               /**< (pointer to) number of elements moved */
   );

void find_nosplit_rewards(
   CONST_UNITS           units,              /**< units */
   int                   num_cols_y,         /**< number of actions */
   const double*         data_y,             /**< gammas, data_y+(d*num_rows) points to values for reward d */
   int                   num_rows,           /**< number of rows in the data */
   double*               nosplit_rewards     /**< space for computed no split rewards */
   );

/* find best action and its associated reward for a set of units */
void find_best_reward(
   CONST_UNITS           units,              /**< units */
   const double*         data_y,             /**< gammas, data_y+(d*num_rows) points to values for reward d */
   int                   num_rows,           /**< number of units in full dataset */
   int                   num_cols_y,         /**< number of rewards/actions */
   WORKSPACE*            workspace,          /**< workspace */
   double*               best_reward,        /**< (pointer to) best reward */
   int*                  best_action         /**< (pointer to) best action */
   );

/** initialise so that each left_sorted_set (for each covariate) is empty and each 
 *  right_sorted_set is a copy of the sorted set (for that covariate) 
*/
void initialise_sorted_sets(
   CONST_UNITS           sorted_sets,        /**< input units */
   int                   depth,              /**< depth of associated node */
   int                   num_cols_x,         /**< number of covariates */
   WORKSPACE*            workspace,          /**< workspace */
   UNITS*                left_units,         /**< pointer to output left units */
   UNITS*                right_units         /**< pointer to output right units */
   );

/* make a 'shallow' copy of source sorted sets */
UNITS shallow_copy_sorted_sets(
   CONST_UNITS           sources,            /**< source sorted sets */
   int                   nsets               /**< number of sources */
   );

void free_sorted_sets(
   SORTED_SET**          sorted_sets,        /**< sorted sets */
   int                   nsets               /**< number of sorted sets */
   );

void shallow_free_sorted_sets(
   SORTED_SET**          sorted_sets,        /**< sorted sets */
   int                   nsets               /**< number of sorted sets */
   );

/** find units with same covariate value starting from a given index */
int next_shallow_split(
   const SORTED_SET**    right_sorted_sets,  /**< sorted set */
   int                   p,                  /**< covariate to split on */
   int                   start,              /**< starting index */
   const double*         data_xp,            /**< values for covariate to split on */
   double*               splitval,           /**< (pointer to) found value to split on */
   int**                 elts,               /**< (pointer to) the elements moved */
   int*                  nelts               /**< (pointer to) number of elements moved */
   );

#ifdef __cplusplus
}
#endif
