#ifdef __cplusplus
extern "C" {
#endif

#include "workspace.h"

struct sorted_set;
typedef struct sorted_set SORTED_SET;  /**< Trie structure for storing scored parent sets */

SORTED_SET** make_initial_sorted_sets(
   const double*         data_x,             /**< covariates, data_x+(j*num_rows) points to values for covariate j */
   int                   num_rows,           /**< number of units in full dataset */
   int                   num_cols_x          /**< number of covariates */
   );


/** Determine whether a set is 'pure'.
 * A pure set is one where each unit has the same best action
 * @return 1 if the set is pure, else 0
 */
int is_pure(
   const SORTED_SET**    sorted_sets,        /**< sorted sets, representing a common unsorted set */
   const int*            best_actions        /**< best_actions[i] is the best action for unit i */
   );

/** get common size of sorted sets */
int get_size(
   const SORTED_SET**    sorted_sets         /**< sorted sets */
   );

/** find next splitting value ('splitval') for covariate p (if any) and move units from right sorted set for p to left so
 * that x[p] <= splitval for all units on left and x[p] > splitval on right. 
 * Return 1 if a split found, else 0
*/
int next_split(
   SORTED_SET**          left_sorted_sets,   /**< left sorted sets, one for each covariate */
   SORTED_SET**          right_sorted_sets,  /**< right sorted sets, one for each covariate */ 
   int                   p,                  /**< covariate to split on */
   const double*         data_xp,            /**< values for covariate to split on */
   int                   num_cols_x,         /**< number of covariates */
   WORKSPACE*            workspace,          /**< workspace */
   double*               splitval,           /**< (pointer to) found value to split on */
   int**                 elts,               /**< (pointer to) the elements moved */
   int*                  nelts               /**< (pointer to) number of elements moved */
   )


/* find best action and its associated reward for a set of units */
void find_best_reward(
   const SORTED_SETS**   sorted_sets,        /**< sorted sets */
   const double*         data_y,             /**< gammas, data_y+(d*num_rows) points to values for reward d */
   int                   num_rows,           /**< number of units in full dataset */
   int                   num_cols_y,         /**< number of rewards/actions */
   WORKSPACE*            workspace,          /**< workspace */
   double*               best_reward,        /**< (pointer to) best reward */
   int*                  best_action,        /**< (pointer to) best action */
   )


#ifdef __cplusplus
}
#endif
