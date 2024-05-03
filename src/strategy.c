/** @file strategy.c
 *  @brief Solving strategy
 *  @author James Cussens
 */

#include "strategy.h"
#include "simple_set.h"
#include "stdlib.h"

enum datatype
{
   SORTED_SET_TYPE = 0,
   SIMPLE_SET_TYPE = 1,
   UNDECIDED_TYPE = 2
};
typedef enum datatype DATATYPE;
   
struct strategy
{
   DATATYPE              datatype;           /**< type for data (currently either
                                                sorted sets (policytree style) or simple set */
   int                   find_reward_ub;     /**< whether to compute an upper bound on reward for each (sub-) dataset */
   int                   find_dummy_split_reward;     /**< whether to compute rewards for dummy splits */
   int                   use_last_rewards;   /**< whether to use last rewards */
   int                   use_cutoffs;        /**< whether to use cutoffs */
   int                   use_cache;          /**< whether to use the cache */
};

/** return an uninitialised strategy */
STRATEGY* get_unint_strategy(
   void
   )
{
   STRATEGY* strategy = (STRATEGY*) malloc(sizeof(STRATEGY));
   strategy->datatype = UNDECIDED_TYPE;
   strategy->find_reward_ub = -1;
   strategy->find_dummy_split_reward = -1;
   strategy->use_last_rewards = -1;
   strategy->use_cutoffs = -1;
   strategy->use_cache = -1;
   return strategy;
}
   
/** are we using sorted sets for datasets? */
int using_sorted_sets(
   const STRATEGY*       strategy            /**< solving strategy */
   )
{
   return (int) (strategy->datatype == SORTED_SET_TYPE);
}

/** decide to use sorted sets */
void use_sorted_sets(
   STRATEGY*             strategy            /**< solving strategy */
   )
{
   strategy->datatype = SORTED_SET_TYPE;
}

/** decide to use simple sets */
void use_simple_sets(
   STRATEGY*             strategy            /**< solving strategy */
   )
{
   strategy->datatype = SIMPLE_SET_TYPE;
}

void decide_datatype(
   STRATEGY*             strategy,           /**< solving strategy */
   const double*         data_x,             /**< covariates, data_x+(j*num_rows) points to values for covariate j */
   int                   num_rows,           /**< number of units in full dataset */
   int                   num_cols_x          /**< number of covariates */
   )
{
   int p;
   SIMPLE_SET* simple_set = simple_set_make_units(data_x, num_rows, num_cols_x);
   int nfewkeyvals = 0;
   const int threshold = 30;
   
   for( p = 0; p < num_cols_x; p++ )
      if( nkeyvals((const SIMPLE_SET*) simple_set, p) < threshold )
         nfewkeyvals++;

   /* if most covariates don't have too many distinct values go for simple set,
    * else sorted sets
    */
   if( nfewkeyvals > num_cols_x/2 )
      strategy->datatype = SIMPLE_SET_TYPE;
   else
      strategy->datatype = SORTED_SET_TYPE;

   simple_set_free_units(simple_set, num_cols_x);
}

/** return whether we wish to compute an upper bound on reward for each (sub-) dataset */
int find_reward_ub(
   const STRATEGY*       strategy            /**< solving strategy */
   )
{
   return strategy->find_reward_ub;
}

/** set whether we wish to compute an upper bound on reward for each (sub-) dataset */
void set_find_reward_ub(
   STRATEGY*             strategy,           /**< solving strategy */
   int                   val                 /**< 0 for no, 1 for yes */
   )
{
   strategy->find_reward_ub = val;
}

/** return whether we wish to find dummy split rewards */
int find_dummy_split_reward(
   const STRATEGY*       strategy            /**< solving strategy */
   )
{
   return strategy->find_dummy_split_reward;
}

/** set whether we wish to find dummy split rewards */
void set_find_dummy_split_reward(
   STRATEGY*             strategy,           /**< solving strategy */
   int                   val                 /**< 0 for no, 1 for yes */
   )
{
   strategy->find_dummy_split_reward = val;
}

/** return whether we wish to use last rewards rewards */
int use_last_rewards(
   const STRATEGY*       strategy            /**< solving strategy */
   )
{
   return strategy->use_last_rewards;
}

/** set whether we wish to use last rewards */
void set_use_last_rewards(
   STRATEGY*             strategy,           /**< solving strategy */
   int                   val                 /**< 0 for no, 1 for yes */
   )
{
   strategy->use_last_rewards = val;
}

/** return whether we wish to use cutoffs */
int use_cutoffs(
   const STRATEGY*       strategy            /**< solving strategy */
   )
{
   return strategy->use_cutoffs;
}

/** set whether we wish to use cutoffs */
void set_use_cutoffs(
   STRATEGY*             strategy,           /**< solving strategy */
   int                   val                 /**< 0 for no, 1 for yes */
   )
{
   strategy->use_cutoffs = val;
}

/** return whether we wish to use the cache */
int use_cache(
   const STRATEGY*       strategy            /**< solving strategy */
   )
{
   return strategy->use_cache;
}

/** set whether we wish to use the cache */
void set_use_cache(
   STRATEGY*             strategy,           /**< solving strategy */
   int                   val                 /**< 0 for no, 1 for yes */
   )
{
   strategy->use_cache = val;
}

