/** @file strategy.c
 *  @brief Solving strategy
 *  @author James Cussens
 */

#include "strategy.h"
#include "stdlib.h"

enum datatype
{
   SORTED_SET = 0,
   SIMPLE_SET = 1
};
typedef enum datatype DATATYPE;
   
struct strategy
{
   DATATYPE              datatype;           /**< type for data (currently either
                                                sorted sets (policytree style) or simple set */
};

/** return an uninitialised strategy */
STRATEGY* get_unint_strategy(
   void
   )
{
   return (STRATEGY*) malloc(sizeof(STRATEGY));
}
   
/** are we using sorted sets for datasets? */
int using_sorted_sets(
   const STRATEGY*       strategy            /**< solving strategy */
   )
{
   return (int) (strategy->datatype == SORTED_SET);
}

/** decide to use sorted sets */
void use_sorted_sets(
   STRATEGY*             strategy            /**< solving strategy */
   )
{
   strategy->datatype = SORTED_SET;
}

/** decide to use simple sets */
void use_simple_sets(
   STRATEGY*             strategy            /**< solving strategy */
   )
{
   strategy->datatype = SIMPLE_SET;
}


