/** @file cache.h
 *  @brief Functions for the cache
 *  @author James Cussens
 */
#ifndef __CACHE_H__
#define __CACHE_H__

#include "type_all.h"

#ifdef __cplusplus
extern "C" {
#endif

struct cache;
typedef struct cache CACHE;

/** search the cache for an optimal tree of given depth for the given set 
 * @return 1 if tree is in the cache, else 0
 */
int search_cache(
   const CACHE*          cache,              /**< cache */
   int                   nelts,              /**< number of elts in set */
   const ELEMENT*        elts,               /**< the set */
   int                   depth,              /**< depth of tree */
   NODE*                 tree                /**< if optimal tree is in cache, tree is set to be that optimal tree */
   );

/** add an optimal tree of given depth for the given set to the cache */
void add_to_cache(
   CACHE*                cache,              /**< cache */
   int                   nelts,              /**< number of elts in set */
   const ELEMENT*        elts,               /**< the set */
   int                   depth,              /**< depth of tree */
   const NODE*           tree                /**< optimal tree for set */
   );

/** make (an empty) cache 
 * @return an empty cache
 */
CACHE* make_cache(
   int                   num_rows            /**< number of units */
   );

/** free the cache */
void free_cache(
   CACHE*                cache               /**< cache */
   );


#ifdef __cplusplus
}
#endif

#endif
