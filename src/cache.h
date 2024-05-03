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

/** get the number of ints required to represent a set of units */
int getnints(
   const CACHE*          cache               /**< cache */
   );

/** search cache for an optimal tree */
int search_cache(
   const CACHE*          cache,              /**< cache */
   int                   nelts,              /**< number of elts in set */
   const ELEMENT*        elts,               /**< the set */
   int                   depth,              /**< depth of tree */
   NODE*                 tree                /**< if in cache, tree is set to optimal tree for set */
   );

/** add an optimal tree to cache */
void add_to_cache(
   CACHE*                cache,              /**< cache */
   int                   nelts,              /**< number of elts in set */
   const ELEMENT*        elts,               /**< the set */
   int                   depth,              /**< depth of tree */
   const NODE*           tree                /**< optimal tree for set */
   );

/** make cache */
CACHE* make_cache(
   int                   num_rows            /**< number of units */
   );

/** free cache */
void free_cache(
   CACHE*                cache               /**< cache */
   );


#ifdef __cplusplus
}
#endif

#endif
