#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "cache.h"
#include "tree.h"
#include "assert.h"

#define NSLOTS 100

struct entry
{
   int                   nelts;              /**< number of elts in set */
   ELEMENT*              elts;               /**< the set */
   int                   depth;              /**< depth of tree */
   NODE*                 tree;               /**< optimal tree for set */
};
typedef struct entry ENTRY;

struct cache
{
   ENTRY***              slots;              /** slots */
   int*                  nentries;           /** number of entries in each slot */
};

/** get slot for given set and depth */
static
int get_slot(
   int                   nelts,              /**< number of elts in set */
   const ELEMENT*        elts,               /**< the set */
   int                   depth               /**< depth of tree */
   )
{
   return (nelts+depth) % NSLOTS;
}

/** is this the entry for give set and depth? */
static
int match(
   int                   nelts,              /**< number of elts in set */
   const ELEMENT*        elts,               /**< the set */
   int                   depth,              /**< depth of tree */
   const ENTRY*          entry               /**< entry */
   )
{
   int i;
   
   if( depth != entry->depth || nelts != entry->nelts )
      return 0;

   for( i = 0; i < nelts; i++ )
      if( elts[i] != entry->elts[i] )
         return 0;

   return 1;
}

/** create an entry */
static
ENTRY* make_entry(
   int                   nelts,              /**< number of elts in set */
   const ELEMENT*        elts,               /**< the set */
   int                   depth,              /**< depth of tree */
   const NODE*           tree                /**< optimal tree for set */
   )
{
   ENTRY* entry = (ENTRY*) malloc(sizeof(ENTRY));
   NODE* tree_cp = make_tree(depth);
   ELEMENT* elts_cp = (ELEMENT*) malloc(nelts*sizeof(ELEMENT));
   
   tree_copy(tree, tree_cp);
   memcpy(elts_cp, elts, nelts*sizeof(ELEMENT));

   entry->nelts = nelts;
   entry->elts = elts_cp;
   entry->depth = depth;
   entry->tree = tree_cp;

   return entry;
}

/** free an entry */
static
void free_entry(
   ENTRY* entry
   )
{
   free(entry->elts);
   tree_free(entry->tree);
   free(entry);
}

/** search cache for an optimal tree */
int search_cache(
   const CACHE*          cache,              /**< cache */
   int                   nelts,              /**< number of elts in set */
   const ELEMENT*        elts,               /**< the set */
   int                   depth,              /**< depth of tree */
   NODE*                 tree                /**< if in cache, tree is set to optimal tree for set */
   )
{
   int slot = get_slot(nelts, elts, depth);
   int i;
   for( i = 0; i < cache->nentries[slot]; i++ )
   {
      if( match(nelts, elts, depth, (const ENTRY*) cache->slots[slot][i]) )
      {
         tree_copy(cache->slots[slot][i]->tree,tree);
         return 1;
      }
   }
   return 0;
}


/** add an optimal tree to cache */
void add_to_cache(
   CACHE*                cache,              /**< cache */
   int                   nelts,              /**< number of elts in set */
   const ELEMENT*        elts,               /**< the set */
   int                   depth,              /**< depth of tree */
   const NODE*           tree                /**< optimal tree for set */
)
{
   int slot = get_slot(nelts, elts, depth);
   ENTRY* entry = make_entry(nelts, elts, depth, tree);
   int idx = cache->nentries[slot];

   assert(cache != NULL);
   assert(nelts >= 0);
   assert(elts != NULL);
   assert(depth >= 0);
   assert(tree != NULL);

   assert(slot >= 0);
   assert(slot < NSLOTS);
   assert(entry != NULL);
   assert(idx >= 0);

   if( idx == 0 )
      cache->slots[slot] = (ENTRY**) malloc(sizeof(ENTRY*));
   else
      cache->slots[slot] = (ENTRY**) realloc(cache->slots[slot], (idx+1)*sizeof(ENTRY*));

   assert(cache->slots[slot] != NULL);
   cache->slots[slot][idx] = entry;
   (cache->nentries[slot])++;
}

/** make cache */
CACHE* make_cache(
   void
   )
{
   CACHE* cache = (CACHE*) malloc(sizeof(CACHE));
   
   cache->slots = (ENTRY***) malloc(NSLOTS * sizeof(ENTRY**));
   cache->nentries = (int*) calloc(NSLOTS, sizeof(int));

   return cache;
}

/** free cache */
void free_cache(
   CACHE*                cache               /**< cache */
   )
{
   int slot;
   int i;
   
   for( slot = 0; slot < NSLOTS; slot++ )
   {
      for(i = 0; i < cache->nentries[slot]; i++ )
      {
         free_entry(cache->slots[slot][i]);
      }
      if( cache->nentries[slot] > 0 )
         free(cache->slots[slot]);
   }
   free(cache->slots);
   free(cache->nentries);
   free(cache);
}
   
