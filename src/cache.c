#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "cache.h"
#include "tree.h"
#include "assert.h"
#include <limits.h>
#include <stdint.h>

#define NSLOTS 100
#define BLOCKSIZE 5

/* following macros adapted from https://c-faq.com/misc/bitsets.html */

#define NBITS (CHAR_BIT * sizeof(uint32_t))       /**< number of bits in a uint32_t int - will be 32 */
#define BITMASK(b) (1 << ((b) % NBITS))           /**< given b, sets only bth bit in a slot (all others zero) */
#define BITSLOT(b) ((b) / NBITS)                  /**< finds correct slot for an integer */
#define BITSET(a, b) ((a)[BITSLOT(b)] |= BITMASK(b))     /**< given integer b and bitset a, sets correct bit for b */
#define BITCLEAR(a, b) ((a)[BITSLOT(b)] &= ~BITMASK(b))  /**< given integer b and bitset a, clears correct bit for b */
#define BITTEST(a, b) ((a)[BITSLOT(b)] & BITMASK(b))     /**< given integer b and bitset a, test whether correct bit for b is set */
#define BITNSLOTS(nb) ((nb + NBITS - 1) / NBITS)         /**< compute number of uint32_t ints needed for 'universe' of size nb */


struct entry
{
   uint32_t*             key;                /**< bitset representation of set */
   int                   nelts;              /**< number of elements in the set */
   int                   depth;              /**< depth of tree */
   NODE*                 tree;               /**< optimal tree for set */
};
typedef struct entry ENTRY;

struct cache
{
   ENTRY***              slots;              /** slots */
   int*                  nentries;           /** number of entries in each slot */
   int*                  sizes;              /** available space for each slot */
   int                   nints;              /** number of uint32_t ints required for a bitset representation */
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

/** is this the entry for given set and depth? */
static
int match(
   const uint32_t*       key,                /**< bitset key for set */
   int                   nints,              /**< number of ints for a bitset */
   int                   nelts,              /**< number of elts in the set */
   int                   depth,              /**< depth of tree */
   const ENTRY*          entry               /**< entry */
   )
{
   int i;
   
   if( depth != entry->depth || nelts != entry->nelts )
      return 0;

   for( i = 0; i < nints; i++ )
      if( key[i] != entry->key[i] )
         return 0;

   return 1;
}

/** create an entry */
static
ENTRY* make_entry(
   const CACHE*          cache,              /**< cache */
   int                   nelts,              /**< number of elts in set */
   const ELEMENT*        elts,               /**< the set */
   int                   depth,              /**< depth of tree */
   const NODE*           tree                /**< optimal tree for set */
   )
{
   ENTRY* entry = (ENTRY*) malloc(sizeof(ENTRY));
   NODE* tree_cp = make_tree(depth);
   uint32_t* key = (uint32_t*) calloc(cache->nints, sizeof(uint32_t));
   int i;
   
   for( i = 0; i < nelts; i++)
      BITSET(key,elts[i]);
   
   tree_copy(tree, tree_cp);

   entry->key = key;
   entry->nelts = nelts;
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
   free(entry->key);
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
   uint32_t* key = (uint32_t*) calloc(cache->nints, sizeof(uint32_t));
   int result = 0;
   
   for( i = 0; i < nelts; i++)
      BITSET(key,elts[i]);
   
   for( i = 0; i < cache->nentries[slot]; i++ )
   {
      if( match((const uint32_t*) key, cache->nints, nelts, depth, (const ENTRY*) cache->slots[slot][i]) )
      {
         tree_copy(cache->slots[slot][i]->tree,tree);
         result = 1;
         break;
      }
   }

   free(key);
   return result;
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
   ENTRY* entry = make_entry(cache, nelts, elts, depth, tree);

   assert(cache != NULL);
   assert(nelts >= 0);
   assert(elts != NULL);
   assert(depth >= 0);
   assert(tree != NULL);

   assert(slot >= 0);
   assert(slot < NSLOTS);
   assert(entry != NULL);

   if( cache->nentries[slot] == cache->sizes[slot] )
   {
      /* need more space */
      cache->sizes[slot] += BLOCKSIZE;
      cache->slots[slot] = (ENTRY**) realloc(cache->slots[slot], (cache->sizes[slot])*sizeof(ENTRY*));
   }
   cache->slots[slot][(cache->nentries[slot])++] = entry;
}

/** make cache */
CACHE* make_cache(
   int                   num_rows            /**< number of units */
   )
{
   int slot;
   CACHE* cache = (CACHE*) malloc(sizeof(CACHE));
   
   cache->slots = (ENTRY***) malloc(NSLOTS * sizeof(ENTRY**));
   cache->nentries = (int*) calloc(NSLOTS, sizeof(int));
   cache->sizes = (int*) malloc(NSLOTS * sizeof(int));
   cache->nints = BITNSLOTS(num_rows);

   /* assign initial space for BLOCKSIZE entries in each slot */
   for( slot = 0; slot < NSLOTS; slot++ )
   {
      cache->slots[slot] = (ENTRY**) malloc(BLOCKSIZE * sizeof(ENTRY*));
      cache->sizes[slot] = BLOCKSIZE;
   }
   
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
      free(cache->slots[slot]);
   }
   free(cache->slots);
   free(cache->nentries);
   free(cache->sizes);
   free(cache);
}
   
