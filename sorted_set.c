#include "sorted_set.h"
#include <assert.h>
#include <string.h>
#include <stdlib.h>

#define MIN(a,b) (((a)<(b))?(a):(b))

/** 
 * sorted set, policy tree style
*/
struct sorted_set
{
  int*                   elements;           /**< the sorted set */
  int                    n;                  /**< size of subset */
  int*                   key;                /**< sorting key: if i < j then key[elements[i]] < key[elements[j]]  */
};

/** merge a pair of ordered subsequences to get a single ordered subsequence
 * left run is indices[ileft:iright-1]
 * right run is indices[iright:iend-1]
 */
static
void bottomupmerge(
   const int*            indices,            /**< global array */
   int                   ileft,              /**< index of first element of left subsequence */
   int                   iright,             /**< index of first element of right subsequence */
   int                   iend,               /**< iend-1 is index of last element of right subsequence */
   int*                  tmp,                /**< on output tmp[ileft:iend-1] has the merged subsequence */
   const double*         data_x_p            /**< data_x_p[i] is ordering key for i */
  )
{
  int i = ileft;
  int j = iright;
  int k;

  for( k = ileft; k < iend; k++ )
  {
    if( i < iright && (j >= iend || data_x_p[indices[i]] <= data_x_p[indices[j]]) )
      tmp[k] = indices[i++];
    else
      tmp[k] = indices[j++];
  }
}


/** sort an array by bottom up merge sort */
static
void bottomupmergesort(
   int*                  indices,            /**< array to sort */
   int*                  tmp,                /**< on output, sorted array */              
   int                   n,                  /**< length of array */
   const double*         data_x_p            /**< data_x_p[i] is ordering key for i */
  )
{
  int width;
  int i;
  const size_t size = n*sizeof(int);
  
  for( width = 1; width < n; width *= 2)
  {
    for( i = 0; i < n; i += 2*width)
      bottomupmerge(indices, i, MIN(i+width,n), MIN(i+2*width,n), tmp, data_x_p);

    memcpy(indices,tmp,size);
  }
}


/** insert an element into a simple sorted set, assuming it is not already there */
static
void insert_element(
   SIMPLE_SORTED_SET*    sorted_set,         /**< sorted set */
   int                   elt                 /**< element to insert */ 
  )
{
  unsigned int left;
  unsigned int right;
  int mid;
  int eltval;
  const int* key = sorted_set->key;
  int n;
  int* set;  
  
  assert(sorted_set != NULL);
  assert(elt >= 0);
  
  n = sorted_set->n;
  set = sorted_set->elements;  
  eltval = key[elt];

  assert(key != NULL);
  assert(set != NULL);
  assert(n >= 0);
  
  /* if set not empty look for an index 'mid' such that
   * key[set[mid]] < eltval < key[set[mid+1]]
   * can then insert elt by shifting set[mid+1], set[mid+2],,, and then
   * doing set[mid+1] = elt
   */
  
  if( n == 0)
  {
     set[0] = elt;
     sorted_set->n = 1;
     return;
  }
  else if( eltval < key[set[0]] )
  {
    /* insert at beginning */
    mid = -1;
  }
  else if ( eltval > key[set[n-1]] )
  {
    /* insert at end */
    mid = n-1;
  }
  else
  {
     /* since set is ordered there will be a value mid such that
      * key[set[mid]] < eltval < key[set[mid+1]]
      * where -1 < mid < n-1 (i.e. mid+1 is the index of an existing element )
      */

     assert(n > 1);
     left = 0;                /* smallest possible value for mid */
     right = n - 2;           /* largest possible value for mid */

     while( left <= right )
     {
        mid = (left + right)/2;  /* left <= mid <= right */
        
        if( key[set[mid]] > eltval )
        {
           /* key[set[mid]] too big, update biggest possible value for mid */
           right = mid - 1;
        }
        else if ( key[set[mid+1]] < eltval )
        {
           /* key[set[mid+1]] too small, update smallest possible value for mid*/
           left = mid + 1;
        }
        else
           break;
     }
  }
  
  /* move, if necessary */
  if( mid < n )
    memmove(set+mid+2,set+mid+1,(n-mid-1)*sizeof(int));
  
  /* insert at mid+1 */
  set[mid+1] = elt;
  sorted_set->n++;
}

/** remove an element from a simple sorted set, assuming it is already there 
 * @return 1 if the element is a member of the set, else -1
 */
static
int remove_element(
   SIMPLE_SORTED_SET*    sorted_set,         /**< sorted set */
   int                   elt                 /**< element to remove */
  )
{

   int left;
   int right;
   int mid;
   int midsave;
   int midval;
   int eltval;
   const int* rank = sorted_set->rank;
   int n;
   int* set;  
   int retval = -1;
   
   assert(sorted_set != NULL);
   
   n = sorted_set->n;
   set = sorted_set->elements;  
   eltval = rank[elt];

   left = 0;
   right = n-1;

   /* binary search for element */
   while( left <= right )
   {
      mid = (left+right)/2;
      midval = rank[set[mid]];
      
      if( midval < eltval )
         left = mid+1;
      else if ( midval > eltval )
         right = mid - 1;
      else
      {
         retval = 1;
         break;
      }
   }

   if( retval == 1 )
   {
      /* unless the last element, have to shift elements */
      if( mid < n-1)
         memmove(set+mid,set+mid+1,(n-mid-1)*sizeof(int));
      sorted_set->n--;
   }
   
   return retval;
}


/** insert elements into a sorted set, assuming they are not already there */
static
void insert_elements(
   SIMPLE_SORTED_SET*    sorted_set,         /**< sorted set */
   int                   nelts,              /**< number of elements to insert */ 
   int*                  elts                /**< elements to insert */ 
  )
{
   int i;

   /* just insert one at a time */
   for(i = 0; i < nelts; i++)
      insert_element(sorted_set, elts[i]);
}

/** add elements to the end of a sorted set, assuming they are not already there */
static
void add_elements_at_end(
   SIMPLE_SORTED_SET*    sorted_set,         /**< sorted set */
   int                   nelts,              /**< number of elements to insert */ 
   int*                  elts                /**< elements to insert */ 
  )
{
   int i;

   for(i = 0; i < nelts; i++)
      sorted_set->elements[sorted_set->n+i] = elts[i];

   sorted_set->n += nelts;
}


/** remove elements from a sorted set, assuming they are in the set */
static
void remove_elements(
   SIMPLE_SORTED_SET*    sorted_set,         /**< sorted set */
   int                   nelts,              /**< number of elements to remove */ 
   int*                  elts                /**< elements to insert */ 
  )
{
   int i;
   
   /* just remove one at a time */
   for(i = 0; i < nelts; i++)
      remove_element(sorted_set, elts[i]);
}

/** remove elements at the start of a sorted set */
static
void remove_elements_at_start(
   SIMPLE_SORTED_SET*    sorted_set,         /**< sorted set */
   int                   nelts               /**< number of elements to remove from start */ 
  )
{
   memmove(sorted_set->elements,sorted_set->elements+nelts,(sorted_set->n) - nelts);
   sorted_set->n -= nelts;
}


/** Determine whether a set is 'pure'.
 * A pure set is one where each unit has the same best action
 * @return 1 if the set is pure, else 0
 */
int is_pure(
   const SORTED_SET**    sorted_sets,        /**< sorted sets, representing a common unsorted set */
   const int*            best_actions        /**< best_actions[i] is the best action for unit i */
  )
{
  int i;
  int best_action;
  SORTED_SET* sorted_set;
  
  assert(sorted_sets != NULL);
  assert(sorted_sets[0] != NULL);
  assert(best_actions != NULL);
  
  /* arbitrarily use the first sorted set */
  sorted_set = sorted_sets[0];
  
  if( sorted_set->n > 0)
     best_action = best_actions[sorted_set->elements[0]];
  
  for( i = 1; i < sorted_set->n; i++)
     if( best_action != best_actions[sorted_set->elements[i]] )
        return 0;

  return 1;
}

/** get common size of sorted sets */
int get_size(
   SORTED_SET**          sorted_sets         /**< sorted sets */
   )
{
   assert(sorted_sets != NULL);
   assert(sorted_sets[0] != NULL);
   
   return (sorted_sets[0])->n;
}

/* find best action and its associated reward for a set of units */
void find_best_reward(
   const SORTED_SETS**   sorted_sets,        /**< sorted sets, representing a common unsorted set */
   const double*         data_y,             /**< gammas, data_y+(d*num_rows) points to values for reward d */
   int                   num_rows,           /**< number of units in full dataset */
   int                   num_cols_y,         /**< number of rewards/actions */
   WORKSPACE*            workspace,          /**< workspace */
   double*               best_reward,        /**< (pointer to) best reward */
   int*                  best_action,        /**< (pointer to) best action */
   )
{
   double* rewards;
   int d;
   SORTED_SET sorted_set;
   int i;
  
   assert(sorted_sets != NULL);
   assert(sorted_sets[0] != NULL);
   assert(data_y != NULL);
   assert(workspace != NULL);
   assert(best_reward != NULL);
   assert(best_action != NULL);
   
   rewards = get_double_array(workspace);
   
   /* get reward for each action */
   sorted_set = sorted_sets[0];
   for( d = 0; d < num_cols_y; d++ )
   {
      double* dyelt = data_y + d*num_rows;
      
      rewards[d] = 0.0;
      for( i = 0; i < sorted_set->n; i++)
         rewards[d] += dyelt[sorted_set->elements[i]];
   }
   
   *best_reward = rewards[0];
   *best_action = 0;
   for( d = 1; d < num_cols_y; d++ )
   {
      if( rewards[d] > *best_reward )
      {
         *best_reward = rewards[d];
         *best_action = d;
      }
   }
}

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
{
   SORTED_SET* left_sorted_setp;
   SORTED_SET* right_sorted_setp;
   int nmoved;
   int pp;

   assert(left_sorted_sets != NULL);
   assert(right_sorted_sets != NULL);
   assert(p > 0 && p < num_cols_x);
   assert(data_xp != NULL);
   assert(num_cols_x > 0);
   assert(workspace != NULL);
   assert(splitval != NULL);
   assert(elts != NULL);
   assert(nelts != NULL);
   
   left_sorted_setp = left_sorted_sets[p];
   right_sorted_setp = right_sorted_sets[p];

   /* nothing to move from right to left */
   if( right_sorted_setp->n == 0 )
      return 0;

   /* splitting value is just first covariate value on the right */
   *splitval = data_xp[right_sorted_setp->elements[0]]; 

   /* if left is non-empty the last element on left must have a strictly lower pth covariate value
      that first element on right */
   assert(left_sorted_setp->n == 0 || data_xp[left_sorted_setp->elements[(left_sorted_setp->n)-1]] < *splitval); 

   /* find any additional units on right with *splitval as covariate value */
   nmoved = 1;
   while( nmoved < right_sorted_setp->n && data_xp[right_sorted_setp->elements[nmoved]] == *splitval )
      nmoved++;

   /* update left and right sorted sets of units for covariates other than splitting covariate */
   for( pp = 0; pp < num_cols_x; pp++)
   {
      if( pp != p )
      {
         insert_elements(left_sorted_sets[pp], nmoved, right_sorted_setp);
         remove_elements(right_sorted_sets[pp], nmoved, right_sorted_setp);
      }
   }

   /* update left and right sorted sets for splitting covariate */
   add_elements_at_end(left_sorted_setp, nmoved, right_sorted_setp);
   remove_elements_at_start(right_sorted_setp, nmoved);

   /* return set of moved units */
   *nelts = nmoved;
   *elts = left_sorted_setp->elements + left_sorted_setp->n - nmoved;

   return 1;
   
}

/* make a 'shallow' copy of source sorted sets */
SORTED_SET** shallow_copy_sorted_sets(
   SORTED_SET**          sources,            /**< source sorted sets */
   int                   nsets              /**< number of sources */
   )
{
   SORTED_SET** targets;

   int p;

   targets = (SORTED_SET**) malloc(nsets*sizeof(SORTED_SET*));

   for( p = 0; p < nsets; p++ )
   {
      SORTED_SET* source = sources[p];
      SORTED_SET* target;
      
      target = (SORTED_SET*) malloc(sizeof(SORTED_SET));
      target->elements = (int*) malloc(source->n*sizeof(int));
      target->n = source->n;
      /* just copy pointer */
      target->key = source->key;

      memcpy(target->elements, source->elements, (source->n)*sizeof(int));

      targets[p] = target;
   }

   return targets;
}


/** make a sorted set from elements 0,...,num_indices-1 
 * this function allocates all the necessary memory for the sorted set
 * @return The sorted set
 */
static
SORTED_SET* make_sorted_set(
   int                   num_indices,        /**< size of set */
   const double*         data_xx,            /**< data_xx[i] is ordering key for element i */
   int*                  tmp                 /**< temporary storage for ordering of size at least num_indices */
  )
{
  SORTED_SET* sorted_set;
  int i;
  int* elements;
  int* key;

  sorted_set = (SORTED_SET*) malloc(sizeof(SORTED_SET));
  elements = (int* ) malloc(num_indices*sizeof(int));
  key = (int* ) malloc(num_indices*sizeof(int));
  for( i = 0; i < num_indices; i++)
    elements[i] = i;

  /* sort the elements using data_xx values as key */
  bottomupmergesort(elements, tmp, num_indices, data_xx);

  for( i = 0; i < num_indices; i++)
    key[elements[i]] = i;
  
  sorted_set->elements = elements;
  sorted_set->n = num_indices;
  sorted_set->key = key;
  
  return sorted_set;
}

/** initialise so that each left_sorted_set (for each covariate) is empty and each 
 *  right_sorted_set is a copy of the sorted set (for that covariate) 
*/
void initialise_sorted_sets(
   SORTED_SET**          sorted_sets,        /**< input sorted sets */
   int                   depth,              /**< depth of associated node */
   int                   num_cols_x,         /**< number of covariates */
   WORKSPACE*            workspace,          /**< workspace */
   SORTED_SET***         left_sorted_sets,   /**< pointer to output left sets */
   SORTED_SET***         right_sorted_sets   /**< pointer to output right sets */
   )
{
   int p;
   SORTED_SET** lefts;
   SORTED_SET** rights;
   /* get common size of sorted sets */
   int n = (sorted_sets[0])->n;
   
   lefts = get_left_sorted_sets(workspace,depth);
   rights = get_right_sorted_sets(workspace,depth);

   /* do not need to set key, since this is fixed */
   for( p = 0; p < num_cols_x; p++ )
   {
      (lefts[p])->n = 0;
      memcpy((rights[p])->elements, (sorted_sets[p])->elements, n*sizeof(int));
      (rights[p])->n = n;
   }

   *left_sorted_sets = lefts;
   *right_sorted_sets = rights;
}

SORTED_SET** make_initial_sorted_sets(
   const double*         data_x,             /**< covariates, data_x+(j*num_rows) points to values for covariate j */
   int                   num_rows,           /**< number of units in full dataset */
   int                   num_cols_x          /**< number of covariates */
   )
{

   SORTED_SET** initial_sorted_set;
   int* tmp_indices;
   int p;
   
   initial_sorted_sets = (SORTED_SET**) malloc(num_cols_x*sizeof(SORTED_SET*));

   tmp_indices = (int*) malloc(num_rows*sizeof(int));

   /* make initial sorted sets */
   for( p = 0; p < num_cols_x; p++)
      initial_sorted_sets[p] = make_sorted_set(num_rows,data_x+p*num_rows,tmp_indices);
   
   free(tmp_indices);

   return initial_sorted_sets;
}

static
void free_sorted_set(
   SORTED_SET*           sorted_set          /**< sorted set */
   )
{
   free(sorted_set->elements);
   free(sorted_set->key);
   free(sorted_set);
}

static
void shallow_free_sorted_set(
   SORTED_SET*           sorted_set          /**< sorted set */
   )
{
   free(sorted_set->elements);
   free(sorted_set);
}


void free_sorted_sets(
   SORTED_SET**          sorted_sets,        /**< sorted sets */
   int                   nsets               /**< number of sorted sets */
   )
{
   int p;

   for(p = 0; p < nsets; p++ )
      free_sorted_set(sorted_sets[p]);
   free(sorted_sets);
}

void shallow_free_sorted_sets(
   SORTED_SET**          sorted_sets,        /**< sorted sets */
   int                   nsets               /**< number of sorted sets */
   )
{
   int p;

   for(p = 0; p < nsets; p++ )
      shallow_free_sorted_set(sorted_sets[p]);
   free(sorted_sets);
}

/** make a very shallow copy of a sorted set */
void very_shallow_copy(
   SORTED_SET*           source,             /**< source sorted set */
   SORTED_SET*           target              /**< target sorted set */
   )
{
   target->elements = sorted_set->elements;
   target->key = sorted_set->key;
   target->n = sorted_set->n;
}

/** remove elements from the start of sorted set with the same covariate value */
int next_shallow_split(
   SORTED_SET*           right_sorted_set,   /**< sorted set */ 
   const double*         data_xp,            /**< values for covariate to split on */
   double*               splitval,           /**< (pointer to) found value to split on */
   int**                 elts,               /**< (pointer to) the elements moved */
   int*                  nelts               /**< (pointer to) number of elements moved */
   )
{
   /* nothing to move from right to left */
   if( right_sorted_set->n == 0 )
      return 0;

   /* splitting value is just first covariate value on the right */
   *splitval = data_xp[right_sorted_set->elements[0]]; 

   /* find any additional units on right with *splitval as covariate value */
   nmoved = 1;
   while( nmoved < right_sorted_set->n && data_xp[right_sorted_set->elements[nmoved]] == *splitval )
      nmoved++;

   /* record what's removed */
   *nelts = nmoved;
   *elts = right_sorted_set->elements;

   /* update set */
   /* NB elements pointers is moved, so original value must be available 
      somewhere to allow eventual freeing */
   right_sorted_set->elements += nmoved;
   right_sorted_set->n -= nmoved;
   
   return 1;
}

void find_nosplit_rewards(
   const SORTED_SET**    sorted_sets,        /**< sorted sets */
   const double*         data_y,             /**< gammas, data_y+(d*num_rows) points to values for reward d */
   int                   num_rows,           /**< number of rows in the data */
   double*               nosplit_rewards     /**< space for computed no split rewards */
   )
{
   int d;
   int i;
   double* dyelt;
   int elt;
   SORTED_SET* sorted_set = sorted_sets[0];
   
   /* find reward for each action if no split were done */
   for( d = 0; d < num_cols_y; d++ )
      nosplit_rewards[d] = 0.0;
   for( i = 0; i < sorted_set->n; i++ )
   {
      elt = sorted_set->elements[i];
      dyelt = data_y + elt;
      for( d = 0; d < num_cols_y; d++ )
      {
         nosplit_rewards[d] += *dyelt;
         dyelt += num_rows;
      }
   }
}

SORTED_SET* make_uninitialised_sorted_set(
   void
   )
{
   return (SORTED_SET*) malloc(sizeof(SORTED_SET));
}
