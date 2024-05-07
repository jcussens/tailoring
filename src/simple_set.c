/** @file simple_set.c
 *  @brief Implements a set of units as a single array
 *  @author James Cussens
 */
#include "units.h"
#include "workspace.h"
#include "simple_set.h"
#include "strategy.h"
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#ifdef PRINTING_ALLOWED
#include <stdio.h> /* only for debugging */
#endif
#define MIN(a,b) (((a)<(b))?(a):(b))         /**< compute minimum of 2 values */
#define MAX(a,b) (((a)>(b))?(a):(b))         /**< compute maximum of 2 values */

#define MAXKEYVALS 30                        /**< maximum number of key values to use counting sort (otherwise use radix sort) */

typedef int KEY;

/** 
 * simple set
*/
struct simple_set
{
   ELEMENT*              elements;           /**< pointer to space for the set */
   int                   size;               /**< size of space for set */
   int                   start;              /**< elements[start] is the first element of the set */ 
   int                   n;                  /**< number of elements in set */
   KEY**                 keys;               /**< let data_xp be covariate values for covariate p then 
                                              * (1) if data_xp[i] < data_xp[j] then keys[p][i] < keys[p][j] and
                                              * (2) if data_xp[i] == data_xp[j] then keys[p][i] == keys[p][j] */
   int*                  nkeyvals;           /**< nkeyvals[p] is the number of distinct data_xp values */
};

/** counting sort */
static
void counting_sort(
   const ELEMENT*        a,                  /**< input array to sort */
   int                   n,                  /**< length of input array */
   ELEMENT*              b,                  /**< will contain sorted array */
   int                   k,                  /**< number of distinct key values */
   const KEY*            key,                /**< key[elem] is the key value for elem */
   ELEMENT*              tmp                 /**< tmp array, at least as long as b */
   )
{
   int i;
   int count[MAXKEYVALS] = {0};

   assert( a != NULL );
   assert( n >= 0 );
   assert( b != NULL );
   assert( k >= 0 && k <= MAXKEYVALS );
   assert( key != NULL );
   
   /* get count for each key value */
   for(i = 0; i < n; i++ )
   {
      tmp[i] = key[a[i]];
      count[tmp[i]]++;
   }

   /* alter so that c[i] is number of key vals <= i */ 
   for(i = 1; i < k; i++ )
      count[i] += count[i-1];

   /* put each element of a into correct place in b */
   for(i = n-1; i >=0; i--)
      b[(count[tmp[i]]--)-1] = a[i];
}

/** counting sort on a particular digit */
static
void counting_sort_radix(
   const ELEMENT*        a,                  /**< input array to sort */
   int                   n,                  /**< length of input array */
   ELEMENT*              b,                  /**< will contain sorted array */
   int*                  tmp2,               /**< temporary working space */
   const KEY*            key,                /**< key[elem] is the key value for elem */
   int                   exp                 /**< = 1, 10, 100, .. indicates which digit to sort on */
   )
{
   int i;
   int count[10] = {0};

   assert( a != NULL );
   assert( n >= 0 );
   assert( b != NULL );
   assert( key != NULL );
   
   /* get count for each key value */
   for(i = 0; i < n; i++ )
   {
      tmp2[i] = (key[a[i]]/exp)%10;
      count[tmp2[i]]++;
   }

   /* alter so that c[i] is number of key vals <= i */ 
   for(i = 1; i < 10; i++ )
      count[i] += count[i-1];

   /* put each element of a into correct place in b */
   for(i = n-1; i >=0; i--)
      b[(count[tmp2[i]]--)-1] = a[i];
}

/** radix sort */
static
void radix_sort(
   const ELEMENT*        a,                  /**< input array to sort */
   int                   n,                  /**< length of input array */
   ELEMENT*              b,                  /**< will contain sorted array */
   ELEMENT*              tmp,                /**< tmp array, at least as long as b */
   int*                  tmp2,               /**< tmp array, at least as long as b */
   const KEY*            key,                /**< key[elem] is the key value for elem */
   int                   maxkey              /**< maximal value of key[elem] */
   )
{
   int exp;
   ELEMENT* tmp3;

   /* sort a on least significant digit and put result in b */
   counting_sort_radix(a, n, b, tmp2, key, 1);
   
   for( exp = 10; maxkey / exp  > 0; exp *= 10 )
   {
      counting_sort_radix((const ELEMENT*) b, n, tmp, tmp2, key, exp);
      tmp3 = tmp;
      tmp = b;
      b = tmp3;
   }
}

#if 0
static
void swap(
   ELEMENT*              b,                  /**< array to be sorted */
   int                   i,
   int                   j
   )
{
   int tmp;
   
   tmp = b[i];
   b[i] = b[j];
   b[j] = tmp;
}

static
void quicksort(
   ELEMENT*              b,                  /**< array to be sorted */
   int                   left,               /**< length of array to be sorted */
   int                   right,              /**< length of array to be sorted */
   const KEY*            key                 /**< key[elem] is key value for elem */
   )
{
   int i;
   int last;

   if( left >= right )
      return;
   swap(b, left, (left+right)/2);
   last = left;
   for( i = left+1; i <= right; i++)
      if( key[b[i]] < key[b[left]] )
         swap(b, ++last, i);
   swap(b, left, last);
   quicksort(b, left, last-1, key);
   quicksort(b, last+1, right, key);
}
#endif


/** sort an array of integers according to key values */
static
void sort_units(
   const ELEMENT*        a,                  /**< array to be sorted */
   int                   n,                  /**< length of array to be sorted */
   const KEY*            key,                /**< key[elem] is key value for elem */
   int                   nkeyvals,           /**< key values are 0,...,nkeyvals-1 */
   ELEMENT*              tmp,                /**< temporary array of length at least n */
   int*                  tmp2,               /**< temporary array of length at least n */
   ELEMENT*              b                   /**< output array of length at least n */
   )
{
   if( nkeyvals <= MAXKEYVALS )
      counting_sort(a, n, b, nkeyvals, key, tmp);
   else
   {
      /* int i; */
      /* for(i = 0; i < n; i++) */
      /*    b[i] = a[i]; */
      /* quicksort(b, 0, n-1, key); */
      radix_sort(a, n, b, tmp, tmp2, key, nkeyvals-1); 
   }
}



/** merge a pair of ordered subsequences to get a single ordered subsequence
 * left run is indices[ileft:iright-1]
 * right run is indices[iright:iend-1]
 */
static
void bottomupmerge(
   const ELEMENT*        indices,            /**< global array */
   int                   ileft,              /**< index of first element of left subsequence */
   int                   iright,             /**< index of first element of right subsequence */
   int                   iend,               /**< iend-1 is index of last element of right subsequence */
   ELEMENT*              tmp,                /**< on output tmp[ileft:iend-1] has the merged subsequence */
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
   ELEMENT*              indices,            /**< array to sort */
   ELEMENT*              tmp,                /**< on output, sorted array */              
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

/** return array of elements and size of array */
void simple_set_elements(
   const SIMPLE_SET*     simple_set,         /**< set */
   ELEMENT**             elts,               /**< (pointer to) elements */
   int*                  nelts               /**< (pointer to) number of elements */
   )
{
   *elts = simple_set->elements+simple_set->start;
   *nelts = simple_set->n;
}

/** (for debugging) check that a simple_set is valid. If p!=-1 then check it is ready for splitting on covariate p */
int simple_set_units_ok(
   const SIMPLE_SET*     simple_set,         /**< set */
   int                   p,                  /**< if !=-1 then units must be ready for splitting on covariate p */
   const double*         data_x,             /**< covariates, data_x+(j*num_rows) points to values for covariate j */
   int                   num_rows,           /**< number of units in full dataset */
   int                   num_cols_x          /**< number of covariates */
   )
{

   assert( simple_set != NULL );
   assert( num_rows >= 1 );
   assert( num_cols_x >= 0 );

   if( simple_set->start < 0 )
      return 0;

   if( simple_set->n < 0 )
      return 0;

   if( simple_set->size < 0 )
      return 0;

   if( simple_set->size < simple_set->n )
      return 0;

   if( !(simple_set->start < simple_set->size) )
      return 0;

   if( num_cols_x > 0 )
   {
      int pp;
      
      if( simple_set->keys == NULL )
         return 0;
      if( simple_set->nkeyvals == NULL )
         return 0;

      for( pp = 0; pp < num_cols_x; pp++ )
      {
         int i;
         
         if(simple_set->keys[pp] == NULL )
            return 0;
         if(simple_set->nkeyvals[pp] < 1 )
            return 0;

         for( i = 0; i < num_rows; i++ )
         {
            if( !(simple_set->keys[pp][i] < simple_set->nkeyvals[pp]) )
               return 0;
         }
      }
   }

   /* now check ordered for covariate p */
   if( p != -1 )
   {
      int idx;
      const double* data_xp = data_x + p*num_rows;
      
      for( idx = simple_set->start + 1; idx < simple_set->start + simple_set->n; idx++)
         if( data_xp[simple_set->elements[idx-1]] > data_xp[simple_set->elements[idx]] )
            return 0;
   }
   
   return 1;
}

/** get the reward for a set if all units in the set were assigned their best action *
 * @return the reward for a set if all units in the set were assigned their best action *
 */
double simple_set_get_reward_ub(
   const SIMPLE_SET*     simple_set,         /**< set */
   const double*         data_y,             /**< gammas, data_y+(d*num_rows) points to values for reward d */
   int                   num_rows,           /**< number of units in full dataset */
   const int*            best_actions        /**< best_actions[i] is the best action for unit i */
   )
{
   int i;
   double ub = 0.0;

   for(i = simple_set->start; i < simple_set->start + simple_set->n; i++)
   {
      ub += *(data_y + num_rows*best_actions[simple_set->elements[i]] + simple_set->elements[i]);
      /* printf("%d %d %d %g\n", i, simple_set->elements[i], best_actions[simple_set->elements[i]], ub); */
   }     
   return ub;
}

/** Determine whether a set is 'pure'.
 * A pure set is one where each unit has the same best action
 * @return 1 if the set is pure, else 0
 */
int simple_set_is_pure(
   const SIMPLE_SET*     simple_set,         /**< set */
   const int*            best_actions        /**< best_actions[i] is the best action for unit i */
  )
{
  int i;
  int best_action;
  
  assert( simple_set != NULL );
  assert( best_actions != NULL );
  
  if( simple_set->n > 0)
     best_action = best_actions[simple_set->elements[simple_set->start]];
  
  for( i = simple_set->start+1; i < simple_set->start+simple_set->n; i++)
     if( best_action != best_actions[simple_set->elements[i]] )
        return 0;

  return 1;
}

/** get size of set */
int simple_set_get_size(
   const SIMPLE_SET*     simple_set          /**< set */
   )
{
   assert(simple_set != NULL);
   
   return simple_set->n;
}

/** find best action and its associated reward for a set of units */
void simple_set_find_best_action(
   const SIMPLE_SET*     simple_set,         /**< simple set */
   const double*         data_y,             /**< gammas, data_y+(d*num_rows) points to values for reward d */
   int                   num_rows,           /**< number of units in full dataset */
   int                   num_cols_y,         /**< number of rewards/actions */
   WORKSPACE*            workspace,          /**< workspace */
   double*               best_reward,        /**< (pointer to) best reward */
   int*                  best_action         /**< (pointer to) best action */
   )
{
   double* rewards;
   int d;
   int i;
  
   assert( simple_set != NULL );
   assert( data_y != NULL );
   assert( workspace != NULL );
   assert( best_reward != NULL );
   assert( best_action != NULL );
   
   rewards = get_rewards_space(workspace);
   
   /* get reward for each action */
   for( d = 0; d < num_cols_y; d++ )
   {
      const double* dyelt = data_y + d*num_rows;
      
      rewards[d] = 0.0;
      for( i = simple_set->start; i < simple_set->start + simple_set->n; i++)
         rewards[d] += dyelt[simple_set->elements[i]];
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

/** given that left_set and right_set are sorted according to covariate p, find next split (if any)  and associated split value
 * if there is a split, then both left_set and right_set are updated accordingly.
 * If splitting on a binary variable then a different approach is taken
 * @return 1 if there is a next split, otherwise 0
 */
int simple_set_next_split(
   const STRATEGY*       strategy,           /**< tree-building strategy */
   const SIMPLE_SET*     simple_set,         /**< input set */
   SIMPLE_SET*           left_set,           /**< left set sorted on covariate p */
   SIMPLE_SET*           right_set,          /**< right set sorted on covariate p */ 
   int                   p,                  /**< covariate to split on */
   const double*         data_xp,            /**< values for covariate to split on */
   int                   num_cols_x,         /**< number of covariates */
   double*               splitval,           /**< (pointer to) found value to split on */
   ELEMENT**             elts,               /**< (pointer to) elements moved */
   int*                  nelts,              /**< (pointer to) number of elements moved */
   int                   splitcount          /**< number of previous splits */
   )
{
   int leftend;
   
   assert(left_set != NULL);
   assert(right_set != NULL);
   assert(p >= 0 && p < num_cols_x);
   assert(data_xp != NULL);
   assert(num_cols_x > 0);
   assert(splitval != NULL);

   if( exploit_binaryvars(strategy) && (simple_set->nkeyvals[p] == 2) )
   {
      int i;
      const KEY* keys = simple_set->keys[p];
      
      assert( splitcount == 0 || splitcount == 1 );

      /* since can only be one split for a binary variable, just return if it's been done */
      if( splitcount == 1 )
         return 0;

      /* check left and right set initialised as expected */
      assert(left_set->n == 0);
      assert(left_set->start == 0);
      assert(left_set->elements != NULL);
      assert(right_set->n == 0);
      assert(right_set->start == 0);
      assert(right_set->elements != NULL);
      
      /* just put lower key value on left and higher on right */
      for( i = simple_set->start; i < simple_set->start + simple_set->n; i++)
      {
         ELEMENT element = simple_set->elements[i];
         
         if( keys[element] == 0 )
            left_set->elements[(left_set->n)++] = element;
         else
            right_set->elements[(right_set->n)++] = element;
      }

      *nelts = left_set->n;
      *elts = left_set->elements;

      if( *nelts > 0 )
      {
         *splitval = data_xp[left_set->elements[0]];
         return 1;
      }
      else
      {
         return 0;
      }
   }
   
   *nelts = 0;

   /* nothing to move from right to left */
   if( right_set->n == 0 )
      return 0;

   /* splitting value is just first covariate value on the right */
   *splitval = data_xp[right_set->elements[right_set->start]]; 

   /* leftend is the position just beyond any elements of left */
   leftend = left_set->start + left_set->n;
   
   /* if left is non-empty the last element on left must have a strictly lower pth covariate value
      that first element on right */
   assert(left_set->n == 0 || data_xp[left_set->elements[leftend-1]] < *splitval); 

   /* move all units on right with *splitval as covariate value to left */
   while( right_set->n > 0 && data_xp[right_set->elements[right_set->start]] == *splitval )
   {
      left_set->elements[leftend++] = right_set->elements[right_set->start++];
      left_set->n++;
      right_set->n--;
      (*nelts)++;
   }

   /* if all moved from right to left this is not a split */
   if( right_set->n == 0 )
      return 0;

   *elts = left_set->elements + leftend - *nelts;
   
   return 1;
   
}

/** make a 'shallow' copy of source set */
SIMPLE_SET* simple_set_shallow_copy_units(
   const SIMPLE_SET*     source,            /**< source  */
   int                   num_cols_x         /**< number of covariates */
   )
{
   SIMPLE_SET* target;

   target = (SIMPLE_SET*) malloc(sizeof(SIMPLE_SET));
   target->elements = (ELEMENT*) malloc(source->size * sizeof(ELEMENT));
   memcpy(target->elements, source->elements, (source->size)*sizeof(int));
   target->size = source->size;
   target->start = source->start;
   target->n = source->n;
   target->keys = source->keys;
   target->nkeyvals = source->nkeyvals;

   return target;
}

/** get integer key values after sorting units on a double key function, elements with the same double key value
 * get the same integer key value. also compute the number of distinct keys (== number of distinct double key values)
 * @return key where key[elem] is key value for elem
 */
static
KEY* get_key(
   ELEMENT*              elements,           /**< elements (unordered) */
   const double*         data_xp,            /**< covariate values for some covariate for sorting */
   int                   num_rows,           /**< number of units in full dataset */
   ELEMENT*              tmp,                /**< temporary working space */
   int*                  nkeyvals            /**< number of distinct data_xp values */
   )
{
   int i;
   KEY* keysp;
   int keyval;
   
   /* sort the elements using data_xp values as key */
   bottomupmergesort(elements, tmp, num_rows, data_xp);

   keysp = (KEY*) malloc(num_rows * sizeof(KEY));

   keyval = 0;
   keysp[elements[0]] = keyval;
   *nkeyvals = 1;

   for(i = 1; i < num_rows; i++ )
   {
      if( data_xp[elements[i]] == data_xp[elements[i-1]] )
      {
         keysp[elements[i]] = keyval;
      }
      else
      {
         keysp[elements[i]] = ++keyval;
         (*nkeyvals)++;
      }
   }

   return keysp;
}




/** initialise so that left_simple_set is empty and  
 *  right_simple_set is a copy of the simple_set sorted on covariate p
*/
void simple_set_initialise_units(
   const STRATEGY*       strategy,           /**< tree-building strategy */
   const SIMPLE_SET*     simple_set,         /**< input set */
   int                   p,                  /**< splitting covariate */
   int                   depth,              /**< depth of associated node */
   int                   num_cols_x,         /**< number of covariates */
   WORKSPACE*            workspace,          /**< workspace */
   SIMPLE_SET**          left_simple_set,    /**< pointer to output left set */
   SIMPLE_SET**          right_simple_set    /**< pointer to output right set */
   )
{
   SIMPLE_SET* left = NULL;
   SIMPLE_SET* right = NULL;

   /* TODO change names of 2 functions */
   left = (SIMPLE_SET*) get_left_sorted_sets(workspace,depth);
   right = (SIMPLE_SET*) get_right_sorted_sets(workspace,depth);

   assert( left != NULL );
   assert( right != NULL );

   /* make left empty and ensure any additions added at start */
   left->n = 0;
   left->start = 0;

   right->start = 0;

   if( exploit_binaryvars(strategy) && simple_set->nkeyvals[p] == 2 )
   {
      /* do basically nothing */
      right->n = 0;
   }
   else
   {
      right->n = simple_set->n;

      /* make right a copy of input, but with elements ordered according to covariate p */
      sort_units(simple_set->elements + simple_set->start, simple_set->n, simple_set->keys[p], simple_set->nkeyvals[p],
         ((SIMPLE_SET*) get_tmpunits(workspace))->elements, get_tmp2(workspace), right->elements);
   }
   *left_simple_set = left;
   *right_simple_set = right;
}


/** initialise so that right_simple_set is a copy of the input set and we are ready to look for depth=1 splits using covariate p */
void simple_set_shallow_initialise_units(
   const STRATEGY*       strategy,           /**< tree-building strategy */
   const SIMPLE_SET*     simple_set,         /**< input set */
   int                   p,                  /**< splitting covariate */
   int                   num_cols_x,         /**< number of covariates */
   WORKSPACE*            workspace,          /**< workspace */
   SIMPLE_SET**          right_simple_set    /**< pointer to output right set */
   )
{

   SIMPLE_SET* right;

   right = (SIMPLE_SET*) get_right_sorted_sets(workspace,1);

   assert( right != NULL );

   if( !(exploit_binaryvars(strategy) && simple_set->nkeyvals[p] == 2) )
   {
      /* make right a copy of input, but with elements ordered according to covariate p */
      right->n = simple_set->n;
      right->start = 0;
      sort_units(simple_set->elements + simple_set->start, simple_set->n, simple_set->keys[p], simple_set->nkeyvals[p],
         ((SIMPLE_SET*) get_tmpunits(workspace))->elements, get_tmp2(workspace), right->elements);
   }

   *right_simple_set = right;
}


#ifdef PRINTING_ALLOWED
/** print out a simple set (for debugging only) */
void print_simple_set(
   const SIMPLE_SET*     simple_set,         /**< units */
   const double*         data_x,             /**< covariates, data_x+(j*num_rows) points to values for covariate j */
   int                   num_rows,           /**< number of units in full dataset */
   int                   num_cols_x          /**< number of covariates */
   )
{
   int i;
   int p;

   for(i = simple_set->start; i < simple_set->start + simple_set->n; i++)
      printf("%d,",simple_set->elements[i]);
   printf("\n");
   
   for( p = 0; p < num_cols_x; p++)
   {
      const double* data_xp = data_x+p*num_rows;
      
      printf("%d:\n", p);
      for(i = simple_set->start; i < simple_set->start + simple_set->n; i++)
         printf("(%g,%d) ", data_xp[simple_set->elements[i]], (simple_set->keys[p])[simple_set->elements[i]]);
      printf("\n");
   }
}
#endif

/** make a set of units from covariate data 
 * @return the set of units
 */
SIMPLE_SET* simple_set_make_units(
   const double*         data_x,             /**< covariates, data_x+(j*num_rows) points to values for covariate j */
   int                   num_rows,           /**< number of units in full dataset */
   int                   num_cols_x          /**< number of covariates */
   )
{

   ELEMENT* tmp_indices;
   int p;
   SIMPLE_SET* initial_simple_set; 
   int i;
   
   tmp_indices = (ELEMENT*) malloc(num_rows*sizeof(ELEMENT));
   
   initial_simple_set = (SIMPLE_SET*) malloc(sizeof(SIMPLE_SET));
   initial_simple_set->elements = (ELEMENT*) malloc(num_rows * sizeof(ELEMENT));
   for( i = 0; i < num_rows; i++ )
      initial_simple_set->elements[i] = i;
   initial_simple_set->size = num_rows;
   initial_simple_set->start = 0;
   initial_simple_set->n = num_rows;
   initial_simple_set->keys = (int**) malloc(num_cols_x * sizeof(int*));
   initial_simple_set->nkeyvals = (int*) malloc(num_cols_x * sizeof(int));
   for( p = 0; p < num_cols_x; p++)
   {
      int nkeyvals;
      
      initial_simple_set->keys[p] = get_key(initial_simple_set->elements, data_x+p*num_rows, num_rows, tmp_indices, &nkeyvals);
      initial_simple_set->nkeyvals[p] = nkeyvals;
   }
   
   free(tmp_indices);

   return initial_simple_set;
}

/** free a set of units */
void simple_set_free_units(
   SIMPLE_SET*           simple_set,         /**< set */
   int                   num_cols_x          /**< number of covariates */
   )
{
   int p;

   for( p = 0; p < num_cols_x; p++)
   {
      free(simple_set->keys[p]);
   }
   free(simple_set->keys);
   free(simple_set->nkeyvals);
   free(simple_set->elements);
   free(simple_set);
}

/** free a shallow copy of a set of units */
void simple_set_shallow_free_units(
   SIMPLE_SET*           simple_set,         /**< set */
   int                   num_cols_x          /**< number of covariates */
   )
{

   free(simple_set->elements);
   free(simple_set);
}


/** find units with same covariate value starting from a given index */
int simple_set_next_shallow_split(
   const STRATEGY*       strategy,           /**< tree-building strategy */
   const SIMPLE_SET*     right_set,          /**< set */
   int                   p,                  /**< covariate to split on */
   int                   start,              /**< starting index */
   const double*         data_xp,            /**< values for covariate to split on */
   double*               splitval,           /**< (pointer to) found value to split on */
   ELEMENT**             elts,               /**< (pointer to) the elements moved */
   int*                  nelts,              /**< (pointer to) number of elements moved */
   int                   splitcount          /**< number of previous splits */
   )
{

   int idx;

   if( exploit_binaryvars(strategy) && (right_set->nkeyvals[p] == 2) )
   {
      /* do nothing: left rewards will be computed directly from full dataset */
      return (splitcount == 0);
   }

   idx = right_set->start + start;
   
   /* nothing to move from right to left */
   if( !(idx < right_set->start + right_set->n) )
      return 0;

   /* splitting value is just starting covariate value on the right */
   *splitval = data_xp[right_set->elements[idx]]; 

   /* find any additional units on right with *splitval as covariate value */
   while( idx < right_set->start + right_set->n && data_xp[right_set->elements[idx]] == *splitval )
      idx++;

   /* if all moved from right to left this is not a split */
   if( idx == right_set->start + right_set->n )
      return 0;

   /* record what's moved from right to left */
   *nelts = idx - (right_set->start + start);
   *elts = right_set->elements + right_set->start + start;

   return 1;
}

/** for each action find (total) reward if that action applied to each unit in a set of units */
void simple_set_find_nosplit_rewards(
   const SIMPLE_SET*     simple_set,         /**< simple set */
   int                   num_cols_y,         /**< number of actions */
   const double*         data_y,             /**< gammas, data_y+(d*num_rows) points to values for reward d */
   int                   num_rows,           /**< number of rows in the data */
   double*               nosplit_rewards     /**< space for computed no split rewards */
   )
{
   int d;
   int i;
   const double* dyelt;

   assert( simple_set != NULL );
   assert( nosplit_rewards != NULL );

   /* find reward for each action if no split were done */
   for( d = 0; d < num_cols_y; d++ )
      nosplit_rewards[d] = 0.0;
   for( i = simple_set->start; i < simple_set->start + simple_set->n; i++ )
   {
      dyelt = data_y + simple_set->elements[i];
      for( d = 0; d < num_cols_y; d++ )
      {
         nosplit_rewards[d] += *dyelt;
         dyelt += num_rows;
      }
   }
}

/** get number of distinct values of a covariate */
int nkeyvals(
   const SIMPLE_SET*     simple_set,         /**< simple set */
   int                   i                   /**< covariate */
   )
{
   return simple_set->nkeyvals[i];
}

/** update left rewards using single split on given covariate, return whether there is,
 * in fact, a split
 */
int simple_set_update_left_rewards_from_full(
   const SIMPLE_SET*     simple_set,         /**< units */
   int                   p,                  /**< covariate to split on */
   double*               left_rewards,       /**< rewards for each action for a 'left' set of units */
   const double*         data_xp,            /**< values for covariate to split on */
   const double*         data_y,             /**< data_y[d*num_rows+elt] is the reward for action d for unit elt */
   int                   num_rows,           /**< number of units in full dataset */
   int                   num_cols_y,         /**< number of actions */
   double*               splitval            /**< (pointer to) found value to split on */
   )
{

   int i;
   const KEY* keys = simple_set->keys[p];
   ELEMENT elt;
   int foundkey0 = 0;
   int foundkey1 = 0;
   int firstkey0 = -1;
   
   if( num_cols_y == 2 )
   {
      const double* data_y1 = data_y + num_rows;
      for( i = simple_set->start; i < simple_set->start + simple_set->n; i++)
      {
         elt = simple_set->elements[i];
         
         if( keys[elt] == 0 )
         {
            left_rewards[0] += data_y[elt];
            left_rewards[1] += data_y1[elt];
            if( foundkey0 == 0 )
            {
               foundkey0 = 1;
               firstkey0 = i;
            }
         }
         else
         {
            foundkey1 = 1;
         }
      }
   }
   else
   {
      int d;
      for( i = simple_set->start; i < simple_set->start + simple_set->n; i++)
      {
         elt = simple_set->elements[i];
         
         if( keys[elt] == 0 )
         {
            for( d = 0; d < num_cols_y; d++ )
               left_rewards[d] += data_y[d*num_rows+elt];
            if( foundkey0 == 0 )
            {
               foundkey0 = 1;
               firstkey0 = i;
            }
         }
         else
         {
            foundkey1 = 1;
         }
      }
   }

   if( foundkey0 && foundkey1 )
   {
      *splitval = data_xp[simple_set->elements[firstkey0]];
      return 1;
   }
   else
   {
      return 0;
   }
}
