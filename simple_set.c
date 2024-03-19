#include "units.h"
#include "workspace.h"
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h> /* only for debugging */

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

#define MAXKEYVALS 30

/** 
 * simple set
*/
struct simple_set
{
   int*                  elements;           /**< pointer to space for the set */
   int                   size;               /**< size of space for set */
   int                   start;              /**< elements[start] is the first element of the set */ 
   int                   n;                  /**< number of elements in set */
   int**                 keys;               /**< let data_xp be covariate values for covariate p then 
                                              * (1) if data_xp[i] < data_xp[j] then keys[p][i] < keys[p][j] and
                                              * (2) if data_xp[i] == data_xp[j] then keys[p][i] == keys[p][j] */
   int*                  nkeyvals;           /** nkeyvals[p] is the number of distinct data_xp values */
};

static
void counting_sort(
   const int*            a,                  /**< input array to sort */
   int                   n,                  /**< length of input array */
   int*                  b,                  /**< will contain sorted array */
   int                   k,                  /**< number of distinct key values */
   const int*            key,                /**< key[elem] is the key value for elem */
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
      count[key[a[i]]]++;

   /* alter so that c[i] is number of key vals <= i */ 
   for(i = 1; i < k; i++ )
      count[i] += count[i-1];

   /* put each element of a into correct place in b */
   for(i = n-1; i >=0; i--)
      b[(count[key[a[i]]]--)-1] = a[i];
}


static
void counting_sort_radix(
   const int*            a,                  /**< input array to sort */
   int                   n,                  /**< length of input array */
   int*                  b,                  /**< will contain sorted array */
   const int*            key,                /**< key[elem] is the key value for elem */
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
      count[(key[a[i]]/exp)%10]++;

   /* alter so that c[i] is number of key vals <= i */ 
   for(i = 1; i < 10; i++ )
      count[i] += count[i-1];

   /* put each element of a into correct place in b */
   for(i = n-1; i >=0; i--)
      b[(count[(key[a[i]]/exp)%10]--)-1] = a[i];
}

static
void radix_sort(
   const int*            a,                  /**< input array to sort */
   int                   n,                  /**< length of input array */
   int*                  b,                  /**< will contain sorted array */
   int*                  tmp,                /**< tmp array, at least as long as b */
   const int*            key,                /**< key[elem] is the key value for elem */
   int                   maxkey              /**< maximal value of key[elem] */
   )
{
   int exp;

   /* sort a on least significant digit and put result in b */
   counting_sort_radix(a, n, b, key, 1);
   
   for( exp = 10; exp / maxkey > 0; exp *= 10 )
   {
      counting_sort_radix((const int*) b, n, tmp, key, exp);
      b = tmp;
   }
}

/** sort an array of integers according to key values */
static
void sort_units(
   const int*            a,                  /**< array to be sorted */
   int                   n,                  /**< length of array to be sorted */
   const int*            key,                /**< key[elem] is key value for elem */
   int                   nkeyvals,           /**< key values are 0,...,nkeyvals-1 */
   int*                  tmp,                /**< temporary array of length at least n */
   int*                  b                   /**< output array of length at least n */
   )
{
   if( nkeyvals <= MAXKEYVALS )
      counting_sort(a, n, b, nkeyvals, key);
   else
      radix_sort(a, n, b, tmp, key, nkeyvals-1);
}

/** comparison function for sorting integers in ascending order */
static
int cmpfunc (
   const void* a,
   const void* b
   )
{
   return ( *(int*)a - *(int*)b );
}


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

/** merge a pair of ordered subsequences to get a single ordered subsequence
 * left run is indices[ileft:iright-1]
 * right run is indices[iright:iend-1]
 */
static
void bottomupmergeint(
   const int*            indices,            /**< global array */
   int                   ileft,              /**< index of first element of left subsequence */
   int                   iright,             /**< index of first element of right subsequence */
   int                   iend,               /**< iend-1 is index of last element of right subsequence */
   int*                  tmp,                /**< on output tmp[ileft:iend-1] has the merged subsequence */
   const int*            key                 /**< key[i] is ordering key for i */
  )
{
  int i = ileft;
  int j = iright;
  int k;

  for( k = ileft; k < iend; k++ )
  {
    if( i < iright && (j >= iend || key[indices[i]] <= key[indices[j]]) )
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

/** sort an array by bottom up merge sort */
static
void bottomupmergesortint(
   int*                  indices,            /**< array to sort */
   int                   n,                  /**< length of array */
   const key*            key                 /**< key[i] is ordering key for i */
  )
{
  int width;
  int i;
  const size_t size = n*sizeof(int);
  /* TODO: preallocate this (or use a different sorting algorithm) */
  int* tmp = (int*) malloc(size);
  
  for( width = 1; width < n; width *= 2)
  {
    for( i = 0; i < n; i += 2*width)
      bottomupmergeint(indices, i, MIN(i+width,n), MIN(i+2*width,n), tmp, key);

    memcpy(indices,tmp,size);
  }
  free(tmp);
}



int units_ok(
   const SIMPLE_SET*     simple_set,         /**< set */
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
      int p;
      
      if( simple_set->keys == NULL )
         return 0;
      if( simple_set->nkeyvals == NULL )
         return 0;

      for( p = 0; p < num_cols_x; p++ )
      {
         int i;
         
         if(simple_set->keys[p] == NULL )
            return 0;
         if(simple_set->nkeyvals[p] < 1 )
            return 0;

         for( i = 0; i < num_rows; i++ )
         {
            if( !(simple_set->keys[p][i] < simple_set->nkeyvals[p]) )
               return 0;
         }
      }
   }
   
   return 1;
}
   
/** Determine whether a set is 'pure'.
 * A pure set is one where each unit has the same best action
 * @return 1 if the set is pure, else 0
 */
int is_pure(
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

/** get common size of sorted sets */
int get_size(
   const SIMPLE_SET*     simple_set          /**< set */
   )
{
   assert(simple_set != NULL);
   
   return simple_set->n;
}

/* find best action and its associated reward for a set of units */
void find_best_action(
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

/** given that left_set and right_set are sorted according to covariate p, find next split and associated split value
 * @return 1 if there is a next split, otherwise 0
 */
int next_split(
   SIMPLE_SET*           left_set,           /**< left set sorted on covariate p */
   SIMPLE_SET*           right_set,          /**< right set sorted on covariate p */ 
   int                   p,                  /**< covariate to split on */
   const double*         data_xp,            /**< values for covariate to split on */
   int                   num_cols_x,         /**< number of covariates */
   double*               splitval            /**< (pointer to) found value to split on */
   )
{
   int leftend;

   assert(left_set != NULL);
   assert(right_set != NULL);
   assert(p >= 0 && p < num_cols_x);
   assert(data_xp != NULL);
   assert(num_cols_x > 0);
   assert(splitval != NULL);

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
   while( right_set->n > 0 && data_xp[right_set->_elements[right_set->start]] == *splitval )
   {
      left->elements[leftend++] = right_set->elements[right_set->start++];
      left->n++;
      right->n--;
   }

   /* if all moved from right to left this is not a split */
   if( right_set->n == 0 )
      return 0;
   
   return 1;
   
}

/* make a 'shallow' copy of source sorted sets */
SIMPLE_SET* shallow_copy_units(
   const SIMPLE_SET*     source,            /**< source  */
   int                   num_cols_x         /**< number of covariates */
   )
{
   SIMPLE_SET* target;
   int i;

   target = (SIMPLE_SET*) malloc(sizeof(SIMPLE_SET));
   target->elements = (int*) malloc(source->size * sizeof(int));
   memcpy(target->elements, source->elements, (source->size)*sizeof(int));
   target->size = source->size;
   target->start = source->start;
   target->n = source->n;
   target->keys = source->keys;
   target->nkeyvsals = source->nkeyvals;

   return target;
}


static
int* get_key(
   int*                  elements,           /**< elements (unordered) */
   const double*         data_xp,            /**< covariate values for some covariate */
   int                   num_rows,           /**< number of units in full dataset */
   int*                  tmp,                /**< temporary working space */
   int*                  nkeyvals            /**< number of distinct data_xp values */
   )
{
   int i;
   int* keysp;
   int keyval;
   
   /* sort the elements using data_xp values as key */
   bottomupmergesort(elements, tmp, num_rows, data_xp);

   keysp = (int*) malloc(num_rows * sizeof(int));

   keyval = 0;
   keysp[0] = keyval;
   *nkeyvals = 1;

   for(i = 1; i < num_rows; i++ )
   {
      if( data_xp[elements[i]] == data_xp[elements[i-1]] )
      {
         keysp[i] = keyval;
      }
      else
      {
         keysp[i] = ++keyval;
         (*nkeyvals)++;
      }
   }
   return keysp;
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

  if( data_xx != NULL )
  {
     /* sort the elements using data_xx values as key */
     bottomupmergesort(elements, tmp, num_indices, data_xx);
  }
  
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
void initialise_units(
   const SIMPLE_SET*     simple_set,         /**< input set */
   int                   p,                  /**< splitting covariate */
   int                   depth,              /**< depth of associated node */
   int                   num_cols_x,         /**< number of covariates */
   WORKSPACE*            workspace,          /**< workspace */
   SIMPLE_SET**          left_simple_set,    /**< pointer to output left set */
   SIMPLE_SET**          right_simple_set    /**< pointer to output right set */
   )
{
   int pp;
   SIMPLE_SET* left = NULL;
   SIMPLE_SET* right = NULL;
   
   left = get_left_sorted_sets(workspace,depth);
   right = get_right_sorted_sets(workspace,depth);

   assert( left != NULL );
   assert( right != NULL );

   /* make left empty and ensure any additions added at start */
   left->n = 0;
   left->start = 0;

   /* make right a copy of input, but with elements ordered according to covariate p */
   right->n = simple_set->n;
   right->start = 0;
   sort_units(simple_set->elements + simple_set->start, simple_set->n, simple_set->keys[p], simple_set->nkeyvals[p],
      get_tmpunits(workspace)->elements, right->elements);
   
   *left_simple_set = left;
   *right_simple_set = right;
}

SIMPLE_SET* make_units(
   const double*         data_x,             /**< covariates, data_x+(j*num_rows) points to values for covariate j */
   int                   num_rows,           /**< number of units in full dataset */
   int                   num_cols_x          /**< number of covariates */
   )
{

   int* tmp_indices;
   int p;
   SIMPLE_SET* initial_simple_set; 
   
   tmp_indices = (int*) malloc(num_rows*sizeof(int));
   
   initial_simple_set = (SIMPLE_SET*) malloc(sizeof(SIMPLE_SET));
   initial_simple_set->elements = (int*) malloc(num_rows * sizeof(int));
   for( i = 0; i < num_rows: i++ )
      initial_simple_set->elements[i] = i;
   initial_simple_set->size = num_rows;
   initial_simple_set->start = 0;
   initial_simple_set->n = num_rows;
   initial_simple_set->keys = (int**) malloc(num_cols_x * sizeof(int*));
   initial_simple_set->nkeyvals = (int*) malloc(num_cols_x * sizeof(int));
   for( p = 0; p < num_cols_x; p++)
   {
      initial_simple_set->keys[p] = get_key(initial_simple_set->elements, data_x+p*num_rows, num_rows, tmp_indices, &nkeyvals);
      initial_simple_set->nkeyvals[p] = nkeyvals;
   }
   
   free(tmp_indices);

   return initial_simple_set;
}

void free_units(
   SIMPLE_SET*           simple_set,         /**< set */
   int                   num_cols_x          /**< number of covariates */
   )
{
   int p;

   for( p = 0; p < num_cols_x; p++)
   {
      free(simple_set->keys[p]);
      free(simple_set->nkeyvals[p]);
   }
   free(simple_set->keys);
   free(simple_set->nkeyvals);
   free(simple_set->elements);
   free(simple_set);
}

void shallow_free_units(
   SIMPLE_SET*           simple_set,         /**< set */
   int                   num_cols_x          /**< number of covariates */
   )
{

   free(simple_set->elements);
   free(simple_set);
}


/** find units with same covariate value starting from a given index */
int next_shallow_split(
   const SIMPLE_SET*     right_set,          /**< set */
   int                   p,                  /**< covariate to split on */
   int                   start,              /**< starting index */
   const double*         data_xp,            /**< values for covariate to split on */
   double*               splitval,           /**< (pointer to) found value to split on */
   int**                 elts,               /**< (pointer to) the elements moved */
   int*                  nelts               /**< (pointer to) number of elements moved */
   )
{

   int idx = right_set->start + start;
   
   /* nothing to move from right to left */
   if( !(idx < right_set->start + right_set->n) )
      return 0;

   /* splitting value is just starting covariate value on the right */
   *splitval = data_xp[right_set->elements[idx]]; 

   /* find any additional units on right with *splitval as covariate value */
   while( idx < right_set->start + right_set->n && data_xp[right_set->elements[idx]] == *splitval )
      idx++

   /* if all moved from right to left this is not a split */
   if( idx == right_set->start + right_set->n )
      return 0;

   /* record what's moved from right to left */
   *nelts = idx - (right_set->start + start);
   *elts = right_set->elements + right_set->start + start;

   return 1;
}

/** for each action find (total) reward if that action applied to each unit in a set of units */
void find_nosplit_rewards(
   const SIMPLE_SET*     simple_set,         /**< sorted sets */
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

