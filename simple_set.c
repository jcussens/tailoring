#include "units.h"
#include "workspace.h"
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h> /* only for debugging */

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

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
   SORTED_SET*           sorted_set,         /**< sorted set */
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
   SORTED_SET*           sorted_set,         /**< sorted set */
   int                   elt                 /**< element to remove */
  )
{

   int left;
   int right;
   int mid;
   int midval;
   int eltval;
   const int* key = sorted_set->key;
   int n;
   int* set;  
   int retval = -1;
   
   assert(sorted_set != NULL);
   
   n = sorted_set->n;
   set = sorted_set->elements;  
   eltval = key[elt];

   left = 0;
   right = n-1;

   /* binary search for element */
   while( left <= right )
   {
      mid = (left+right)/2;
      midval = key[set[mid]];
      
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
   SORTED_SET*           sorted_set,         /**< sorted set */
   int                   nelts,              /**< number of elements to insert */ 
   const int*            elts                /**< elements to insert */ 
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
   SORTED_SET*           sorted_set,         /**< sorted set */
   int                   nelts,              /**< number of elements to insert */ 
   const int*            elts                /**< elements to insert */ 
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
   SORTED_SET*           sorted_set,         /**< sorted set */
   int                   nelts,              /**< number of elements to remove */ 
   const int*            elts                /**< elements to insert */ 
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
   SORTED_SET*           sorted_set,         /**< sorted set */
   int                   nelts               /**< number of elements to remove from start */ 
  )
{
   memmove(sorted_set->elements,sorted_set->elements+nelts,(sorted_set->n-nelts)*sizeof(int));
   sorted_set->n -= nelts;
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
   double*               splitval            /**< (pointer to) found value to split on */
   )
{
   SORTED_SET* left_sorted_setp;
   SORTED_SET* right_sorted_setp;
   int nmoved;
   int pp;
   const int* right_sorted_setp_elements;
   /* int i; */
   
   assert(left_sorted_sets != NULL);
   assert(right_sorted_sets != NULL);
   assert(p >= 0 && p < num_cols_x);
   assert(data_xp != NULL);
   assert(num_cols_x > 0);
   assert(splitval != NULL);
   assert((elts != NULL && nelts != NULL) || (elts == NULL && nelts == NULL));
   
   left_sorted_setp = left_sorted_sets[p];
   right_sorted_setp = right_sorted_sets[p];
   right_sorted_setp_elements = right_sorted_setp->elements;

   /* nothing to move from right to left */
   if( right_sorted_setp->n == 0 )
      return 0;

   /* splitting value is just first covariate value on the right */
   *splitval = data_xp[right_sorted_setp_elements[0]]; 

   /* if left is non-empty the last element on left must have a strictly lower pth covariate value
      that first element on right */
   assert(left_sorted_setp->n == 0 || data_xp[left_sorted_setp->elements[(left_sorted_setp->n)-1]] < *splitval); 

   /* find any additional units on right with *splitval as covariate value */
   nmoved = 1;
   while( nmoved < right_sorted_setp->n && data_xp[right_sorted_setp_elements[nmoved]] == *splitval )
      nmoved++;

   /* if all moved from right to left this is not a split */
   if( nmoved == right_sorted_setp->n )
      return 0;
   
   /* update left and right sorted sets of units for covariates other than splitting covariate */
   for( pp = 0; pp < num_cols_x; pp++)
   {
      if( pp != p )
      {
         insert_elements(left_sorted_sets[pp], nmoved, right_sorted_setp_elements);
         remove_elements(right_sorted_sets[pp], nmoved, right_sorted_setp_elements);
      }
   }

   /* update left and right sorted sets for splitting covariate */
   add_elements_at_end(left_sorted_setp, nmoved, right_sorted_setp_elements);
   remove_elements_at_start(right_sorted_setp, nmoved);
   
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
   memcpy(right->elements, simple_set->elements + simple_set->start, simple_set->n*sizeof(int));
   sort_elements(right->elements, right->n, right->keys[p], right->nkeyvals[p]); 
   
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
   const SORTED_SET**    right_sorted_sets,  /**< sorted set */
   int                   p,                  /**< covariate to split on */
   int                   start,              /**< starting index */
   const double*         data_xp,            /**< values for covariate to split on */
   double*               splitval,           /**< (pointer to) found value to split on */
   int**                 elts,               /**< (pointer to) the elements moved */
   int*                  nelts               /**< (pointer to) number of elements moved */
   )
{
   const SORTED_SET* right_sorted_setp = right_sorted_sets[p];
   int idx;
   
   /* nothing to move from right to left */
   if( !(start < right_sorted_setp->n) )
      return 0;

   /* splitting value is just starting covariate value on the right */
   *splitval = data_xp[right_sorted_sets[p]->elements[start]]; 

   /* find any additional units on right with *splitval as covariate value */
   idx = start + 1;
   while( idx < right_sorted_setp->n && data_xp[right_sorted_setp->elements[idx]] == *splitval )
      idx++;

   /* if all moved from right to left this is not a split */
   if( idx == right_sorted_setp->n )
      return 0;

   /* record what's removed */
   *nelts = idx - start;
   *elts = (right_sorted_setp->elements) + start;

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

