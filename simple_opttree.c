#include "simple_opttree.h"
#include "sorted_set.h"


/**
 * For each unit find and record the best and worst actions for that unit
 */
static
void store_best_worst_actions(
   const double*         data_y,             /**< rewards (column major) */
   int                   num_rows,           /**< number of units */
   int                   num_cols_y,         /**< number of rewards */
   int**                 best_actions,       /**< (pointer to) best actions */
   int**                 worst_actions       /**< (pointer to) worst actions */
  )
{
   int i;
   int best_action;
   int worst_action;
   int d;
   double best_reward;
   double worst_reward;

   assert(num_cols_y > 1);

   *best_actions = (int*) malloc(num_rows*sizeof(int));
   *worst_actions = (int*) malloc(num_rows*sizeof(int));  

   for( i = 0; i < num_rows; i++ )
   {
      /* there are always at least two actions
         initialise best to action 0
         and worst to action 1, so that even
         if we get equal rewards the best and worst action must differ
         (it is just simpler to assume that best and worst differ
      */
      
      best_action = 0;
      best_reward = data_y[i];
      worst_action = 1;
      worst_reward = data_y[num_rows+i];
    
      for( d = 0; d < num_cols_y; d++ )
      {
         if( data_y[d*num_rows+i] > best_reward )
         {
            best_action = d;
            best_reward = data_y[d*num_rows+i];
         }
         if( data_y[d*num_rows+i] < worst_reward )
         {
            worst_action = d;
            worst_reward = data_y[d*num_rows+i];
         }
      }
      (*best_actions)[i] = best_action;
      (*worst_actions)[i] = worst_action;
      assert(worst_action != best_action);
      assert(worst_reward <= best_reward);
   }
}


/** Find an optimal depth=1 (ie at most one single split) tree for a set of units 
 */
static
void level_one_learning(
   NODE*                 node,               /**< uninitialised tree to be populated with optimal tree */
   SORTED_SET**          sorted_sets,        /**< sorted sets for the units, one for each covariate */
   const double*         data_x,             /**< covariates, data_x+(j*num_rows) points to values for covariate j */
   const double*         data_y,             /**< gammas, data_y+(d*num_rows) points to values for reward d */
   int                   num_rows,           /**< number of units in full dataset */
   int                   num_cols_x,         /**< number of covariates */
   int                   num_cols_y,         /**< number of rewards/actions */
   const int*            best_actions,       /**< best_actions[i] is the best action for unit i */
   const int*            worst_actions,      /**< worst_actions[i] is the worst action for unit i */
   WORKSPACE*            workspace,          /**< workspace */
   int*                  perfect             /**< *perfect=1 iff each unit assigned its best action, else *perfect=0 */
   )
{
   int p;
   double best_reward;
   int first_reward = 1;
   double* nosplit_rewards = get_rewards(workspace);

   /* get reward for each action if no split were done */
   find_nosplit_rewards(sorted_sets, data_y, num_rows, nosplit_rewards);
   
   /* consider each covariate for splitting */
   for( p = 0; !(*perfect) && p < num_cols_x; p++)
   {
      const double* data_xp = data_x+(p*num_rows);
      const SORTED_SET* sorted_set = sorted_sets[p];
      SORTED_SET* left_sorted_set = light_emtpy_copy(sorted_set);
      SORTED_SET* right_sorted_set = light_full_copy(sorted_set);

      double this_reward;

      int left_perfect;
      int right_perfect;
      
      int* elts;
      int nelts;

      /* consider each split (x[p] <= splitval, x[p] > splitval) of the data */
      while( !(*perfect) && next_light_split(left_sorted_set, right_sorted_set, p, data_xp, &splitval, elts, &nelts) )
      {
         
         /* find optimal tree for left data set */
         find_best_split(left_child(node), depth-1, left_sorted_sets, min_node_size, data_x, data_y,
            num_rows, num_cols_x, num_cols_y, best_actions, worst_actions, workspace, &left_perfect); 

         /* find optimal tree for right data set */
         find_best_split(right_child(node), depth-1, right_sorted_sets, min_node_size, data_x, data_y,
            num_rows, num_cols_x, num_cols_y, best_actions, worst_actions, workspace, &right_perfect); 

         /* tree is perfect if and only if both left and right tree are perfect */
         *perfect = left_perfect && right_perfect;

         /* get reward for this split */
         this_reward = get_reward(left_child(node)) + get_reward(right_child(node));

         /* if best so far, update */
         if( first_reward || this_reward > best_reward ) 
         {
            best_reward = this_reward;
            record_split(node, p, splitval, best_reward);
            record_best_tree(workspace, node, depth);

            first_reward = 0;
         }
      }
   }

   /* set node to best tree */
   retrieve_best_tree(workspace, node, depth);
   
}

/** on return, `node` will be the root of an optimal tree of depth `depth` for the data
 * represented by `sorted_sets`
 */
static
void find_best_split(
   NODE*                 node,               /**< uninitialised tree to be populated with optimal tree */
   int                   depth,              /**< depth of tree */
   const SORTED_SET**    sorted_sets,        /**< sorted sets for the units, one for each covariate */
   int                   min_node_size,      /**< smallest terminal node size */
   const double*         data_x,             /**< covariates, data_x+(j*num_rows) points to values for covariate j */
   const double*         data_y,             /**< gammas, data_y+(d*num_rows) points to values for reward d */
   int                   num_rows,           /**< number of units in full dataset */
   nt                    num_cols_x,         /**< number of covariates */
   int                   num_cols_y,         /**< number of rewards/actions */
   const int*            best_actions,       /**< best_actions[i] is the best action for unit i */
   const int*            worst_actions,      /**< worst_actions[i] is the worst action for unit i */
   WORKSPACE*            workspace,          /**< working space */
   int*                  perfect             /**< *perfect=1 iff each unit assigned its best action, else *perfect=0 */
  )
{

   int p;
   int pure;
   double best_reward;
   int first_reward = 1;
   
   /* determine whether the dataset is pure, i.e. whether each unit has same best action */
   pure = is_pure(sorted_sets, best_actions);

   /* if this dataset is pure we can make a perfect leaf */
   if( pure )
      *perfect = 1;

   /* nothing further to do if depth limit reached or if too few datapoints for splitting or if dataset is pure */
   if( depth == 0 || get_size(sorted_sets) <= min_node_size || pure )
   {
      int best_action;

      /* find best action and its associated reward */
      find_best_reward(sorted_sets, data_y, num_rows, num_cols_y, workspace, &best_reward, &best_action);

      /* make node a leaf with found best action and associated reward */
      make_leaf(node, best_reward, best_action);

      /* all done */
      return;
   }

   /* use specialised routine if only one split allowed */
   if( depth == 1 )
   {
      level_one_learning(node, sorted_sets, data_x, data_y,
         num_rows, num_cols_x, num_cols_y, best_actions, worst_actions, temp_space, perfect); 
      
      return;
   }

   /* consider each covariate for splitting */
   for( p = 0; !(*perfect) && p < num_cols_x; p++)
   {
      const double* data_xp = data_x+(p*num_rows);

      SORTED_SET** left_sorted_sets;
      SORTED_SET** right_sorted_sets;

      double this_reward;

      int left_perfect;
      int right_perfect;

      int* elts;
      int nelts;
      
      /* initialise so that each left_sorted_set (for each covariate) is empty and each 
         right_sorted_set is a copy of the input sorted set (for that covariate) */
      initialise_sorted_sets(sorted_sets, depth, num_cols_x, workspace, &left_sorted_sets, &right_sorted_sets);

      /* consider each split (x[p] <= splitval, x[p] > splitval) of the data */
      while( !(*perfect) && next_split(left_sorted_sets, right_sorted_sets, p, data_xp, num_cols_x, workspace, 
            &splitval, &elts, &nelts) )
      {
         
         /* find optimal tree for left data set */
         find_best_split(left_child(node), depth-1, left_sorted_sets, min_node_size, data_x, data_y,
            num_rows, num_cols_x, num_cols_y, best_actions, worst_actions, workspace, &left_perfect); 

         /* find optimal tree for right data set */
         find_best_split(right_child(node), depth-1, right_sorted_sets, min_node_size, data_x, data_y,
            num_rows, num_cols_x, num_cols_y, best_actions, worst_actions, workspace, &right_perfect); 

         /* tree is perfect if and only if both left and right tree are perfect */
         *perfect = left_perfect && right_perfect;

         /* get reward for this split */
         this_reward = get_reward(left_child(node)) + get_reward(right_child(node));

         /* if best so far, update */
         if( first_reward || this_reward > best_reward ) 
         {
            best_reward = this_reward;
            record_split(node, p, splitval, best_reward);
            record_best_tree(workspace, node, depth);

            first_reward = 0;
         }
      }
   }

   /* set node to best tree */
   retrieve_best_tree(workspace, node, depth);
}



/** Find an optimal policy tree of given maximal depth 
 * @return An optimal policy tree
 */
NODE* tree_search_simple(
  int                    depth,              /**< (maximum) depth of returned tree */
  int                    min_node_size,      /**< smallest terminal node size */
  const double*          data_x,             /**< covariates, data_x+(j*num_rows) points to values for covariate j */
  const double*          data_y,             /**< gammas, data_y+(d*num_rows) points to values for reward d */
  int                    num_rows,           /**< number of units in full dataset */
  int                    num_cols_x,         /**< number of covariates */
  int                    num_cols_y          /**< number of rewards/actions */
  )
{
   NODE* tree;
   SORTED_SET** initial_sorted_sets;
   WORKSPACE* workspace;
   int* best_actions;
   int* worst_actions;
   int perfect;
   
   /* make tree of right depth, values for its nodes will be overwritten with those defining an optimal tree */
   tree = make_tree(depth);

   /* for each of the num_cols_x covariates, use the values of that covariate to sort the
    * (indices of the) datapoints and store the result   
    */
   initial_sorted_sets = make_initial_sorted_sets(data_x, num_rows, num_cols_x);

   /* create working spaces of various sorts (trees, sorted sets, arrays of rewards, etc) */
   workspace = make_workspace(depth, initial_sorted_sets, num_rows, num_cols_x, num_cols_y);

   /* compute and store best and worst actions for each unit */
   store_best_worst_actions(data_y, num_rows, num_cols_y, &best_actions, &worst_actions);

   /* record that we have yet to find a 'perfect' tree for this data */
   perfect = 0;

   /* find the optimal tree */
   find_best_split(tree, depth, initial_sorted_sets, min_node_size, data_x, data_y,
      num_rows, num_cols_x, num_cols_y, best_actions, worst_actions, workspace, &perfect); 

   /* free memory */
   free_workspace(workspace, depth, num_cols_x);
   free_sorted_sets(initial_sorted_sets, num_cols_x);

   /* remove any nodes below leaves, and merge leaves with the same action */
   fix_tree(tree);

   /* return optimal tree */
   return tree;
}
