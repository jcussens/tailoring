#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include "opttree.h"
#define MAX_LINE_LENGTH 2000
#define MAX_COVARIATE_NAME_LENGTH 50

/** count number of space separate fields in a line 
 * @return The number of space separated fields
 */
static
int getnfields(
   char*                 line                /**< line */
  )
{
  int i;
  int nfields = 0;
  
  for( i = 0; i < strlen(line); i++ )
  {
    if( !isspace(line[i]) && (line[i+1] == '\0' || isspace(line[i+1]) ) )
      nfields++;
  }
  return nfields;
}


/** prune tree: if a subtree has the same action for all leaves, replace with a
 * single leaf with that action 
*/
static
void prune_tree(
   NODE*                 root                /**< root node */
   )
{
   if( root->left_child != NULL )
      prune_tree( root->left_child );

   if( root->right_child != NULL )
      prune_tree( root->right_child );

   /* just delete dummy nodes */
   if( root->left_child != NULL && root->left_child->index == -1
      && root->left_child->action_id == -1 )
   {
      free(root->left_child);
      root->left_child = NULL;
   }

   if( root->right_child != NULL && root->right_child->index == -1
      && root->right_child->action_id == -1 )
   {
      free(root->right_child);
      root->right_child = NULL;
   }

   if( root->left_child != NULL  && root->right_child != NULL 
      && root->left_child->index == -1  && root->right_child->index == -1
      && root->left_child->action_id == root->right_child->action_id )
   {
      root->index = -1;
      root->reward = root->left_child->reward + root->right_child->reward;
      root->action_id = root->left_child->action_id;
      free(root->left_child);
      free(root->right_child);
      root->left_child = NULL;
      root->right_child = NULL;
   }
   
   
   
   
}
/** breadth-first index for a node, return -1 if not found */
static
int bfidx(
   NODE*                 root,                /**< root node */
   NODE*                 node,                /**< node to find bf index for */
   int                   count                /**< how much to add to index */
   )
{

   int idx;
   
   if( node == root )
      return count;
   else if( root->left_child != NULL )
   {
      idx = bfidx(root->left_child, node, count+1);
      if( idx != -1)
         return idx;
   }
   else if( root->right_child != NULL )
   {
      idx = bfidx(root->right_child, node, count+2);
      if( idx != -1)
         return idx;
   }
   else
   {
      return -1;
   }
}

/** get number of nodes in a tree */
static
int nnodes(
   NODE*                 tree                /**< policy tree to print */
   )
{
   if( tree->index != -1)
      return 1 + nnodes(tree->left_child) + nnodes(tree->right_child);
   else
      return 1;
}

/** prints a policy tree in policytree style to standard output */
static
void print_tree_policytree_rec(
   NODE*                 tree,               /**< policy tree to print */
   char**                covnames,           /**< covnames[j] is the name of covariate j */
   int                   level               /**< level of the tree */
  )
{

   /* indent */
   int i;

   for(i = 0; i < level; i++)
      printf("  ");
   
   if( tree->index != -1)
   {
      printf("split_variable: %s  split_value: %g\n", covnames[tree->index], tree->value);
      print_tree_policytree_rec(tree->left_child, covnames, level+1);
      print_tree_policytree_rec(tree->right_child, covnames, level+1);
   }
   else
   {
      printf("* action: %d\n", tree->action_id + 1);
   }
}

/** prints a policy tree in policytree style to standard output */
static
void print_tree_policytree(
   NODE*                 tree,               /**< policy tree to print */
   char**                covnames,           /**< covnames[j] is the name of covariate j */
   int                   depth,              /**< depth of the tree */
   int                   nactions,           /**< number of actions */
   char**                actionnames         /**< covnames[j] is the name of action j */
  )
{
   int i;
   
   printf("policy_tree object\n");
   printf("Tree depth:   %d\n", depth);
   printf("Actions: ");
   for(i = 0; i < nactions; i++)
      printf(" %d: %s", i+1, actionnames[i]);
   printf("\nVariable splits:\n");

   print_tree_policytree_rec(tree, covnames, 0);   
}

/** prints a policy tree to standard output */
static
void print_tree(
   NODE*                 tree,               /**< policy tree to print */
   char**                covnames            /**< covnames[j] is the name of covariate j */
  )
{
  printf("node = %p\n", (void*) tree);
  printf("reward = %g\n", tree->reward);
  if( tree->index != -1)
  {
    printf("covariate = %s\n", covnames[tree->index]);
    printf("value = %g\n", tree->value);
    printf("left_child = %p\n", (void*) tree->left_child);
    printf("right_child = %p\n", (void*) tree->right_child);
  }
  else
  {
     printf("action_id = %d\n", tree->action_id);
  }
  printf("\n");

  if( tree->index != -1)
  {
    print_tree(tree->left_child,covnames);
    print_tree(tree->right_child,covnames);
  }
}


/** 
 * call like ./a.out filename nactions [depth]
 * file should have a single header line and then each subsequent line is covariates followed by nactions reward values.
 * assume white space separator with no white space in names of covariates
 */
int main(
   int                   argc,               /**< number of command line arguments */
   char*                 argv[]              /**< command line arguments */
  )
{

  FILE* file;
  int depth = 3;            /** (maximum) depth of returned tree */
  int split_step = 0;       /** consider splits every split_step'th possible split */
  int min_node_size = 1;    /** smallest terminal node size */
  int num_rows;
  int num_cols_x;
  int num_cols_y;
  double* data_x;
  double* data_y;
  int status;
  int row;
  int col;

  char line[MAX_LINE_LENGTH];
  char* lineptr;
  char tmpstr[MAX_COVARIATE_NAME_LENGTH];
  int nfields;
  int i;
  char** covnames;
  char** actionnames;
  int nlines;
  int offset;

  NODE* tree;
  int method;

  
  if( argc < 3 )
  {
    printf("Need to supply at least a filename, the number of actions and method indicator.\n");
    return 1;
  }

  file = fopen(argv[1],"r");

  if( file == NULL )
  {
    printf("Could not open %s.\n", argv[1]);
    return 1;
  }    

  num_cols_y = atoi(argv[2]);

  method = atoi(argv[3]);
  if(method != 1 && method != 2 )
  {
    printf("Incorrect method number: %d, should be 1 or 2.\n", method);
    return 1;
  }    

  if( argc > 4)
    depth = atoi(argv[4]);

  
  /* count lines in file, ignoring any lines composed entirely of white space*/
  nlines = 0;
  while( fgets(line,MAX_LINE_LENGTH,file) != NULL)
  {
    for(i = 0; line[i] != '\0'; i++)
    {
      if( !isspace(line[i]) )
        {
          if( nlines == 0 )
            nfields = getnfields(line);
          else if(nfields != getnfields(line))
          {
            printf("Error: line %d has %d fields but should have %d\n", nlines+1, getnfields(line), nfields);
            return 1;
          }
          nlines++;
          break;
        }
    }
  }
  rewind(file);
  num_rows = nlines - 1;
  
  /* read in header */

  if( fgets(line,MAX_LINE_LENGTH,file) == NULL )
  {
    printf("Could not read header line.\n");
    return 1;
  }
  
  /* get number of names in header line */
  nfields = 0;
  for( i = 0; i < strlen(line); i++ )
  {
    if( !isspace(line[i]) && (line[i+1] == '\0' || isspace(line[i+1]) ) )
      nfields++;
  }
  /* assume names for actions also in header line, but don't read them */
  num_cols_x = nfields - num_cols_y;

  /* read in covariate names */
  covnames = (char**) malloc(num_cols_x*sizeof(char*));
  actionnames = (char**) malloc(num_cols_y*sizeof(char*));
  lineptr = line;
  for(i = 0; i < num_cols_x+num_cols_y; i++)
  {
    status = sscanf(lineptr, "%s%n", tmpstr, &offset);
    if( status != 1 )
    {
      if( i < num_cols_x )
        printf("Error reading covariate name with index %d from %s.\n", i, line);
      else
        printf("Error reading action name with index %d from %s.\n", i-num_cols_x, line);
      return 1;
    }

    /* add an 'X' if name starts with a digit (to be like policytree) */
    if( isdigit(*tmpstr) )
    {
       memmove(tmpstr+1,tmpstr,(strlen(tmpstr)+1));
       tmpstr[0] = 'X';
    }
    
    if( i < num_cols_x )
    {
      covnames[i] = (char*) malloc((strlen(tmpstr)+1)*sizeof(char));
      strcpy(covnames[i],tmpstr);
    }
    else
    {
      actionnames[i-num_cols_x] = (char*) malloc((strlen(tmpstr)+1)*sizeof(char));
      strcpy(actionnames[i-num_cols_x],tmpstr);
    }
    lineptr += offset;
  }

  data_x = (double *) malloc(num_rows*num_cols_x * sizeof(double));
  data_y = (double *) malloc(num_rows*num_cols_y * sizeof(double));

  for(row = 0; row < num_rows; row++)
  {
    for(col = 0; col < num_cols_x; col++)
    {
      status = fscanf(file, "%lg", data_x+col*num_rows+row);
      if( status != 1)
      {
        printf("Error reading in value for covariate %s with index %d for datapoint %d.\n",covnames[col],col,row);
        return 1;
      }
    }
    for(col = 0; col < num_cols_y; col++)
    {
      status = fscanf(file, "%lg", data_y+col*num_rows+row);
      if( status != 1)
      {
        printf("Error reading in value for action %d for datapoint %d.\n",col,row);
        return 1;
      }
    }
  }

  fclose(file);

  if( method == 1 )
     tree = tree_search_jc_policytree(depth, split_step, min_node_size, data_x, data_y, num_rows, num_cols_x, num_cols_y);
  else
     tree = tree_search_jc_discretedata(depth, split_step, min_node_size, data_x, data_y, num_rows, num_cols_x, num_cols_y);

  /* printf("Actions: "); */
  /* for(i = 0; i < num_cols_y; i++) */
  /*   printf("%d: %s ",i,actionnames[i]); */
  /* printf("\n"); */
  
  /* print_tree(tree,covnames); */

  prune_tree(tree);
  print_tree_policytree(tree, covnames, depth, num_cols_y, actionnames);
  
  for(i = 0; i < num_cols_x; i++)
    free(covnames[i]);
  free(covnames);

  for(i = 0; i < num_cols_y; i++)
    free(actionnames[i]);
  free(actionnames);

  tree_free(tree);
  
  free(data_x);
  free(data_y);

  return 1-status;
}
