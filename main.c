#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include "opttree.h"
#define MAX_LINE_LENGTH 2000
#define MAX_COVARIATE_NAME_LENGTH 50


static
int getnfields(
  char* line
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


static
void print_tree(
  NODE* tree,
  char** covnames
  )
{
  printf("node = %p\n", (void*) tree);
  if( tree-> index != -1)
    printf("covariate = %s\n", covnames[tree->index]);
  printf("value = %g\n", tree->value);
  printf("reward = %g\n", tree->reward);
  printf("action_id = %d\n", tree->action_id);
  printf("left_child = %p\n", (void*) tree->left_child);
  printf("right_child = %p\n", (void*) tree->right_child);
  printf("\n");

  if(tree->left_child != NULL)
    print_tree(tree->left_child,covnames);

  if(tree->right_child != NULL)
    print_tree(tree->right_child,covnames);
}


/** 
 * call like ./a.out filename nactions [depth]
 * file should have a single header line and then each subsequent line is covariates followed by nactions reward values.
 * assume white space separator with no white space in names of covariates
 */
int main(
  int argc,
  char* argv[]
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
  
  if( argc < 2 )
  {
    printf("Need to supply at least a filename and the number of actions.\n");
    return 1;
  }

  file = fopen(argv[1],"r");

  if( file == NULL )
  {
    printf("Could not open %s.\n", argv[1]);
    return 1;
  }    

  num_cols_y = atoi(argv[2]);
  
  if( argc > 3)
    depth = atoi(argv[3]);

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
  
  tree = tree_search(depth, split_step, min_node_size, data_x, data_y, num_rows, num_cols_x, num_cols_y);

  printf("Actions: ");
  for(i = 0; i < num_cols_y; i++)
    printf("%d: %s ",i,actionnames[i]);
  printf("\n");
  
  print_tree(tree,covnames);

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
