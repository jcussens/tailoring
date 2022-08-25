#include <stdio.h>
#include <stdlib.h>
#include "opttree.h"

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
  
  if( argc == 1 )
  {
    printf("Need to supply at least a filename.\n");
    return 1;
  }

  file = fopen(argv[1],"r");

  if( file == NULL )
  {
    printf("Could not open %s.\n", argv[1]);
    return 1;
  }    

  if( argc > 2)
    depth = atoi(argv[2]);

  status = fscanf(file, "%d", &num_rows);
  if( status != 1 )
  {
    printf("Could not read number of rows.\n");
    return 1;
  }

  status = fscanf(file, "%d", &num_cols_x);
  if( status != 1 )
  {
    printf("Could not read number of covariates.\n");
    return 1;
  }

  status = fscanf(file, "%d", &num_cols_y);
  if( status != 1 )
  {
    printf("Could not read number of actions.\n");
    return 1;
  }

  data_x = (double *) malloc(num_rows*num_cols_x * sizeof(double));
  data_y = (double *) malloc(num_rows*num_cols_y * sizeof(double));

  for(row = 0; row < num_rows; row++)
  {
    for(col = 0; col < num_cols_x; col++)
    {
      status = fscanf(file, "%lg", data_x+col*num_rows+row);
      if( status != 1)
        break;
    }
    for(col = 0; col < num_cols_y; col++)
    {
      status = fscanf(file, "%lg", data_y+col*num_rows+row);
      if( status != 1)
        break;
    }
  }

  fclose(file);
  
  if( status == 1)
  {
    tree_search(depth, split_step, min_node_size, data_x, data_y, num_rows, num_cols_x, num_cols_y);
  }
  else
  {
    printf("Error reading row=%d, col=%d or col=%d.\n", row, col, col+num_cols_x);
  }

  free(data_x);
  free(data_y);

  return 1-status;
}
