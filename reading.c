#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
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

  assert(line!=NULL);
  
  for( i = 0; i < strlen(line); i++ )
  {
    if( !isspace(line[i]) && (line[i+1] == '\0' || isspace(line[i+1]) ) )
      nfields++;
  }
  return nfields;
}

int readfile(
   char*                 filename,           /**< name of file with data */
   int                   num_cols_y,         /**< num_cols_y is the number of actions */
   double**              data_x,             /**< *data_x will be covariate data (column major) */
   double**              data_y,             /**< *data_y will be reward data (column major) */
   int*                  num_rows,           /**< *num_rows will be number of rows (=individuals) in data */
   int*                  num_cols_x,         /**< *num_cols_x will be number of covariates */
   char***               covnames,           /**< (*covnames)[i] will be the name of the ith covariate */
   char***               actionnames         /**< (*actionnames)[i] will be the name of the ith action */
   )
{
   FILE* file;
   char line[MAX_LINE_LENGTH];
   int nlines;
   int i;
   int nfields;
   char* lineptr;
   int status;

   int row;
   int col;

   assert(filename != NULL);
   assert(num_cols_y >= 0);
   assert(data_x != NULL);
   assert(data_y != NULL);
   assert(num_rows != NULL);
   assert(num_cols_x != NULL);
   assert(covnames != NULL);
   assert(actionnames != NULL);

   file = fopen(filename,"r");

   if( file == NULL )
   {
      printf("Could not open %s.\n", filename);
      return 1;
   }    
  
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
   *num_rows = nlines - 1;
  
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
   
   /* compute number of covariates */

   *num_cols_x = nfields - num_cols_y;

  /* read in covariate and action names */
  
   *covnames = (char**) malloc(*num_cols_x*sizeof(char*));
   *actionnames = (char**) malloc(num_cols_y*sizeof(char*));
   lineptr = line;
   for( i = 0; i < *num_cols_x+num_cols_y; i++ )
   {
      int offset;
      char tmpstr[MAX_COVARIATE_NAME_LENGTH];

      status = sscanf(lineptr, "%s%n", tmpstr, &offset);
      if( status != 1 )
      {
         if( i < *num_cols_x )
            printf("Error reading covariate name with index %d from %s.\n", i, line);
         else
            printf("Error reading action name with index %d from %s.\n", i-*num_cols_x, line);
         return 1;
      }

      /* add an 'X' if name starts with a digit (to be like policytree) */
      if( isdigit(*tmpstr) )
      {
         memmove(tmpstr+1,tmpstr,(strlen(tmpstr)+1));
         tmpstr[0] = 'X';
      }
    
      if( i < *num_cols_x )
      {
         (*covnames)[i] = (char*) malloc((strlen(tmpstr)+1)*sizeof(char));
         strcpy((*covnames)[i],tmpstr);
      }
      else
      {
         (*actionnames)[i-*num_cols_x] = (char*) malloc((strlen(tmpstr)+1)*sizeof(char));
         strcpy((*actionnames)[i-*num_cols_x],tmpstr);
      }
      lineptr += offset;
   }

   *data_x = (double *) malloc((*num_rows) * (*num_cols_x) * sizeof(double));
   *data_y = (double *) malloc((*num_rows) * num_cols_y * sizeof(double));

   for(row = 0; row < *num_rows; row++)
   {
      for(col = 0; col < *num_cols_x; col++)
      {
         status = fscanf(file, "%lg", (*data_x) + col * (*num_rows) + row);
         if( status != 1)
         {
            printf("Error reading in value for covariate %s with index %d for datapoint %d.\n", (*covnames)[col], col, row);
            return 1;
         }
      }
      for(col = 0; col < num_cols_y; col++)
      {
         status = fscanf(file, "%lg", (*data_y) + col * (*num_rows) + row);
         if( status != 1)
         {
            printf("Error reading in value for action %d for datapoint %d.\n", col, row);
            return 1;
         }
      }
   }

   fclose(file);

   return 0;
   
}
