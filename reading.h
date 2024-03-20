/** @file reading.h
 *  @brief Functions for reading in data
 *  @author James Cussens
 */
#ifndef __READING_H__
#define __READING_H__

#ifdef __cplusplus
extern "C" {
#endif

/** read data from a given file name and number of actions
 * @return 0 if all is well, else 1
 */
int readfile(
   char*                 filename,           /**< name of file with data */
   int                   num_cols_y,         /**< num_cols_y is the number of actions */
   double**              data_x,             /**< *data_x will be covariate data (column major) */
   double**              data_y,             /**< *data_y will be reward data (column major) */
   int*                  num_rows,           /**< *num_rows will be number of rows (=individuals) in data */
   int*                  num_cols_x,         /**< *num_cols_x will be number of covariates */
   char***               covnames,           /**< (*covnames)[i] will be the name of the ith covariate */
   char***               actionnames         /**< (*actionnames)[i] will be the name of the ith action */
   );

#ifdef __cplusplus
}
#endif

#endif
