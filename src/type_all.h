/** @file type_all.h
 *  @brief Type definitions
 *  @author James Cussens
 */
#ifndef __TYPE_ALL_H__
#define __TYPE_ALL_H__

#ifdef __cplusplus
extern "C" {
#endif

/* #define USE_SORTED_SET  */
/* #define USE_SIMPLE_SET */
/* #define PRINTING_ALLOWED */

#include "versiongit.h"

struct workspace;
typedef struct workspace WORKSPACE;  /**< Work space structure */

typedef void* UNITS;
typedef const void* CONST_UNITS; 

struct node;
typedef struct node NODE;  /**< Policy tree structure */

struct strategy;
typedef struct strategy STRATEGY;

typedef unsigned int ELEMENT;

#ifdef __cplusplus
}
#endif

#endif
