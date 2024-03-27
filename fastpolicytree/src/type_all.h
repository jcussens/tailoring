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
#define USE_SIMPLE_SET 
/* #define PRINTING_ALLOWED */

#include "versiongit.h"

struct workspace;
typedef struct workspace WORKSPACE;  /**< Work space structure */

#ifdef USE_SORTED_SET
struct sorted_set;
typedef struct sorted_set SORTED_SET;  /**< Sorted sets */
typedef SORTED_SET** UNITS;
typedef const SORTED_SET** CONST_UNITS; /**< Constant sorted sets */
#endif

#ifdef USE_SIMPLE_SET
struct simple_set;
typedef struct simple_set SIMPLE_SET; /**< Simple sets */
typedef SIMPLE_SET* UNITS;
typedef const SIMPLE_SET* CONST_UNITS; /**< Constant simple sets */
#endif

struct node;
typedef struct node NODE;  /**< Policy tree structure */

typedef unsigned int ELEMENT;

#ifdef __cplusplus
}
#endif

#endif
