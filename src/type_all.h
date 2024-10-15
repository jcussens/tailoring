/** @file type_all.h
 *  @brief Type definitions
 *  @author James Cussens
 */
#ifndef __TYPE_ALL_H__
#define __TYPE_ALL_H__

#ifdef __cplusplus
extern "C" {
#endif

#if defined __has_include
#  if __has_include ("versiongit.h")
#    include "versiongit.h"
#  endif
#endif


struct workspace;
typedef struct workspace WORKSPACE;  /**< Work space structure */

typedef void* UNITS;                 /**< Datasets independent of how represented */
typedef const void* CONST_UNITS;     /**< Constant datasets independent of how represented */

struct node;
typedef struct node NODE;  /**< Policy tree structure */

struct strategy;
typedef struct strategy STRATEGY;  /**< Records solving strategy */

typedef unsigned int ELEMENT;      /**< A member of a set of UNITS */

#ifdef __cplusplus
}
#endif

#endif
