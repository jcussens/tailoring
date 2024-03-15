#ifndef __TYPE_ALL_H__
#define __TYPE_ALL_H__

#ifdef __cplusplus
extern "C" {
#endif

/* #define USE_SORTED_SET */
#define USE_SIMPLE_SET

struct workspace;
typedef struct workspace WORKSPACE;  /**< Work space structure */

#ifdef USE_SORTED_SET
struct sorted_set;
typedef struct sorted_set SORTED_SET;  
typedef SORTED_SET** UNITS;
typedef const SORTED_SET** CONST_UNITS;
#endif

#ifdef USE_SIMPLE_SET
struct simple_set;
typedef struct simple_set SIMPLE_SET;
typedef SIMPLE_SET* UNITS;
typedef const SIMPLE_SET* CONST_UNITS;
#endif

struct node;
typedef struct node NODE;  /**< Policy tree structure */

#ifdef __cplusplus
}
#endif

#endif
