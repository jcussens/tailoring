#ifndef __TYPE_ALL_H__
#define __TYPE_ALL_H__

#ifdef __cplusplus
extern "C" {
#endif

struct workspace;
typedef struct workspace WORKSPACE;  /**< Work space structure */

struct sorted_set;
typedef struct sorted_set SORTED_SET;  /**< Trie structure for storing scored parent sets */

struct node;
typedef struct node NODE;  /**< Policy tree structure */

#ifdef __cplusplus
}
#endif

#endif
