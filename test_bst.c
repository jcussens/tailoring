#include "bst.h"
#include <stdio.h>

void list_elts(
     BST* bst
   )
{
   BST* next;
   int elt;
   
   next = find_minimum(bst,&elt);
   if( next != NULL )
      printf("%d\n",elt);
  
   while( next != NULL )
   {
      next = find_next(next,&elt);
      if( next != NULL )
         printf("%d\n",elt);
   }
   
   printf("\n\n");
}
   
int main(
  int argc,
  char* argv[]
  )
{
  BST* bst;

  int elts[] = {4,6,7,2,3,1,9,5};

  bst = make_bst(elts,8,0,NULL);

  list_elts(bst);
  
  delete(bst,7,2);
  delete(bst,1,5);

  list_elts(bst);
  
  insert(bst,7,2);
  delete(bst,4,0);
  
  list_elts(bst);

  delete_bst(bst);
}
