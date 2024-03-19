SRC = main.c  reading.c  simple_opttree.c  tree.c  workspace.c

sorted_set : $(SRC) sorted_set.c
	gcc -g -DUSE_SORTED_SET -DNDEBUG -Wall -O2 -o sorted_set $(SRC) sorted_set.c

simple_set : $(SRC) simple_set.c
	gcc -g -DUSE_SIMPLE_SET -DNDEBUG -Wall -O2 -o simple_set $(SRC) simple_set.c

sorted_set_dbg : $(SRC) sorted_set.c
	gcc -g -DUSE_SORTED_SET -Wall -o sorted_set_dbg $(SRC) sorted_set.c

simple_set_dbg : $(SRC) simple_set.c
	gcc -g -DUSE_SIMPLE_SET -Wall -o simple_set_dbg $(SRC) simple_set.c
