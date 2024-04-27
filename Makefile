SRC = main.c  reading.c  simple_opttree.c  tree.c  workspace.c simple_set.c sorted_set.c strategy.c units.c \
	reading.h  simple_opttree.h  tree.h  workspace.h simple_set.h sorted_set.h strategy.h units.h type_all.h versiongit.h
SRC2 = simple_opttree.c  tree.c  workspace.c type_all.h
SRCDIR = src
MAINSRC = $(addprefix $(SRCDIR)/,$(SRC))
MAINSRC2 = $(addprefix $(SRCDIR)/,$(SRC2))
TESTS = $(SRCDIR)/tests.c

fpt : $(MAINSRC) 
	gcc -g -DPRINTING_ALLOWED -DNDEBUG -Wall -O2 -o fpt $(MAINSRC) 

fpt_dbg : $(MAINSRC) 
	gcc -g -DPRINTING_ALLOWED -Wall -o fpt_dbg $(MAINSRC) 

tests: $(MAINSRC2) $(SIMPLE_SET) $(TESTS)
	gcc -g -DUSE_SIMPLE_SET -DPRINTING_ALLOWED -Wall -o tests $(MAINSRC2) $(SIMPLE_SET) $(TESTS)
