SRC = main.c  reading.c  simple_opttree.c  tree.c  workspace.c type_all.h versiongit.h
SRC2 = simple_opttree.c  tree.c  workspace.c type_all.h
SRCDIR = src
MAINSRC = $(addprefix $(SRCDIR)/,$(SRC))
MAINSRC2 = $(addprefix $(SRCDIR)/,$(SRC2))
SORTED_SET = $(SRCDIR)/sorted_set.c
SIMPLE_SET = $(SRCDIR)/simple_set.c
TESTS = $(SRCDIR)/tests.c

sorted_set : $(MAINSRC) $(SORTED_SET)
	gcc -g -DUSE_SORTED_SET -DPRINTING_ALLOWED -DNDEBUG -Wall -O2 -o sorted_set $(MAINSRC) $(SORTED_SET)

simple_set : $(MAINSRC) $(SIMPLE_SET)
	gcc -g -DUSE_SIMPLE_SET -DPRINTING_ALLOWED -DNDEBUG -Wall -O2 -o simple_set $(MAINSRC) $(SIMPLE_SET)

sorted_set_dbg : $(MAINSRC) $(SORTED_SET)
	gcc -g -DUSE_SORTED_SET -DPRINTING_ALLOWED -Wall -o sorted_set_dbg $(MAINSRC) $(SORTED_SET)

simple_set_dbg : $(MAINSRC) $(SIMPLE_SET)
	gcc -g -DUSE_SIMPLE_SET -DPRINTING_ALLOWED -Wall -o simple_set_dbg $(MAINSRC) $(SIMPLE_SET)

tests: $(MAINSRC2) $(SIMPLE_SET) $(TESTS)
	gcc -g -DUSE_SIMPLE_SET -DPRINTING_ALLOWED -Wall -o tests $(MAINSRC2) $(SIMPLE_SET) $(TESTS)
