CC = gcc
SRCDIR = src
FLAGS = -g -DPRINTING_ALLOWED -Wall -O2
ifneq ($(OPT),dbg)
	FLAGS += -DNDEBUG
endif
SOURCES = $(wildcard $(SRCDIR)/*.c)
OBJECTS = $(SOURCES:%.c=%.o)

fpt: $(OBJECTS)
	$(CC) $(FLAGS) -o $@ $(OBJECTS)

include $(SOURCES:.c=.d)

# From GNU Makefile manual, but with -MT added
$(SRCDIR)/%.d: $(SRCDIR)/%.c
	@set -e; rm -f $@; \
	 $(CC) -MM $< -MT $(<:.c=.o) > $@.$$$$; \
	 sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	 rm -f $@.$$$$

$(SRCDIR)/%.o:	$(SRCDIR)/%.c
	$(CC) $(FLAGS) -c $< -o $@

.PHONY : clean
clean :
	rm fpt $(OBJECTS)

.PHONY : doc
doc :
	doxygen

