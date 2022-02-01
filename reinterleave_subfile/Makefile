# standard optimisation
#CFLAGS = -O
# vectorize loops
CFLAGS=-O -D_FILE_OFFSET_BITS=64
CFLAGS_VEC=-ftree-vectorize -ftree-vectorizer-verbose=2
# debugging
#CFLAGS = -g
# profiling
#CFLAGS = -pg

# set default install directory
prefix = $(HOME)
INSTALL=install

TARGETS=reinterleave_subfile

all: $(TARGETS)

install: $(TARGETS)
	mkdir -p $(prefix)/bin
	$(INSTALL) $(TARGETS) $(prefix)/bin

uninstall:
	$(foreach target,$(TARGETS),rm -f $(prefix)/bin/$(target);)

reinterleave_subfile: reinterleave_subfile.c
	$(CC) $(CFLAGS) -Wall -D_FILE_OFFSET_BITS=64 $^ -o $@

clean:
	rm -f reinterleave_subfile
