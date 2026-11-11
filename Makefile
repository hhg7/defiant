# Makefile for Defiant
#
#	make              optimized build            -> ./defiant
#	make debug        -O0 -g3 build              -> ./defiant_debug  (gdb, valgrind)
#	make asan         ASan + UBSan build         -> ./defiant_asan
#	make check        smoke test of ./defiant
#	make valgrind     smoke test of ./defiant_debug under valgrind
#	make install      to $(DESTDIR)$(PREFIX)/bin  (PREFIX defaults to /usr/local)
#	make uninstall
#	make clean
#
# Common overrides:
#	make CC=clang
#	make OPT=-O2
#	make OPENMP=1                  # adds -fopenmp
#	make install PREFIX=$$HOME/.local

PROG    := defiant
SRC     := defiant.c
# every header sits next to defiant.c and is compiled into the one translation
# unit, so any header change must trigger a rebuild
HEADERS := $(wildcard *.h)

PREFIX  ?= /usr/local
BINDIR  ?= $(PREFIX)/bin
INSTALL ?= install

# gnu11, not c11: getline, strtok_r, strdup, strcasecmp and ssize_t are needed
CSTD    := -std=gnu11
WARN    := -Wall -Wextra -pedantic -Wconversion
OPT     ?= -O3
CFLAGS  ?= $(OPT)
# -lm has to come after the source file on the link line, hence LDLIBS
LDLIBS  := -lm

ifeq ($(OPENMP),1)
OMPFLAG := -fopenmp
endif

DEBUG_CFLAGS := -O0 -g3
# -Wno-format-overflow: with sanitizers on, gcc reports "null destination pointer" for
# sprintf calls whose buffer comes from a checked malloc; these are false positives
ASAN_CFLAGS  := -O1 -g3 -fno-omit-frame-pointer -fsanitize=address,undefined -Wno-format-overflow

.PHONY: all debug asan check valgrind install uninstall clean
.DELETE_ON_ERROR:

all: $(PROG)

$(PROG): $(SRC) $(HEADERS)
	$(CC) $(CSTD) $(WARN) $(CPPFLAGS) $(CFLAGS) $(OMPFLAG) $(LDFLAGS) -o $@ $(SRC) $(LDLIBS)

debug: $(PROG)_debug

$(PROG)_debug: $(SRC) $(HEADERS)
	$(CC) $(CSTD) $(WARN) $(CPPFLAGS) $(DEBUG_CFLAGS) $(OMPFLAG) $(LDFLAGS) -o $@ $(SRC) $(LDLIBS)

asan: $(PROG)_asan

$(PROG)_asan: $(SRC) $(HEADERS)
	$(CC) $(CSTD) $(WARN) $(CPPFLAGS) $(ASAN_CFLAGS) $(OMPFLAG) $(LDFLAGS) -o $@ $(SRC) $(LDLIBS)

check: $(PROG)
	./$(PROG) -version
	./$(PROG) -h > /dev/null
	@echo "smoke test passed"

valgrind: $(PROG)_debug
	valgrind --error-exitcode=1 --leak-check=full --track-origins=yes ./$(PROG)_debug -version

install: $(PROG)
	$(INSTALL) -d $(DESTDIR)$(BINDIR)
	$(INSTALL) -m 0755 $(PROG) $(DESTDIR)$(BINDIR)/$(PROG)

uninstall:
	$(RM) $(DESTDIR)$(BINDIR)/$(PROG)

clean:
	$(RM) $(PROG) $(PROG)_debug $(PROG)_asan
