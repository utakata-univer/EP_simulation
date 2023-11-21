IDIR =.
CC=gcc
CFLAGS=-O3  -std=c89 -Wall -Werror -pedantic -fstrict-aliasing -I$(IDIR)

ODIR=.

LIBS=-lm

_DEPS = ranlxd.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = ranlxd.o gauss_rand.c
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))


$(ODIR)/%.o: %.c $(DEPS)
		$(CC) -c -o $@ $< $(CFLAGS)

gauss_rand: $(OBJ)
		gcc -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean

clean:
		rm -f $(ODIR)/*.o  core 


