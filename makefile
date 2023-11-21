IDIR =.
CC=gcc
CFLAGS=-O3  -std=c89 -Wall -Werror -pedantic -fstrict-aliasing -I$(IDIR)

ODIR=.

LIBS=-lm

_DEPS = phi42.h lattice2.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = hopping2.o phi42.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))


$(ODIR)/%.o: %.c $(DEPS)
	        $(CC) -c -o $@ $< $(CFLAGS)

phi42: $(OBJ)
	        gcc -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean

clean:
	        rm -f $(ODIR)/*.o  core
