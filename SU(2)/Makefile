IDIR =.
CC=gcc
CFLAGS=-O3  -std=c89 -Wall -Werror -pedantic -fstrict-aliasing -I$(IDIR)

ODIR=.

LIBS=-lm

_DEPS = su2.h ranlxd.h lattice.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

#_OBJ = hopping.o  ranlxd.o phi4.o measure.o gauss_rand.o hmc.o
_OBJ = hopping.o  ranlxd.o gauss_rand.o su2.o heatbath.o 
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))


$(ODIR)/%.o: %.c $(DEPS)
		$(CC) -c -o $@ $< $(CFLAGS)

su2: $(OBJ)
		gcc -o $@ $^ $(CFLAGS) $(LIBS)


.PHONY: clean

clean:
		rm -f $(ODIR)/*.o  core 


