TARGET=test_su2

OBJ= test_su2.o lattice_setup.o inits.o rnd.o monte.o

FORTRAN_FLAG   = 
LIB =

all : $(OBJ)
	f77 $(FORTRAN_FLAG) -o $(TARGET).out $(OBJ) $(LIB)

$(TARGET).o : $(TARGET).f parameters.f 
	f77 $(FORTRAN_FLAG) -c $(TARGET).f

lattice_setup.o : lattice_setup.f parameters.f 
	f77 $(FORTRAN_FLAG) -c lattice_setup.f

inits.o : inits.f parameters.f
	f77 $(FORTRAN_FLAG) -c inits.f

rnd.o : rnd.f parameters.f
	f77 $(FORTRAN_FLAG) -c rnd.f

monte.o : monte.f parameters.f
	f77 $(FORTRAN_FLAG) -c monte.f

clean_obj : 
	rm $(OBJ)

.PHONY: clean
clean : 
	-rm *.out *.o

