F90 = ifort#gfortran
FFLAGS = -check all -r8#-fcheck=all -ffree-line-length-none
%.o: %.f90
	$(F90) $< $(FFLAGS) -c  -o $@ 

shellmodel: constants.o  neighbours_lists.o shells.o
	$(F90) -o shellmodel constants.o  neighbours_lists.o shells.o
clean:
	@rm -f *.o
	@rm -f *.mod
	@rm -f shellmodel
