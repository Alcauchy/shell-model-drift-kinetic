DIR1 = /usr/local/hdf5-1.10.6/include
L_LIBS = /usr/local/hdf5-1.10.6/lib
#LIBS = -lhdf -lhdf_fortran -lh5lib -lh5
INC=$(DIR1) 

F9X = /usr/local/hdf5-1.10.6/bin/h5fc#ifort#gfortran#
FFLAGS = -check all  -r8 -I${INC} -L${L_LIBS} #-fcheck=all -ffree-line-length-none#

#VPATH = 
#vpath %.mod ../include
INC_PARAMS=$(foreach d, $(INC), -I$d)
#LIB = -libhdf5
#LIBS = $(LIB)

%.o: %.f90
	$(F9X)   $(FFLAGS)   -c  $<    -o $@ ${LIBS}

shellmodel: constants.o  neighbours_lists.o shells.o IO.o main.o
	$(F9X) -o shellmodel *.o
#constants.o  neighbours_lists.o shells.o main.o
clean:
	@rm -f *.o
	@rm -f *.mod
	@rm -f shellmodel
