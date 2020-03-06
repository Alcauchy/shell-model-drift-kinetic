DIR1 = /usr/local/hdf5-1.10.6/include
L_LIBS = /usr/local/hdf5-1.10.6/lib
#LIBS = -lhdf -lhdf_fortran -lh5lib -lh5
INC=$(DIR1) 

F9X = /usr/local/hdf5-1.10.6/bin/h5fc#ifort#gfortran#
FFLAGS = -check all  -r8 -I${INC} -L${L_LIBS}  -module $(OBJDIR) #-fcheck=all -ffree-line-length-none#

#VPATH = 
#vpath %.mod ../include
INC_PARAMS=$(foreach d, $(INC), -I$d)
#LIB = -libhdf5
#LIBS = $(LIB)
DIRS = 
BASEDIR = .
SRCDIR  = $(BASEDIR)/source
OBJDIR  = $(BASEDIR)/obj
BINDIR  = $(BASEDIR)/bin
DIRS += $(BINDIR)
DIRS += $(OBJDIR)

vpath %.f90 $(SRCDIR)


	
$(OBJDIR)/%.o: $(SRCDIR)/%.f90
	$(F9X)   $(FFLAGS) -I$(OBJDIR)/  -c  $<    -o $@
	
$(OBJDIR):
	mkdir $@

$(BINDIR):
	mkdir $@

#$(OBJDIR):
#	mkdir -p $@
	

.PHONY: all
all: shellmodel
shellmodel: 
	@make --no-print-directory $(BINDIR)/shellmodel

$(BINDIR)/shellmodel: $(OBJDIR) $(BINDIR) $(OBJDIR)/constants.o  $(OBJDIR)/neighbours_lists.o $(OBJDIR)/shells.o $(OBJDIR)/IO.o $(OBJDIR)/main.o
	$(F9X) -o $(BINDIR)/shellmodel $(OBJDIR)/*.o
#constants.o  neighbours_lists.o shells.o main.o
.PHONY: clean
clean:
	@rm -f $(OBJDIR)/*.o
	@rm -f $(OBJDIR)/*.mod
	@rm -f $(BINDIR)/*
