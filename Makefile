# Yippee!  This is my very own makefile.
# Jamison A. Smith, AKA JAS
# April 26, 2007

FORTRAN =	ifort
#FORTRAN =	pathf90
#FORTRAN =	gfortran
#FORTRAN =	g95

F90DOC = ../../bin/f90doc-0.4.0/f90doc

PACKAGE =	CARMA
TGZ =		CARMA.tgz

CPPFLAGS = -g

#CPPFLAGS += -DSINGLE                    # for single precision
CPPFLAGS += -openmp                     # for Open/MP Threading of CARMA

#CPPFLAGS += -DDEBUG                     # for debugg print statements
CPPFLAGS += -O0                      # for running in a debugger

# Add options for the Intel Fortran compiler on the Mac.
CPPFLAGS += -ftz -traceback -fp-model precise -use-asm

# Debug options for the Intel Fortran compiler on the Mac.
CPPFLAGS += -fp-stack-check -check bounds -check uninit -fpe0 -ftrapuv

# Overridning the implicit rules, which would try to use m2c to
# create the .mod.
%.mod : %.o ;
%.o : %.F90 ; $(FORTRAN) $(CPPFLAGS) -c $<
%.html : %.F90 ; $(F90DOC) -cs $<

# Add the directories where the source files are located.
VPATH := ../../source/base ../../tests

# These makefiles have the object lists and dependence information
# for the respective components.
#
# NOTE: In the future it might be nice to generate this dependency
# try automatically.
include ../../source/base/Makefile
include ../../tests/Makefile

# Rules for each executable that could be build.
CARMA.exe : $(CARMA_OBJ) carma_test.o atmosphere_mod.o
	$(FORTRAN) $(CPPFLAGS) -o CARMA.exe carma_test.o atmosphere_mod.o $(CARMA_OBJ)
FALLTEST.exe : $(CARMA_OBJ) carma_falltest.o atmosphere_mod.o
	$(FORTRAN) $(CPPFLAGS) -o FALLTEST.exe carma_falltest.o atmosphere_mod.o $(CARMA_OBJ)
SIGMAFALLTEST.exe : $(CARMA_OBJ) carma_sigmafalltest.o atmosphere_mod.o
	$(FORTRAN) $(CPPFLAGS) -o SIGMAFALLTEST.exe carma_sigmafalltest.o atmosphere_mod.o $(CARMA_OBJ)
COAGTEST.exe : $(CARMA_OBJ) carma_coagtest.o atmosphere_mod.o
	$(FORTRAN) $(CPPFLAGS) -o COAGTEST.exe carma_coagtest.o atmosphere_mod.o $(CARMA_OBJ)
BCOCTEST.exe : $(CARMA_OBJ) carma_bcoctest.o atmosphere_mod.o
	$(FORTRAN) $(CPPFLAGS) -o BCOCTEST.exe carma_bcoctest.o atmosphere_mod.o $(CARMA_OBJ)
BC2GTEST.exe : $(CARMA_OBJ) carma_bc2gtest.o atmosphere_mod.o
	$(FORTRAN) $(CPPFLAGS) -o BC2GTEST.exe carma_bc2gtest.o atmosphere_mod.o $(CARMA_OBJ)
GROWTEST.exe : $(CARMA_OBJ) carma_growtest.o
	$(FORTRAN) $(CPPFLAGS) -o GROWTEST.exe carma_growtest.o atmosphere_mod.o $(CARMA_OBJ)
INITTEST.exe : $(CARMA_OBJ) carma_inittest.o atmosphere_mod.o
	$(FORTRAN) $(CPPFLAGS) -o INITTEST.exe carma_inittest.o atmosphere_mod.o $(CARMA_OBJ)
MIETEST.exe : $(CARMA_OBJ) carma_mietest.o
	$(FORTRAN) $(CPPFLAGS) -o MIETEST.exe carma_mietest.o $(CARMA_OBJ)
NUCTEST.exe : $(CARMA_OBJ) carma_nuctest.o
	$(FORTRAN) $(CPPFLAGS) -o NUCTEST.exe carma_nuctest.o atmosphere_mod.o $(CARMA_OBJ)
SWELLTEST.exe : $(CARMA_OBJ) carma_swelltest.o
	$(FORTRAN) $(CPPFLAGS) -o SWELLTEST.exe carma_swelltest.o atmosphere_mod.o $(CARMA_OBJ)
VDIFTEST.exe : $(CARMA_OBJ) carma_vdiftest.o
	$(FORTRAN) $(CPPFLAGS) -o VDIFTEST.exe carma_vdiftest.o atmosphere_mod.o $(CARMA_OBJ)

# Compile everything.
all : FALLTEST.exe COAGTEST.exe BCOCTEST.exe BC2GTEST.exe GROWTEST.exe INITTEST.exe MIETEST.exe NUCTEST.exe SIGMAFALLTEST.exe SWELLTEST.exe VDIFTEST.exe CARMA.exe

# Compile all of the documentation.
doc : $(CARMA_DOC) $(TEST_DOC)

clean:
	/bin/rm -f *.o *.mod *.exe *.txt *.html

# The Mac creates .DS_Store files that we don't want in the tar file, so
# exclude them.
tar:
	tar --directory ../.. -cvf $(TGZ) --wildcards --exclude .DS_Store \
	  Makefile make-carma.csh run-carma.csh README \
	  source tests bin doc/ChangeLog doc/ChangeLog_template doc/index.html
