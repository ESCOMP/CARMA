# Yippee!  This is my very own makefile.
# Jamison A. Smith, AKA JAS
# April 26, 2007

FORTRAN =	ifort
#FORTRAN =	pgf90
#FORTRAN =	pathf90
#FORTRAN =	gfortran
#FORTRAN =	g95
#FORTRAN =	xlf90

F90DOC = ../../bin/f90doc-0.4.0/f90doc

PACKAGE =	CARMA
TGZ =		CARMA.tar

FFLAGS =
#FFLAGS += -DSINGLE                    # for single precision
#FFLAGS += -DDEBUG                     # for debug print statements


# Add options for the Intel Fortran compiler.
ifeq ($(FORTRAN),ifort)
#  FFLAGS += -ftz -fp-model precise
  FFLAGS += -fp-model precise
  
  # Work around for an incompatibility with some versions of ifort and OSX.
#  FFLAGS += -use-asm

  # Debug options.
  FFLAGS += -g -O0 -traceback -fp-stack-check -check bounds -check uninit -fpe0 -ftrapuv
  
  # Open/MP
  FFLAGS += -openmp
  
  # The no_pie flags also the executable to work with idb.
  LDFLAGS = $(FFLAGS) -no_pie
endif

# Add options for the Portland Group compiler.
ifeq ($(FORTRAN),pgf90)
  FFLAGS  += 

  # Debug options.
#  FFLAGS += -g -O0 -Mbounds

  # Open/MP
#  FFLAGS  += -mp

  LDFLAGS = $(FFLAGS)
endif

# Add options for the g95 compiler.
ifeq ($(FORTRAN),g95)
#  FFLAGS  += -fzero -ffree-line-length-huge
  FFLAGS  += -ffree-line-length-huge

  # Debug options.
#  FFLAGS += -g -fbounds-check -ftrace=full
  
  # Open/MP
  #
  # NOTE: g95 does not support Open/MP directives. This will cause one
  # test (carma_test) to fail to link.

  LDFLAGS = $(FFLAGS)
endif

# Add options for the IBM XL Fortran compiler.
#
# NOTE: It doesn't support float to zero.
ifeq ($(FORTRAN),xlf90)
  FFLAGS  += -q64 -qarch=auto -qspillsize=2500 -g -qfullpath 
 
  # Debug options.
  FFLAGS += -qinitauto=7FF7FFFF -qflttrap=ov:zero:inv:en -C
  
  # Open/MP
#  FFLAGS += -qsmp=omp
#  FFLAGS += -qsmp=omp:noopt

  LDFLAGS = $(FFLAGS)
endif


# Overridning the implicit rules, which would try to use m2c to
# create the .mod.
%.mod : %.o ;
%.o : %.F90 ;
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
CARMA.exe : $(CARMA_OBJ) carma_test.o carma_testutils.o atmosphere_mod.o
	$(FORTRAN) $(LDFLAGS) -o CARMA.exe carma_test.o carma_testutils.o atmosphere_mod.o $(CARMA_OBJ)
FALLTEST.exe : $(CARMA_OBJ) carma_falltest.o carma_testutils.o atmosphere_mod.o
	$(FORTRAN) $(LDFLAGS) -no_pie -o FALLTEST.exe carma_falltest.o carma_testutils.o atmosphere_mod.o $(CARMA_OBJ)
SIGMAFALLTEST.exe : $(CARMA_OBJ) carma_sigmafalltest.o carma_testutils.o atmosphere_mod.o
	$(FORTRAN) $(LDFLAGS) -o SIGMAFALLTEST.exe carma_sigmafalltest.o carma_testutils.o atmosphere_mod.o $(CARMA_OBJ)
COAGTEST.exe : $(CARMA_OBJ) carma_coagtest.o carma_testutils.o atmosphere_mod.o
	$(FORTRAN) $(LDFLAGS) -o COAGTEST.exe carma_coagtest.o carma_testutils.o atmosphere_mod.o $(CARMA_OBJ)
BCOCTEST.exe : $(CARMA_OBJ) carma_bcoctest.o carma_testutils.o atmosphere_mod.o
	$(FORTRAN) $(LDFLAGS) -o BCOCTEST.exe carma_bcoctest.o carma_testutils.o atmosphere_mod.o $(CARMA_OBJ)
BC2GTEST.exe : $(CARMA_OBJ) carma_bc2gtest.o carma_testutils.o atmosphere_mod.o
	$(FORTRAN) $(LDFLAGS) -o BC2GTEST.exe carma_bc2gtest.o carma_testutils.o atmosphere_mod.o $(CARMA_OBJ)
GROWTEST.exe : $(CARMA_OBJ) carma_growtest.o carma_testutils.o atmosphere_mod.o
	$(FORTRAN) $(LDFLAGS) -o GROWTEST.exe carma_growtest.o carma_testutils.o atmosphere_mod.o $(CARMA_OBJ)
GROWCLRTEST.exe : $(CARMA_OBJ) carma_growclrtest.o carma_testutils.o atmosphere_mod.o
	$(FORTRAN) $(LDFLAGS) -o GROWCLRTEST.exe carma_growclrtest.o carma_testutils.o atmosphere_mod.o $(CARMA_OBJ)
GROWINTEST.exe : $(CARMA_OBJ) carma_growintest.o carma_testutils.o atmosphere_mod.o
	$(FORTRAN) $(LDFLAGS) -o GROWINTEST.exe carma_growintest.o carma_testutils.o atmosphere_mod.o $(CARMA_OBJ)
GROWSUBTEST.exe : $(CARMA_OBJ) carma_growsubtest.o carma_testutils.o atmosphere_mod.o
	$(FORTRAN) $(LDFLAGS) -o GROWSUBTEST.exe carma_growsubtest.o carma_testutils.o atmosphere_mod.o $(CARMA_OBJ)
INITTEST.exe : $(CARMA_OBJ) carma_inittest.o carma_testutils.o atmosphere_mod.o
	$(FORTRAN) $(LDFLAGS) -o INITTEST.exe carma_inittest.o carma_testutils.o atmosphere_mod.o $(CARMA_OBJ)
MIETEST.exe : $(CARMA_OBJ) carma_mietest.o carma_testutils.o 
	$(FORTRAN) $(LDFLAGS) -o MIETEST.exe carma_mietest.o carma_testutils.o  $(CARMA_OBJ)
NUCTEST.exe : $(CARMA_OBJ) carma_nuctest.o carma_testutils.o 
	$(FORTRAN) $(LDFLAGS) -o NUCTEST.exe carma_nuctest.o carma_testutils.o atmosphere_mod.o $(CARMA_OBJ)
NUC2TEST.exe : $(CARMA_OBJ) carma_nuc2test.o carma_testutils.o 
	$(FORTRAN) $(LDFLAGS) -o NUC2TEST.exe carma_nuc2test.o carma_testutils.o atmosphere_mod.o $(CARMA_OBJ)
PHEATTEST.exe : $(CARMA_OBJ) carma_pheattest.o carma_testutils.o atmosphere_mod.o
	$(FORTRAN) $(LDFLAGS) -o PHEATTEST.exe carma_pheattest.o carma_testutils.o atmosphere_mod.o $(CARMA_OBJ)
SWELLTEST.exe : $(CARMA_OBJ) carma_swelltest.o carma_testutils.o 
	$(FORTRAN) $(LDFLAGS) -o SWELLTEST.exe carma_swelltest.o carma_testutils.o atmosphere_mod.o $(CARMA_OBJ)
VDIFTEST.exe : $(CARMA_OBJ) carma_vdiftest.o carma_testutils.o 
	$(FORTRAN) $(LDFLAGS) -o VDIFTEST.exe carma_vdiftest.o carma_testutils.o atmosphere_mod.o $(CARMA_OBJ)
DRYDEPTEST.exe : $(CARMA_OBJ) carma_drydeptest.o carma_testutils.o 
	$(FORTRAN) $(LDFLAGS) -o DRYDEPTEST.exe carma_drydeptest.o carma_testutils.o atmosphere_mod.o $(CARMA_OBJ)	
SIGMADRYDEPTEST.exe : $(CARMA_OBJ) carma_sigmadrydeptest.o carma_testutils.o atmosphere_mod.o
	$(FORTRAN) $(LDFLAGS) -o SIGMADRYDEPTEST.exe carma_sigmadrydeptest.o carma_testutils.o atmosphere_mod.o $(CARMA_OBJ)
SCFALLTEST.exe : $(CARMA_OBJ) carma_scfalltest.o carma_testutils.o atmosphere_mod.o
	$(FORTRAN) $(LDFLAGS) -o SCFALLTEST.exe carma_scfalltest.o carma_testutils.o atmosphere_mod.o $(CARMA_OBJ)
SULFATETEST.exe : $(CARMA_OBJ) carma_sulfatetest.o carma_testutils.o atmosphere_mod.o
	$(FORTRAN) $(LDFLAGS) -o SULFATETEST.exe carma_sulfatetest.o carma_testutils.o atmosphere_mod.o $(CARMA_OBJ)

# Compile everything.
all : FALLTEST.exe COAGTEST.exe BCOCTEST.exe BC2GTEST.exe GROWTEST.exe INITTEST.exe \
MIETEST.exe NUCTEST.exe SIGMAFALLTEST.exe SWELLTEST.exe VDIFTEST.exe DRYDEPTEST.exe \
SIGMADRYDEPTEST.exe PHEATTEST.exe SCFALLTEST.exe CARMA.exe GROWSUBTEST.exe \
SULFATETEST.exe NUC2TEST.exe GROWINTEST.exe GROWCLRTEST.exe

# Compile all of the documentation.
doc : $(CARMA_DOC) $(TEST_DOC)

clean:
	/bin/rm -f *.o *.mod *.exe *.txt *.html

# The Mac creates .DS_Store files that we don't want in the tar file, so
# exclude them.
tar:
	tar --directory ../.. -cvf $(TGZ) --exclude .DS_Store --exclude .svn \
	  Makefile make-carma.csh run-carma.csh README run-regress.csh view-bench.csh run-all.csh \
	  source tests bin doc/ChangeLog doc/ChangeLog_template doc/index.html
