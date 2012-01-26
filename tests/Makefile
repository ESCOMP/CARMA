# Sub makefile for test files

# This is intended to be included by another makefile to actually
# build the system. It has all the dependency information for the
# files in the test tree.

# NOTE: In the future (or in parnet models) this could be replaced
# by automatic dependency generation.

TEST_OBJ = carma_bcoctest.o carma_bc2gtest.o carma_coagtest.o carma_drydeptest.o carma_falltest.o \
carma_growtest.o carma_inittest.o carma_mietest.o carma_nuctest.o \
carma_pheattest.o carma_sigmadrydeptest.o carma_sigmafalltest.o carma_swelltest.o \
carma_vdiftest.o carma_test.o atmosphere_mod.o carma_scfalltest.o carma_growsubtest.o \
carma_growintest.o

TEST_DOC = carma_bcoctest.html carma_bc2gtest.html carma_coagtest.html carma_drydeptest.html \
carma_falltest.html carma_growtest.html carma_inittest.html carma_mietest.html carma_nuctest.html \
carma_pheattest.html carma_sigmadrydeptest.html carma_sigmafalltest.html carma_swelltest.html \
carma_vdiftest.html carma_test.html atmosphere_mod.html carma_scfalltest.html carma_growsubtest.html \
carma_growintest.html

atmosphere_mod.o : atmosphere_mod.F90
	$(FORTRAN) $(FFLAGS) -c $<

carma_testutils.o : carma_testutils.F90
	$(FORTRAN) $(FFLAGS) -c $<

carma_bcoctest.o : carma_bcoctest.F90 atmosphere_mod.mod carma_mod.mod carma_constants_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(FFLAGS) -c $<
	
carma_bc2gtest.o : carma_bc2gtest.F90 atmosphere_mod.mod carma_mod.mod carma_constants_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(FFLAGS) -c $<
	
carma_coagtest.o : carma_coagtest.F90 atmosphere_mod.mod carma_mod.mod carma_constants_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(FFLAGS) -c $<
	
carma_drydeptest.o : carma_drydeptest.F90 atmosphere_mod.mod carma_mod.mod carma_constants_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(FFLAGS) -c $<
	
carma_falltest.o : carma_falltest.F90 atmosphere_mod.mod carma_mod.mod carma_constants_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(FFLAGS) -c $<
	
carma_growtest.o : carma_growtest.F90 atmosphere_mod.mod carma_mod.mod carma_constants_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(FFLAGS) -c $<
	
carma_growclrtest.o : carma_growclrtest.F90 atmosphere_mod.mod carma_mod.mod carma_constants_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(FFLAGS) -c $<

carma_growintest.o : carma_growintest.F90 atmosphere_mod.mod carma_mod.mod carma_constants_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(FFLAGS) -c $<

carma_growsubtest.o : carma_growsubtest.F90 atmosphere_mod.mod carma_mod.mod carma_constants_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(FFLAGS) -c $<
	
carma_inittest.o : carma_inittest.F90 atmosphere_mod.mod carma_mod.mod carma_constants_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(FFLAGS) -c $<

carma_mietest.o : carma_mietest.F90 carma_mod.mod carma_constants_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(FFLAGS) -c $<
	
carma_nuctest.o : carma_nuctest.F90 atmosphere_mod.mod carma_mod.mod carma_constants_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(FFLAGS) -c $<
	
carma_nuc2test.o : carma_nuc2test.F90 atmosphere_mod.mod carma_mod.mod carma_constants_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(FFLAGS) -c $<
	
carma_pheattest.o : carma_pheattest.F90 atmosphere_mod.mod carma_mod.mod carma_constants_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(FFLAGS) -c $<
	
carma_scfalltest.o : carma_scfalltest.F90 atmosphere_mod.mod carma_mod.mod carma_constants_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(FFLAGS) -c $<
	
carma_sigmadrydeptest.o : carma_sigmadrydeptest.F90 atmosphere_mod.mod carma_mod.mod carma_constants_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(FFLAGS) -c $<
	
carma_sigmafalltest.o : carma_sigmafalltest.F90 atmosphere_mod.mod carma_mod.mod carma_constants_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(FFLAGS) -c $<
	
carma_swelltest.o : carma_swelltest.F90 atmosphere_mod.mod carma_mod.mod carma_constants_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(FFLAGS) -c $<
	
carma_sulfatetest.o : carma_sulfatetest.F90 atmosphere_mod.mod carma_mod.mod carma_constants_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(FFLAGS) -c $<
	
carma_test.o : carma_test.F90 atmosphere_mod.mod carma_mod.mod carma_constants_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(FFLAGS) -c $<
	
carma_vdiftest.o : carma_vdiftest.F90 atmosphere_mod.mod carma_mod.mod carma_constants_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(FFLAGS) -c $<

