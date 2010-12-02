# Sub makefile for test files

# This is intended to be included by another makefile to actually
# build the system. It has all the dependency information for the
# files in the test tree.

# NOTE: In the future (or in parnet models) this could be replaced
# by automatic dependency generation.

TEST_OBJ = carma_drydeptest.o carma_falltest.o carma_coagtest.o carma_test.o atmosphere_mod.o

TEST_DOC = carma_drydeptest.html carma_falltest.html carma_coagtest.html carma_test.html atmosphere_mod.html

atmosphere_mod.o : atmosphere_mod.F90
	$(FORTRAN) $(CPPFLAGS) -c $<

carma_coagtest.o : carma_coagtest.F90 atmosphere_mod.mod carma_mod.mod carma_constants_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(CPPFLAGS) -c $<
	
carma_bcoctest.o : carma_bcoctest.F90 atmosphere_mod.mod carma_mod.mod carma_constants_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(CPPFLAGS) -c $<
	
carma_bc2gtest.o : carma_bc2gtest.F90 atmosphere_mod.mod carma_mod.mod carma_constants_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(CPPFLAGS) -c $<
	
carma_falltest.o : carma_falltest.F90 atmosphere_mod.mod carma_mod.mod carma_constants_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(CPPFLAGS) -c $<
	
carma_growtest.o : carma_growtest.F90 atmosphere_mod.mod carma_mod.mod carma_constants_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(CPPFLAGS) -c $<
	
carma_inittest.o : carma_inittest.F90 atmosphere_mod.mod carma_mod.mod carma_constants_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(CPPFLAGS) -c $<

carma_mietest.o : carma_mietest.F90 carma_mod.mod carma_constants_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(CPPFLAGS) -c $<
	
carma_nuctest.o : carma_nuctest.F90 atmosphere_mod.mod carma_mod.mod carma_constants_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(CPPFLAGS) -c $<
	
carma_sigmadrydeptest.o : carma_sigmadrydeptest.F90 atmosphere_mod.mod carma_mod.mod carma_constants_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(CPPFLAGS) -c $<
	
carma_sigmafalltest.o : carma_sigmafalltest.F90 atmosphere_mod.mod carma_mod.mod carma_constants_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(CPPFLAGS) -c $<
	
carma_swelltest.o : carma_swelltest.F90 atmosphere_mod.mod carma_mod.mod carma_constants_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(CPPFLAGS) -c $<
	
carma_test.o : carma_test.F90 atmosphere_mod.mod carma_mod.mod carma_constants_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(CPPFLAGS) -c $<
	
carma_vdiftest.o : carma_vdiftest.F90 atmosphere_mod.mod carma_mod.mod carma_constants_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(CPPFLAGS) -c $<

