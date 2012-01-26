# Sub makefile for CARMA files

# This is intended to be included by another makefile to actually
# build the system. It has all the dependency information for the
# files in the CARMA tree.

# NOTE: In the future (or in parnet models) this could be replaced
# by automatic dependency generation and/or by building CARMA as a
# library.

CARMA_OBJ = carma_constants_mod.o carma_precision_mod.o carma_types_mod.o \
carmagroup_mod.o carmagas_mod.o carmaelement_mod.o carmasolute_mod.o \
carmastate_mod.o carma_mod.o coagl.o coagp.o csolve.o setupbins.o \
setupatm.o setupvf.o microslow.o newstate.o prestep.o setupckern.o \
setupcoag.o smallconc.o step.o supersat.o vaporp.o vertadv.o \
vertdif.o vertical.o versol.o rhopart.o psolve.o zeromicro.o \
nsubsteps.o setupgrow.o setupgkern.o setupnuc.o growevapl.o microfast.o \
gsolve.o actdropl.o freezglaerl_murray2010.o growp.o downgxfer.o \
gasexchange.o melticel.o upgxfer.o freezdropl.o \
downgevapply.o evapp.o evap_poly.o evap_mono.o \
evap_ingrp.o tsolve.o miess.o vaporp_h2o_buck1981.o wetr.o \
vaporp_h2o_murphy2005.o maxconc.o setupbdif.o setupvf_std.o \
setupvf_std_shape.o totalcondensate.o versub.o freezaerl_tabazadeh2000.o \
freezaerl_koop2000.o rhoice_heymsfield2010.o setupvf_heymsfield2010.o \
freezaerl_mohler2010.o setupvdry.o calcrs.o vaporp_h2o_goff1946.o \
pheat.o planck.o bhmie.o mie.o vaporp_h2so4_ayers1980.o \
sulfate_utils.o sulfnuc.o sulfnucrate.o hetnucl.o newstate_calc.o

CARMA_DOC = carma_constants_mod.html carma_types_mod.html \
carmagroup_mod.html carmagas_mod.html carmaelement_mod.html carmasolute_mod.html \
carmastate_mod.html carma_mod.html coagl.html coagp.html csolve.html setupbins.html \
setupatm.html setupvf.html microslow.html newstate.html prestep.html setupckern.html \
setupcoag.html smallconc.html step.html supersat.html vaporp.html vertadv.html \
vertdif.html vertical.html versol.html rhopart.html psolve.html zeromicro.html \
nsubsteps.html setupgrow.html setupgkern.html setupnuc.html growevapl.html microfast.html \
gsolve.html actdropl.html freezglaerl_murray2010.html growp.html downgxfer.html \
gasexchange.html melticel.html upgxfer.html freezdropl.html \
downgevapply.html evapp.html evap_poly.html evap_mono.html \
evap_ingrp.html tsolve.html miess.html vaporp_h2o_buck1981.html wetr.html \
vaporp_h2o_murphy2005.html  maxconc.html setupbdif.html setupvf_std.html \
setupvf_std_shape.html totalcondensate.html versub.html freezaerl_tabazadeh2000.html \
freezaerl_koop2000.html rhoice_heymsfield2010.html setupvf_heymsfield2010.html \
freezaerl_mohler2010.html setupvdry.html calcrs.html  vaporp_h2o_goff1946.html \
pheat.html planck.html bhmie.html mie.html vaporp_h2so4_ayers1980.html \
sulfate_utils.html sulfnuc.html sulfnucrate.html hetnucl.html newstate_calc.html

carma_precision_mod.o : carma_precision_mod.F90
	$(FORTRAN) $(FFLAGS) -c $<

carma_enums_mod.o : carma_enums_mod.F90
	$(FORTRAN) $(FFLAGS) -c $<

carma_constants_mod.o : carma_constants_mod.F90 carma_precision_mod.mod 
	$(FORTRAN) $(FFLAGS) -c $<
	
carma_types_mod.o : carma_types_mod.F90 carma_constants_mod.mod carma_precision_mod.mod
	$(FORTRAN) $(FFLAGS) -c $<

carmagroup_mod.o : carmagroup_mod.F90 carma_globaer.h carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod
	$(FORTRAN) $(FFLAGS) -c $<

carmagas_mod.o : carmagas_mod.F90 carma_globaer.h carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod
	$(FORTRAN) $(FFLAGS) -c $<

carmaelement_mod.o : carmaelement_mod.F90 carma_globaer.h carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod
	$(FORTRAN) $(FFLAGS) -c $<

carmasolute_mod.o : carmasolute_mod.F90 carma_globaer.h carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod
	$(FORTRAN) $(FFLAGS) -c $<

carmastate_mod.o : carmastate_mod.F90 carma_globaer.h carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod
	$(FORTRAN) $(FFLAGS) -c $<

carma_mod.o : carma_mod.F90 carma_globaer.h carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod
	$(FORTRAN) $(FFLAGS) -c $<
	

# The following files make use of the preprocessor to map the old CARMA names to the new CARMA strucutre members via carma_globaer.h. Some
# compilers (e.g. Portland Group) have versions of the cpp that they use with Fortran that do not properly handle recursion. Because of them,
# we first invoke a conforming cpp and then compile the Fortran file.
#
# NOTE: The ifort and AIX compilers do not have this problem. 	
actdropl.o : actdropl.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(FFLAGS) -c $<

bhmie.o : bhmie.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(FFLAGS) -c $<

calcrs.o : calcrs.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(FFLAGS) -c $<

coagl.o : coagl.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(FFLAGS) -c $<

coagp.o : coagp.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(FFLAGS) -c $<

csolve.o : csolve.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(FFLAGS) -c $<

detrain.o : detrain.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(FFLAGS) -c $<

downgevapply.o : downgevapply.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(FFLAGS) -c $<

downgxfer.o : downgxfer.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(FFLAGS) -c $<

evap_ingrp.o : evap_ingrp.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(FFLAGS) -c $<

evap_mono.o : evap_mono.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(FFLAGS) -c $<

evap_poly.o : evap_poly.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(FFLAGS) -c $<

evapp.o : evapp.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(FFLAGS) -c $<

freezaerl_tabazadeh2000.o : freezaerl_tabazadeh2000.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(FFLAGS) -c $<

freezaerl_koop2000.o : freezaerl_koop2000.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(FFLAGS) -c $<

freezaerl_mohler2010.o : freezaerl_mohler2010.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(FFLAGS) -c $<

freezglaerl_murray2010.o : freezglaerl_murray2010.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(FFLAGS) -c $<

freezdropl.o : freezdropl.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(FFLAGS) -c $<

gasexchange.o : gasexchange.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(FFLAGS) -c $<

growevapl.o : growevapl.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(FFLAGS) -c $<

growp.o : growp.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(FFLAGS) -c $<

gsolve.o : gsolve.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(FFLAGS) -c $<

hetnucl.o : hetnucl.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(FFLAGS) -c $<

maxconc.o : maxconc.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod	
	$(FORTRAN) $(FFLAGS) -c $<

melticel.o : melticel.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod	
	$(FORTRAN) $(FFLAGS) -c $<

microfast.o : microfast.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod	
	$(FORTRAN) $(FFLAGS) -c $<

microslow.o : microslow.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod	
	$(FORTRAN) $(FFLAGS) -c $<

mie.o : mie.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod	
	$(FORTRAN) $(FFLAGS) -c $<

miess.o : miess.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod	
	$(FORTRAN) $(FFLAGS) -c $<

newstate.o : newstate.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 	
	$(FORTRAN) $(FFLAGS) -c $<

newstate_calc.o : newstate_calc.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 	
	$(FORTRAN) $(FFLAGS) -c $<

nsubsteps.o : nsubsteps.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 	
	$(FORTRAN) $(FFLAGS) -c $<

pheat.o : pheat.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod planck.mod	
	$(FORTRAN) $(FFLAGS) -c $<

planck.o : planck.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 	
	$(FORTRAN) $(FFLAGS) -c $<

prestep.o : prestep.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 	
	$(FORTRAN) $(FFLAGS) -c $<

psolve.o : psolve.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 	
	$(FORTRAN) $(FFLAGS) -c $<

rhoice_heymsfield2010.o : rhoice_heymsfield2010.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 	
	$(FORTRAN) $(FFLAGS) -c $<

rhopart.o : rhopart.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod sulfate_utils.mod wetr.mod	
	$(FORTRAN) $(FFLAGS) -c $<

setupatm.o : setupatm.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(FFLAGS) -c $<

setupbdif.o : setupbdif.F90 carma_globaer.h carma_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(FFLAGS) -c $<
	
setupbins.o : setupbins.F90 carma_globaer.h carma_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(FFLAGS) -c $<
	
setupckern.o : setupckern.F90	carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(FFLAGS) -c $<

setupcoag.o : setupcoag.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 	
	$(FORTRAN) $(FFLAGS) -c $<

setupgkern.o : setupgkern.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod sulfate_utils.mod
	$(FORTRAN) $(FFLAGS) -c $<

setupgrow.o : setupgrow.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(FFLAGS) -c $<

setupnuc.o : setupnuc.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(FFLAGS) -c $<

setupvf.o : setupvf.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(FFLAGS) -c $<

setupvdry.o : setupvdry.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(FFLAGS) -c $<

setupvf_heymsfield2010.o : setupvf_heymsfield2010.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 	
	$(FORTRAN) $(FFLAGS) -c $<

setupvf_std.o : setupvf_std.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(FFLAGS) -c $<

setupvf_std_shape.o : setupvf_std_shape.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(FFLAGS) -c $<

smallconc.o : smallconc.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 	
	$(FORTRAN) $(FFLAGS) -c $<

step.o : step.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 	
	$(FORTRAN) $(FFLAGS) -c $<

sulfate_utils.o : sulfate_utils.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 	
	$(FORTRAN) $(FFLAGS) -c $<

sulfnuc.o : sulfnuc.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 	
	$(FORTRAN) $(FFLAGS) -c $<

sulfnucrate.o : sulfnucrate.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 	
	$(FORTRAN) $(FFLAGS) -c $<

supersat.o : supersat.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 	
	$(FORTRAN) $(FFLAGS) -c $<

totalcondensate.o : totalcondensate.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 	
	$(FORTRAN) $(FFLAGS) -c $<

tsolve.o : tsolve.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 	
	$(FORTRAN) $(FFLAGS) -c $<

upgxfer.o : upgxfer.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(FFLAGS) -c $<

vaporp.o : vaporp.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 	
	$(FORTRAN) $(FFLAGS) -c $<

vaporp_h2so4_ayers1980.o : vaporp_h2so4_ayers1980.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 	
	$(FORTRAN) $(FFLAGS) -c $<

vaporp_h2o_buck1981.o : vaporp_h2o_buck1981.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 	
	$(FORTRAN) $(FFLAGS) -c $<

vaporp_h2o_goff1946.o : vaporp_h2o_goff1946.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 	
	$(FORTRAN) $(FFLAGS) -c $<

vaporp_h2o_murphy2005.o : vaporp_h2o_murphy2005.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 	
	$(FORTRAN) $(FFLAGS) -c $<

versol.o : versol.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 		
	$(FORTRAN) $(FFLAGS) -c $<

versub.o : versub.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 		
	$(FORTRAN) $(FFLAGS) -c $<

vertadv.o : vertadv.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 		
	$(FORTRAN) $(FFLAGS) -c $<
	
vertdif.o : vertdif.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 		
	$(FORTRAN) $(FFLAGS) -c $<

vertical.o : vertical.F90	carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 	
	$(FORTRAN) $(FFLAGS) -c $<

wetr.o : wetr.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 		
	$(FORTRAN) $(FFLAGS) -c $<

zeromicro.o : zeromicro.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 		
	$(FORTRAN) $(FFLAGS) -c $<
