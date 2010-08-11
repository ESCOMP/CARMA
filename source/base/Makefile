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
freezaerl_koop2000.o rhoice_heymsfield2010.o setupvf_heymsfield2010.o

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
freezaerl_koop2000.html rhoice_heymsfield2010.html setupvf_heymsfield2010.html

carma_precision_mod.o : carma_precision_mod.F90
	$(FORTRAN) $(CPPFLAGS) -c $<

carma_enums_mod.o : carma_enums_mod.F90
	$(FORTRAN) $(CPPFLAGS) -c $<

carma_constants_mod.o : carma_constants_mod.F90 carma_precision_mod.mod 
	$(FORTRAN) $(CPPFLAGS) -c $<
	
carma_types_mod.o : carma_types_mod.F90 carma_constants_mod.mod carma_precision_mod.mod
	$(FORTRAN) $(CPPFLAGS) -c $<

carmagroup_mod.o : carmagroup_mod.F90 carma_globaer.h carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod
	$(FORTRAN) $(CPPFLAGS) -c $<

carmagas_mod.o : carmagas_mod.F90 carma_globaer.h carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod
	$(FORTRAN) $(CPPFLAGS) -c $<

carmaelement_mod.o : carmaelement_mod.F90 carma_globaer.h carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod
	$(FORTRAN) $(CPPFLAGS) -c $<

carmasolute_mod.o : carmasolute_mod.F90 carma_globaer.h carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod
	$(FORTRAN) $(CPPFLAGS) -c $<

carmastate_mod.o : carmastate_mod.F90 carma_globaer.h carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod
	$(FORTRAN) $(CPPFLAGS) -c $<

carma_mod.o : carma_mod.F90 carma_globaer.h carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod
	$(FORTRAN) $(CPPFLAGS) -c $<
	
	
actdropl.o : actdropl.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(CPPFLAGS) -c $<

coagl.o : coagl.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(CPPFLAGS) -c $<

coagp.o : coagp.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(CPPFLAGS) -c $<

csolve.o : csolve.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(CPPFLAGS) -c $<

detrain.o : detrain.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(CPPFLAGS) -c $<

downgevapply.o : downgevapply.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(CPPFLAGS) -c $<

downgxfer.o : downgxfer.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(CPPFLAGS) -c $<

evap_ingrp.o : evap_ingrp.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(CPPFLAGS) -c $<

evap_mono.o : evap_mono.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(CPPFLAGS) -c $<

evap_poly.o : evap_poly.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(CPPFLAGS) -c $<

evapp.o : evapp.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(CPPFLAGS) -c $<

freezaerl_tabazadeh2000.o : freezaerl_tabazadeh2000.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(CPPFLAGS) -c $<

freezaerl_koop2000.o : freezaerl_koop2000.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(CPPFLAGS) -c $<

freezglaerl_murray2010.o : freezglaerl_murray2010.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(CPPFLAGS) -c $<

freezdropl.o : freezdropl.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(CPPFLAGS) -c $<

gasexchange.o : gasexchange.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(CPPFLAGS) -c $<

growevapl.o : growevapl.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(CPPFLAGS) -c $<

growp.o : growp.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(CPPFLAGS) -c $<

gsolve.o : gsolve.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(CPPFLAGS) -c $<

maxconc.o : maxconc.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod	
	$(FORTRAN) $(CPPFLAGS) -c $<

melticel.o : melticel.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod	
	$(FORTRAN) $(CPPFLAGS) -c $<

microfast.o : microfast.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod	
	$(FORTRAN) $(CPPFLAGS) -c $<

microslow.o : microslow.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod	
	$(FORTRAN) $(CPPFLAGS) -c $<

miess.o : miess.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod	
	$(FORTRAN) $(CPPFLAGS) -c $<

newstate.o : newstate.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 	
	$(FORTRAN) $(CPPFLAGS) -c $<

nsubsteps.o : nsubsteps.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 	
	$(FORTRAN) $(CPPFLAGS) -c $<

prestep.o : prestep.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 	
	$(FORTRAN) $(CPPFLAGS) -c $<

psolve.o : psolve.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 	
	$(FORTRAN) $(CPPFLAGS) -c $<

rhoice_heymsfield2010.o : rhoice_heymsfield2010.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 	
	$(FORTRAN) $(CPPFLAGS) -c $<

rhopart.o : rhopart.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 	
	$(FORTRAN) $(CPPFLAGS) -c $<

setupatm.o : setupatm.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(CPPFLAGS) -c $<

setupbdif.o : setupbdif.F90 carma_globaer.h carma_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(CPPFLAGS) -c $<
	
setupbins.o : setupbins.F90 carma_globaer.h carma_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(CPPFLAGS) -c $<
	
setupckern.o : setupckern.F90	carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(CPPFLAGS) -c $<

setupcoag.o : setupcoag.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 	
	$(FORTRAN) $(CPPFLAGS) -c $<

setupgkern.o : setupgkern.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(CPPFLAGS) -c $<

setupgrow.o : setupgrow.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(CPPFLAGS) -c $<

setupnuc.o : setupnuc.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(CPPFLAGS) -c $<

setupvf.o : setupvf.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(CPPFLAGS) -c $<

setupvf_heymsfield2010.o : setupvf_heymsfield2010.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 	
	$(FORTRAN) $(CPPFLAGS) -c $<

setupvf_std.o : setupvf_std.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(CPPFLAGS) -c $<

setupvf_std_shape.o : setupvf_std_shape.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(CPPFLAGS) -c $<

smallconc.o : smallconc.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 	
	$(FORTRAN) $(CPPFLAGS) -c $<

step.o : step.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 	
	$(FORTRAN) $(CPPFLAGS) -c $<

supersat.o : supersat.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 	
	$(FORTRAN) $(CPPFLAGS) -c $<

totalcondensate.o : totalcondensate.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 	
	$(FORTRAN) $(CPPFLAGS) -c $<

tsolve.o : tsolve.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 	
	$(FORTRAN) $(CPPFLAGS) -c $<

upgxfer.o : upgxfer.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 
	$(FORTRAN) $(CPPFLAGS) -c $<

vaporp.o : vaporp.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 	
	$(FORTRAN) $(CPPFLAGS) -c $<

vaporp_h2o_buck1981.o : vaporp_h2o_buck1981.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 	
	$(FORTRAN) $(CPPFLAGS) -c $<

vaporp_h2o_murphy2005.o : vaporp_h2o_murphy2005.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 	
	$(FORTRAN) $(CPPFLAGS) -c $<

versol.o : versol.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 		
	$(FORTRAN) $(CPPFLAGS) -c $<

versub.o : versub.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 		
	$(FORTRAN) $(CPPFLAGS) -c $<

vertadv.o : vertadv.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 		
	$(FORTRAN) $(CPPFLAGS) -c $<
	
vertdif.o : vertdif.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 		
	$(FORTRAN) $(CPPFLAGS) -c $<

vertical.o : vertical.F90	carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 	
	$(FORTRAN) $(CPPFLAGS) -c $<

wetr.o : wetr.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 		
	$(FORTRAN) $(CPPFLAGS) -c $<

zeromicro.o : zeromicro.F90 carma_globaer.h carma_mod.mod carmastate_mod.mod carma_types_mod.mod carma_constants_mod.mod carma_enums_mod.mod carma_precision_mod.mod 		
	$(FORTRAN) $(CPPFLAGS) -c $<
