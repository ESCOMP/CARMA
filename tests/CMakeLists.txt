################################################################################
# Test utilities

include(test_util)

################################################################################
# Tests

create_standard_test(NAME version SOURCES atmosphere_mod.F90 carma_versiontest.F90)
create_standard_test(NAME bc2g SOURCES atmosphere_mod.F90 carma_bc2gtest.F90)
create_standard_test(NAME bcoc SOURCES atmosphere_mod.F90 carma_bcoctest.F90)
create_standard_test(NAME ce SOURCES atmosphere_mod.F90 carma_cetest.F90)
create_standard_test(NAME coag SOURCES atmosphere_mod.F90 carma_coagtest.F90)
create_standard_test(NAME drydep SOURCES atmosphere_mod.F90 carma_drydeptest.F90)
create_standard_test(NAME fall SOURCES atmosphere_mod.F90 carma_falltest.F90)
create_standard_test(NAME fractalmicro SOURCES atmosphere_mod.F90 carma_fractalmicrotest.F90)
create_standard_test(NAME fractaloptics SOURCES atmosphere_mod.F90 carma_fractalopticstest.F90)
create_standard_test(NAME growclr SOURCES atmosphere_mod.F90 carma_growclrtest.F90)
create_standard_test(NAME growin SOURCES atmosphere_mod.F90 carma_growintest.F90)
create_standard_test(NAME growsub SOURCES atmosphere_mod.F90 carma_growsubtest.F90)
create_standard_test(NAME grow SOURCES atmosphere_mod.F90 carma_growtest.F90)
create_standard_test(NAME init SOURCES atmosphere_mod.F90 carma_inittest.F90)
create_standard_test(NAME kappawetr SOURCES atmosphere_mod.F90 carma_kappawetrtest.F90)
create_standard_test(NAME mie SOURCES atmosphere_mod.F90 carma_mietest.F90)
create_standard_test(NAME nuc2 SOURCES atmosphere_mod.F90 carma_nuc2test.F90)
create_standard_test(NAME nuc SOURCES atmosphere_mod.F90 carma_nuctest.F90)
create_standard_test(NAME pheat SOURCES atmosphere_mod.F90 carma_pheattest.F90)
create_standard_test(NAME scfall SOURCES atmosphere_mod.F90 carma_scfalltest.F90)
create_standard_test(NAME sigmadrydep SOURCES atmosphere_mod.F90 carma_sigmadrydeptest.F90)
create_standard_test(NAME sigmafall SOURCES atmosphere_mod.F90 carma_sigmafalltest.F90)
create_standard_test(NAME sulfate_ccdon SOURCES atmosphere_mod.F90 carma_sulfate_ccdon_test.F90)
create_standard_test(NAME sulfate_vehkamaki SOURCES atmosphere_mod.F90 carma_sulfate_vehkamaki_test.F90)
create_standard_test(NAME sulfate SOURCES atmosphere_mod.F90 carma_sulfatetest.F90)
create_standard_test(NAME sulfhet_vehkamaki SOURCES atmosphere_mod.F90 carma_sulfhet_vehkamaki_test.F90)
create_standard_test(NAME sulfhet SOURCES atmosphere_mod.F90 carma_sulfhettest.F90)
create_standard_test(NAME swell SOURCES atmosphere_mod.F90 carma_swelltest.F90)
create_standard_test(NAME vdif SOURCES atmosphere_mod.F90 carma_vdiftest.F90)

################################################################################