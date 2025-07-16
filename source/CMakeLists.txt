################################################################################
# CARMA source

# object library
add_library(carma_object OBJECT)

set_target_properties(carma_object PROPERTIES
  Fortran_MODULE_DIRECTORY ${CARMA_MOD_DIR}
)

target_include_directories(
  carma_object
  PUBLIC
    $<BUILD_INTERFACE:${CARMA_MOD_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/base>
    $<INSTALL_INTERFACE:${CARMA_INSTALL_MOD_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

target_link_libraries(
  carma_object
  PRIVATE
    PkgConfig::netcdff
    PkgConfig::netcdfc
    ${LAPACK_LIBRARIES}
    LAPACK::LAPACK
    ${BLAS_LIBRARIES}
)

add_library(carma $<TARGET_OBJECTS:carma_object>)
add_library(musica::carma ALIAS carma)

if(NOT BUILD_SHARED_LIBS)
  set_target_properties(carma_object PROPERTIES POSITION_INDEPENDENT_CODE ON)
endif()

set_target_properties(
  carma
  PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY ${CARMA_LIB_DIR}
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR}
)

target_link_libraries(carma PUBLIC carma_object)

message(
  STATUS
    "CARMA build include directories: ${CARMA_MOD_DIR};${CMAKE_CURRENT_SOURCE_DIR}/base"
)
message(
  STATUS
    "CARMA install include directories: ${CARMA_INSTALL_MOD_DIR};${CMAKE_INSTALL_INCLUDEDIR}"
)
target_include_directories(
  carma
  PUBLIC
    $<BUILD_INTERFACE:${CARMA_MOD_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/base>
    $<INSTALL_INTERFACE:${CARMA_INSTALL_MOD_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

configure_file(version.F90.in ${CMAKE_CURRENT_BINARY_DIR}/version.F90 @ONLY)

target_sources(
  carma_object
  PRIVATE
    base/actdropl.F90
    base/adgaquad_mod.F90
    base/adgaquad_types_mod.F90
    base/bhmie.F90
    base/calcrs.F90
    base/carma_constants_mod.F90
    base/carma_enums_mod.F90
    base/carma_mod.F90
    base/carma_precision_mod.F90
    base/carma_types_mod.F90
    base/carmaelement_mod.F90
    base/carmagas_mod.F90
    base/carmagroup_mod.F90
    base/carmasolute_mod.F90
    base/carmastate_mod.F90
    base/coagl.F90
    base/coagp.F90
    base/coremasscheck.F90
    base/csolve.F90
    base/detrain.F90
    base/downgevapply.F90
    base/downgxfer.F90
    base/evap_ingrp.F90
    base/evap_mono.F90
    base/evap_poly.F90
    base/evapp.F90
    base/fixcorecol.F90
    base/fractal_meanfield_mod.F90
    base/freezaerl_koop2000.F90
    base/freezaerl_mohler2010.F90
    base/freezaerl_tabazadeh2000.F90
    base/freezdropl.F90
    base/freezglaerl_murray2010.F90
    base/gasexchange.F90
    base/growevapl.F90
    base/growp.F90
    base/gsolve.F90
    base/hetnucl.F90
    base/hygroscopicity.F90
    base/maxconc.F90
    base/melticel.F90
    base/microfast.F90
    base/microslow.F90
    base/mie.F90
    base/miess.F90
    base/newstate.F90
    base/newstate_calc.F90
    base/nsubsteps.F90
    base/pheat.F90
    base/planck.F90
    base/prestep.F90
    base/psolve.F90
    base/rhoice_heymsfield2010.F90
    base/rhopart.F90
    base/setupatm.F90
    base/setupbdif.F90
    base/setupbins.F90
    base/setupckern.F90
    base/setupcoag.F90
    base/setupgkern.F90
    base/setupgrow.F90
    base/setupnuc.F90
    base/setupvdry.F90
    base/setupvf.F90
    base/setupvf_heymsfield2010.F90
    base/setupvf_std.F90
    base/setupvf_std_shape.F90
    base/smallconc.F90
    base/step.F90
    base/sulfate_utils.F90
    base/sulfhetnucrate.F90
    base/sulfnuc.F90
    base/sulfnucrate.F90
    base/supersat.F90
    base/totalcondensate.F90
    base/tsolve.F90
    base/upgxfer.F90
    base/vaporp.F90
    base/vaporp_h2o_buck1981.F90
    base/vaporp_h2o_goff1946.F90
    base/vaporp_h2o_murphy2005.F90
    base/vaporp_h2so4_ayers1980.F90
    base/versol.F90
    base/versub.F90
    base/vertadv.F90
    base/vertdif.F90
    base/vertical.F90
    base/wetr.F90
    base/zeromicro.F90
    ${CMAKE_CURRENT_BINARY_DIR}/version.F90
    ../tests/atmosphere_mod.F90
)