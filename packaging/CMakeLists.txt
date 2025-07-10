include(CMakePackageConfigHelpers)

# install the binary and static library
install(
  TARGETS 
   carma carma carma_object
  EXPORT 
    carma_Exports
  LIBRARY DESTINATION ${INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

# install the mod files
install(
  DIRECTORY 
    ${CARMA_MOD_DIR}/
  DESTINATION 
    ${CARMA_INSTALL_MOD_DIR}
  FILES_MATCHING PATTERN "*.mod"
)

# install the cmake config files
set(cmake_config_install_location ${INSTALL_PREFIX}/cmake/carma-${PROJECT_VERSION})

install(
  EXPORT 
    carma_Exports
  DESTINATION 
    ${cmake_config_install_location}
  NAMESPACE musica::
)

configure_package_config_file(
  "${PROJECT_SOURCE_DIR}/cmake/carmaConfig.cmake.in"
  "${PROJECT_BINARY_DIR}/carmaConfig.cmake"
  INSTALL_DESTINATION
    ${cmake_config_install_location}
)

write_basic_package_version_file(
  "${PROJECT_BINARY_DIR}/carmaConfigVersion.cmake"
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion
)

install(
  FILES
    ${PROJECT_BINARY_DIR}/carmaConfig.cmake
    ${PROJECT_BINARY_DIR}/carmaConfigVersion.cmake
  DESTINATION
    ${INSTALL_PREFIX}/cmake/carma-${PROJECT_VERSION}
)

######################################################################
# uninstall target

# https://gitlab.kitware.com/cmake/community/-/wikis/FAQ#can-i-do-make-uninstall-with-cmake
if(NOT TARGET uninstall)
  configure_file(
    "${PROJECT_SOURCE_DIR}/cmake/cmake_uninstall.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
    IMMEDIATE @ONLY)

  add_custom_target(uninstall
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
endif()