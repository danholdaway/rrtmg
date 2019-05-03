################################################################################
# RRTMG 
################################################################################

cmake_minimum_required( VERSION 3.3.2 FATAL_ERROR )

project( rrtmg Fortran )

set( CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake;${CMAKE_MODULE_PATH})

set( ECBUILD_DEFAULT_BUILD_TYPE Release )
set( ENABLE_OS_TESTS           OFF CACHE BOOL "Disable OS tests" FORCE )
set( ENABLE_LARGE_FILE_SUPPORT OFF CACHE BOOL "Disable testing of large file support" FORCE )

include( ecbuild_system NO_POLICY_SCOPE )

ecbuild_requires_macro_version( 2.5 )

################################################################################
# Project
################################################################################

ecbuild_declare_project()

ecbuild_enable_fortran( REQUIRED )

set( RRTMG_LINKER_LANGUAGE Fortran )

################################################################################
# Dependencies
################################################################################

# NetCDF
find_package( NetCDF REQUIRED COMPONENTS Fortran F90 C )
include_directories( ${NETCDF_INCLUDE_DIR} )

################################################################################
# Export package info
################################################################################

list( APPEND RRTMG_TPLS NetCDF)

set( RRTMG_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR} 
                          ${CMAKE_CURRENT_SOURCE_DIR}/src
                          ${CMAKE_CURRENT_SOURCE_DIR}/src/rrtmg_lw/gcm_model/modules
                          ${CMAKE_Fortran_MODULE_DIRECTORY})

set( RRTMG_LIBRARIES rrtmg )

get_directory_property( RRTMG_DEFINITIONS COMPILE_DEFINITIONS )

foreach( _tpl ${RRTMG_TPLS} )
  string( TOUPPER ${_tpl} TPL )
  list( APPEND RRTMG_EXTRA_DEFINITIONS   ${${TPL}_DEFINITIONS}  ${${TPL}_TPL_DEFINITIONS}  )
  list( APPEND RRTMG_EXTRA_INCLUDE_DIRS  ${${TPL}_INCLUDE_DIRS} ${${TPL}_TPL_INCLUDE_DIRS} )
  list( APPEND RRTMG_EXTRA_LIBRARIES     ${${TPL}_LIBRARIES}    ${${TPL}_TPL_LIBRARIES}    )
endforeach()

################################################################################
# Sources
################################################################################

include( rrtmg_compiler_flags )

include_directories( ${RRTMG_INCLUDE_DIRS}
                     ${RRTMG_EXTRA_INCLUDE_DIRS} )

add_subdirectory( mains )
add_subdirectory( src )

################################################################################
# Finalise configuration
################################################################################

# prepares a tar.gz of the sources and/or binaries
ecbuild_install_project( NAME rrtmg )

# print the summary of the configuration
ecbuild_print_summary()

ecbuild_add_test( TARGET rrtmg_coding_norms
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/cpplint.py
                  ARGS --quiet --recursive ${CMAKE_CURRENT_SOURCE_DIR} )

