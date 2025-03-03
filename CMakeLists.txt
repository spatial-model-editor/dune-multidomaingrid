cmake_minimum_required(VERSION 3.16)
project(dune-multidomaingrid CXX)

if(NOT (dune-common_DIR OR dune-common_ROOT OR
      "${CMAKE_PREFIX_PATH}" MATCHES ".*dune-common.*"))
    string(REPLACE  ${CMAKE_PROJECT_NAME} dune-common dune-common_DIR
      ${PROJECT_BINARY_DIR})
endif()

#find dune-common and set the module path
find_package(dune-common)
list(APPEND CMAKE_MODULE_PATH ${dune-common_MODULE_PATH}
  "${PROJECT_SOURCE_DIR}/cmake/modules")

#include the dune macros
include(DuneMacros)

# deactivate global include-directories
dune_policy(SET DP_DEFAULT_INCLUDE_DIRS NEW)

# start a dune project with information from dune.module
dune_project()

# Create the main target with a Dune::MultiDomainGrid alias
dune_add_library(dunemultidomaingrid INTERFACE
  EXPORT_NAME MultiDomainGrid
  LINK_LIBRARIES ${DUNE_LIBS})

# set include directories for main library
dune_default_include_directories(dunemultidomaingrid INTERFACE)


add_subdirectory(doc)
add_subdirectory(dune)
add_subdirectory(test)

# finalize the dune project, e.g. export Dune::MultiDomainGrid target
finalize_dune_project()
