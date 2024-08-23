cmake_minimum_required(VERSION 3.13)
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

# start a dune project with information from dune.module
dune_project()

dune_enable_all_packages(VERBOSE)

dune_add_library(dunemultidomaingrid INTERFACE
  EXPORT_NAME MultiDomainGrid
  LINK_LIBRARIES Dune::Grid)

# disable the upstream testing magic
set(DUNE_TEST_MAGIC OFF)

add_subdirectory(doc)
add_subdirectory(dune)
add_subdirectory(test)

# finalize the dune project, e.g. generating config.h etc.
finalize_dune_project()
