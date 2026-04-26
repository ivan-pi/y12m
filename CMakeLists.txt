# SPDX-License-Identifier: GPL-2.0-only
cmake_minimum_required(VERSION 3.14)

project(y12m
    DESCRIPTION "Sparse Direct Solver"
    LANGUAGES Fortran)

option(BUILD_EXAMPLES "Build the example programs" ON)

# Collect all legacy fixed-form Fortran sources.
file(GLOB LEGACY_SOURCES "${CMAKE_SOURCE_DIR}/src/legacy/*.f")

# Collect free-form sources in the legacy directory (e.g. y12mff.f90).
file(GLOB LEGACY_F90_SOURCES "${CMAKE_SOURCE_DIR}/src/legacy/*.f90")

# Build a static library from the legacy routines plus the Fortran module.
# These sources are pre-existing; no extra diagnostic flags are applied.
add_library(y12m_legacy STATIC ${LEGACY_SOURCES} ${LEGACY_F90_SOURCES}
    "${CMAKE_SOURCE_DIR}/src/y12m.f90")

# Place generated .mod files in a dedicated directory and expose it to
# any target that links against this library.
set(FORTRAN_MOD_DIR "${CMAKE_BINARY_DIR}/include")
set_target_properties(y12m_legacy PROPERTIES
    Fortran_MODULE_DIRECTORY "${FORTRAN_MOD_DIR}")
target_include_directories(y12m_legacy INTERFACE "${FORTRAN_MOD_DIR}")

include(CTest)

if(BUILD_TESTING)
    find_package(LAPACK)
    add_subdirectory(test)
endif()

if(BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()
