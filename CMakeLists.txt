# SPDX-License-Identifier: GPL-2.0-only
cmake_minimum_required(VERSION 3.14)

project(y12m
    DESCRIPTION "Sparse Direct Solver"
    LANGUAGES Fortran)

# Collect all legacy fixed-form Fortran sources.
file(GLOB LEGACY_SOURCES "${CMAKE_SOURCE_DIR}/src/legacy/*.f")

# Build a static library from the legacy routines plus the Fortran module.
# These sources are pre-existing; no extra diagnostic flags are applied.
add_library(y12m_legacy STATIC ${LEGACY_SOURCES}
    "${CMAKE_SOURCE_DIR}/src/y12m.f90")

# Place generated .mod files in a dedicated directory and expose it to
# any target that links against this library.
set(FORTRAN_MOD_DIR "${CMAKE_BINARY_DIR}/include")
set_target_properties(y12m_legacy PROPERTIES
    Fortran_MODULE_DIRECTORY "${FORTRAN_MOD_DIR}")
target_include_directories(y12m_legacy INTERFACE "${FORTRAN_MOD_DIR}")


find_package(LAPACK)

enable_testing()
add_subdirectory(test)
add_subdirectory(examples)
