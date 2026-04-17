# SPDX-License-Identifier: GPL-2.0-only
cmake_minimum_required(VERSION 3.14)

project(y12m_tests Fortran)

# Collect all legacy fixed-form Fortran sources.
file(GLOB LEGACY_SOURCES "${CMAKE_SOURCE_DIR}/src/legacy/*.f")

# src/legacy/y12mhf.f is missing its END statement and does not compile.
# Exclude it and use the complete replacement from compat/.
list(REMOVE_ITEM LEGACY_SOURCES
    "${CMAKE_SOURCE_DIR}/src/legacy/y12mhf.f")
list(APPEND LEGACY_SOURCES
    "${CMAKE_SOURCE_DIR}/compat/y12mhf.f")

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

enable_testing()
add_subdirectory(test)
