# SPDX-License-Identifier: GPL-2.0-only
cmake_minimum_required(VERSION 3.14)

project(y12m
    DESCRIPTION "Sparse Direct Solver"
    LANGUAGES Fortran C)

option(BUILD_EXAMPLES "Build the example programs" ON)
option(WITH_NIST_MMIO "Build the Matrix Market driver (y12m_mm)" OFF)

if(WITH_NIST_MMIO)
    include(cmake/FetchMMIO.cmake)
endif()

# Residual accumulator for y12mff (quad precision or double-double); sets
# Y12MFF_COMPILE_DEFINITIONS and Y12MFF_COMPILE_OPTIONS and prints the
# selection.  Options: Y12M_FORCE_DOUBLE_DOUBLE, Y12M_USE_IEEE_FMA.
include(cmake/Y12mffAccumulator.cmake)

# Collect all legacy fixed-form Fortran sources.
file(GLOB LEGACY_SOURCES "${CMAKE_SOURCE_DIR}/src/legacy/*.f")

# Build a static library from the legacy routines plus the Fortran module.
# These sources are pre-existing; no extra diagnostic flags are applied.
add_library(y12m_legacy STATIC ${LEGACY_SOURCES}
    "${CMAKE_SOURCE_DIR}/src/y12m.F90"
    "${CMAKE_SOURCE_DIR}/src/y12mff.F90")

# The accumulator selection must be identical in y12mff.F90 (the algorithm)
# and y12m.F90 (the y12mff_uses_quad / y12mff_accumulator_kind constants),
# so the macros are applied to the whole target.  The IEEE-compliance flags
# only concern the double-double arithmetic in y12mff.F90.
target_compile_definitions(y12m_legacy PRIVATE ${Y12MFF_COMPILE_DEFINITIONS})
if(Y12MFF_COMPILE_OPTIONS)
    set_source_files_properties("${CMAKE_SOURCE_DIR}/src/y12mff.F90"
        PROPERTIES COMPILE_OPTIONS "${Y12MFF_COMPILE_OPTIONS}")
endif()

# Place generated .mod files in a dedicated directory and expose it to
# any target that links against this library.
set(FORTRAN_MOD_DIR "${CMAKE_BINARY_DIR}/include")
set_target_properties(y12m_legacy PROPERTIES
    Fortran_MODULE_DIRECTORY "${FORTRAN_MOD_DIR}")
target_include_directories(y12m_legacy INTERFACE "${FORTRAN_MOD_DIR}")
target_include_directories(y12m_legacy INTERFACE "${CMAKE_SOURCE_DIR}/include")

include(CTest)

if(BUILD_TESTING)
    find_package(LAPACK)
    add_subdirectory(test)
endif()

if(BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()
