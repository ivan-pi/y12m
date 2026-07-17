# SPDX-License-Identifier: GPL-2.0-only
cmake_minimum_required(VERSION 3.14)

project(y12m
    DESCRIPTION "Sparse Direct Solver"
    LANGUAGES Fortran C)

option(BUILD_EXAMPLES "Build the example programs" ON)
option(WITH_NIST_MMIO "Build the Matrix Market driver (y12m_mm)" OFF)
option(Y12M_FORCE_DOUBLE_DOUBLE
    "Force double-double residual accumulation in y12mff even when the compiler supports quad precision" OFF)

if(WITH_NIST_MMIO)
    include(cmake/FetchMMIO.cmake)
endif()

# ---------------------------------------------------------------------------
# Residual accumulator for y12mff (double precision iterative refinement).
#
# y12mff accumulates residuals in higher than double precision: native quad
# (selected_real_kind(33)) when the compiler provides it, otherwise a
# double-double compensated accumulation built on IEEE_FMA.  Detect the
# capabilities here and fix the strategy at configure time via the
# Y12M_ACCUM_QUAD / Y12M_ACCUM_DD preprocessor macros (see src/y12mff.F90).
# ---------------------------------------------------------------------------
include(CheckFortranSourceCompiles)

check_fortran_source_compiles("
program conftest
integer, parameter :: qp = selected_real_kind(33)
real(qp) :: t
t = 1.0_qp
end program conftest
" Y12M_HAVE_QUAD SRC_EXT F90)

if(Y12M_HAVE_QUAD)
    message(STATUS "Fortran compiler supports quad precision (selected_real_kind(33)): yes")
else()
    message(STATUS "Fortran compiler supports quad precision (selected_real_kind(33)): no")
endif()

if(Y12M_FORCE_DOUBLE_DOUBLE OR NOT Y12M_HAVE_QUAD)
    check_fortran_source_compiles("
program conftest
use, intrinsic :: ieee_arithmetic, only: ieee_fma
real(kind(1.0d0)) :: t
t = ieee_fma(1.0d0, 1.0d0, 1.0d0)
end program conftest
" Y12M_HAVE_IEEE_FMA SRC_EXT F90)
    if(NOT Y12M_HAVE_IEEE_FMA)
        message(FATAL_ERROR
            "y12mff needs the double-double residual accumulator "
            "(quad precision unavailable or Y12M_FORCE_DOUBLE_DOUBLE=ON), "
            "which requires the Fortran 2018 IEEE_FMA intrinsic. "
            "Use a compiler that provides either quad precision or IEEE_FMA "
            "(e.g. gfortran >= 13).")
    endif()
    set(Y12M_ACCUM_MACRO Y12M_ACCUM_DD)
    if(Y12M_FORCE_DOUBLE_DOUBLE)
        message(STATUS "y12mff residual accumulator: double-double (forced by Y12M_FORCE_DOUBLE_DOUBLE)")
    else()
        message(STATUS "y12mff residual accumulator: double-double (quad precision not available)")
    endif()
else()
    set(Y12M_ACCUM_MACRO Y12M_ACCUM_QUAD)
    message(STATUS "y12mff residual accumulator: quad precision (selected_real_kind(33))")
endif()

# Collect all legacy fixed-form Fortran sources.
file(GLOB LEGACY_SOURCES "${CMAKE_SOURCE_DIR}/src/legacy/*.f")

# Build a static library from the legacy routines plus the Fortran module.
# These sources are pre-existing; no extra diagnostic flags are applied.
add_library(y12m_legacy STATIC ${LEGACY_SOURCES}
    "${CMAKE_SOURCE_DIR}/src/y12m.F90"
    "${CMAKE_SOURCE_DIR}/src/y12mff.F90")

# The accumulator selection must be identical in y12mff.F90 (the algorithm)
# and y12m.F90 (the y12mff_has_quad / y12mff_accum_kind constants), so the
# macro is applied to the whole target.
target_compile_definitions(y12m_legacy PRIVATE ${Y12M_ACCUM_MACRO})

# The double-double compensation terms in y12mff.F90 are algebraically zero;
# value-unsafe floating point modes (-ffast-math, -Ofast, ifx's default
# -fp-model=fast) are licensed to eliminate them, silently degrading the
# accumulation to plain double precision.  Pin value-safe semantics for this
# one file regardless of the global flags.
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    set(Y12M_VALUE_SAFE_FLAGS -fno-unsafe-math-optimizations)
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "LLVMFlang|Flang|Clang")
    set(Y12M_VALUE_SAFE_FLAGS -fno-fast-math)
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    if(WIN32)
        set(Y12M_VALUE_SAFE_FLAGS /fp:precise)
    else()
        set(Y12M_VALUE_SAFE_FLAGS -fp-model=precise)
    endif()
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC|PGI")
    set(Y12M_VALUE_SAFE_FLAGS -Kieee)
else()
    set(Y12M_VALUE_SAFE_FLAGS "")
endif()
if(Y12M_VALUE_SAFE_FLAGS)
    set_source_files_properties("${CMAKE_SOURCE_DIR}/src/y12mff.F90"
        PROPERTIES COMPILE_OPTIONS "${Y12M_VALUE_SAFE_FLAGS}")
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
