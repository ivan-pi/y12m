# SPDX-License-Identifier: GPL-2.0-only
# Assisted-by: Claude Code:claude-fable-5
#
# Residual accumulator selection for y12mff (double precision iterative
# refinement, src/y12mff.F90).
#
# y12mff accumulates residuals in higher than double precision: native quad
# (selected_real_kind(33)) when the compiler provides it, otherwise a
# double-double compensated accumulation.  This script detects the compiler
# capabilities, honours the user options, and reports the outcome.  It sets:
#
#   Y12MFF_COMPILE_DEFINITIONS  preprocessor macros for the y12m_legacy
#                               target (Y12M_ACCUM_QUAD or Y12M_ACCUM_DD,
#                               plus Y12M_USE_IEEE_FMA when requested)
#   Y12MFF_COMPILE_OPTIONS      per-source flags for src/y12mff.F90 that pin
#                               IEEE-compliant evaluation of the
#                               double-double arithmetic (empty when the
#                               quad path is selected)
#
# The double-double compensation terms are algebraically zero, so
# value-unsafe modes (-ffast-math/-Ofast, ifx's default -fp-model=fast) are
# licensed to eliminate them, silently degrading the accumulator to plain
# double.  In addition, the default Dekker product splitting is broken by
# FMA contraction (gfortran/flang default -ffp-contract=fast on fma-capable
# hardware): fusing t - a with t = splitter*a skips exactly the rounding
# the split depends on.  Hence value-safe flags whenever the double-double
# path is compiled, and contraction disabled unless IEEE_FMA is used.

option(Y12M_FORCE_DOUBLE_DOUBLE
    "Force double-double residual accumulation in y12mff even when the compiler supports quad precision" OFF)
option(Y12M_USE_IEEE_FMA
    "Use the Fortran 2018 IEEE_FMA intrinsic for the double-double product splitting in y12mff (default: Dekker's algorithm). Worthwhile only when the compiler inlines it to a hardware instruction." OFF)

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

set(Y12MFF_COMPILE_DEFINITIONS "")
set(Y12MFF_COMPILE_OPTIONS "")

if(Y12M_USE_IEEE_FMA)
    check_fortran_source_compiles("
program conftest
use, intrinsic :: ieee_arithmetic, only: ieee_fma
real(kind(1.0d0)) :: t
t = ieee_fma(1.0d0, 1.0d0, 1.0d0)
end program conftest
" Y12M_HAVE_IEEE_FMA SRC_EXT F90)
    if(NOT Y12M_HAVE_IEEE_FMA)
        message(FATAL_ERROR
            "Y12M_USE_IEEE_FMA=ON but the Fortran compiler does not support "
            "the Fortran 2018 IEEE_FMA intrinsic (gfortran needs version 13 "
            "or newer).  Turn the option off to use Dekker's algorithm.")
    endif()
    list(APPEND Y12MFF_COMPILE_DEFINITIONS Y12M_USE_IEEE_FMA)
endif()

if(Y12M_FORCE_DOUBLE_DOUBLE OR NOT Y12M_HAVE_QUAD)
    list(APPEND Y12MFF_COMPILE_DEFINITIONS Y12M_ACCUM_DD)
    if(Y12M_USE_IEEE_FMA)
        set(_y12mff_twoprod "IEEE_FMA")
    else()
        set(_y12mff_twoprod "Dekker splitting")
    endif()
    if(Y12M_FORCE_DOUBLE_DOUBLE)
        message(STATUS "y12mff residual accumulator: double-double via ${_y12mff_twoprod} (forced by Y12M_FORCE_DOUBLE_DOUBLE)")
    else()
        message(STATUS "y12mff residual accumulator: double-double via ${_y12mff_twoprod} (quad precision not available)")
    endif()
    unset(_y12mff_twoprod)

    # Value-safe evaluation for the compensated arithmetic.
    if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
        list(APPEND Y12MFF_COMPILE_OPTIONS -fno-unsafe-math-optimizations)
    elseif(CMAKE_Fortran_COMPILER_ID MATCHES "LLVMFlang|Flang|Clang")
        list(APPEND Y12MFF_COMPILE_OPTIONS -fno-fast-math)
    elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
        if(WIN32)
            list(APPEND Y12MFF_COMPILE_OPTIONS /fp:precise)
        else()
            list(APPEND Y12MFF_COMPILE_OPTIONS -fp-model=precise)
        endif()
    elseif(CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC|PGI")
        list(APPEND Y12MFF_COMPILE_OPTIONS -Kieee)
    endif()

    # No FMA contraction for the Dekker splitting (the IEEE_FMA variant is
    # immune, so contraction may stay enabled there).
    if(NOT Y12M_USE_IEEE_FMA)
        if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU|LLVMFlang|Clang")
            list(APPEND Y12MFF_COMPILE_OPTIONS -ffp-contract=off)
        elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
            if(WIN32)
                list(APPEND Y12MFF_COMPILE_OPTIONS /Qfma-)
            else()
                list(APPEND Y12MFF_COMPILE_OPTIONS -no-fma)
            endif()
        elseif(CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC|PGI")
            list(APPEND Y12MFF_COMPILE_OPTIONS -Mnofma)
        endif()
    endif()
else()
    list(APPEND Y12MFF_COMPILE_DEFINITIONS Y12M_ACCUM_QUAD)
    message(STATUS "y12mff residual accumulator: quad precision (selected_real_kind(33))")
endif()
