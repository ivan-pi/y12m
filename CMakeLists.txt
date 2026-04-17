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

# Build a static library from the legacy sources.
# These are pre-existing files; no extra flags are applied to them.
add_library(y12m_legacy STATIC ${LEGACY_SOURCES})

enable_testing()
add_subdirectory(test)
