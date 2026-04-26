# SPDX-License-Identifier: GPL-2.0-only
# Test programs for the y12m sparse linear algebra package.

#set(TEST_COMPILE_FLAGS -Wall -pedantic)

foreach(tname
    test_y12ma_sp
    test_y12ma_dp
    test_y12mb_mc_md_sp
    test_y12mb_mc_md_dp
    test_y12mf_sp
    test_y12mg_mh
    test_y12mc_z_intent
    test_y12ma_5x5)

  add_executable(${tname} ${tname}.f90)
  target_compile_options(${tname} PRIVATE ${TEST_COMPILE_FLAGS})
  target_link_libraries(${tname} PRIVATE y12m_legacy)
  add_test(NAME ${tname} COMMAND ${tname})

endforeach()

# ---- File-based test driver ----
add_executable(y12m_solve y12m_solve.f90)
target_compile_options(y12m_solve PRIVATE ${TEST_COMPILE_FLAGS})
target_link_libraries(y12m_solve PRIVATE y12m_legacy)

# Drive the two bundled input matrices.  The data files live in the source
# tree; pass absolute paths so the tests work regardless of the build directory.
set(DATA_DIR "${CMAKE_CURRENT_SOURCE_DIR}/data")

foreach(MAT_ID mat0 mat1 mat2 mat3 mat5)
    add_test(
        NAME "y12m_solve_${MAT_ID}"
        COMMAND y12m_solve --verbose "${DATA_DIR}/${MAT_ID}.mtx"
    )
endforeach()

# ---- Random 9x9 double-precision test (requires LAPACK) ----
if(LAPACK_FOUND)
  add_executable(test_y12ma_rand_dp test_y12ma_rand_dp.f90)
  target_compile_options(test_y12ma_rand_dp PRIVATE ${TEST_COMPILE_FLAGS})
  target_link_libraries(test_y12ma_rand_dp PRIVATE y12m_legacy LAPACK::LAPACK)
  add_test(NAME test_y12ma_rand_dp COMMAND test_y12ma_rand_dp)
else()
  message(STATUS
    "LAPACK not found -- skipping test_y12ma_rand_dp (sparse-dense comparison test)")
endif()


# 1. Define the Opt-in Option
option(ENABLE_MMIO_TEST "Build test driver with matrix-market input" OFF)

if(ENABLE_MMIO_TEST)
    set(MMIO_FILE "${CMAKE_CURRENT_BINARY_DIR}/mmio.f")

    # 2. Download the file at configure time
    message(STATUS "Downloading mmio.f from NIST...")
    file(DOWNLOAD
        https://math.nist.gov/MatrixMarket/mmio/f/mmio.f
        "${MMIO_FILE}"
        EXPECTED_HASH SHA256=f3146809daffd47587ecd0fa8af0187354a34bb41d3740cd71aab33e505d99a1
        SHOW_PROGRESS
        STATUS download_status
    )

    list(GET download_status 0 _status_code)
    list(GET download_status 1 _status_message)

    if(_status_code EQUAL 0)
        message(STATUS "mmio.f downloaded successfully.")

        add_library(mmio_lib STATIC "${MMIO_FILE}")
        #add_executable(mm_test_driver main.f90)
        #target_link_libraries(mm_test_driver PRIVATE mmio_lib)
        #target_compile_definitions(mm_test_driver PRIVATE HAS_MMIO)
    else()
        # Non-fatal warning: the rest of the build continues
        message(WARNING "MMIO download failed (${_status_message}). Skipping mm_test_driver.")
    endif()

endif()
