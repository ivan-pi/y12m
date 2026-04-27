# SPDX-License-Identifier: GPL-2.0-only
# Test programs for the y12m sparse linear algebra package.
# Each program uses -Wall -pedantic as requested.

set(TEST_COMPILE_FLAGS -Wall -pedantic)

foreach(tname
    test_y12ma_sp
    test_y12ma_dp
    test_y12mb_mc_md_sp
    test_y12mb_mc_md_dp
    test_y12mf_sp
    test_y12mg_mh
    test_y12ma_5x5)

  add_executable(${tname} ${tname}.f90)
  target_compile_options(${tname} PRIVATE ${TEST_COMPILE_FLAGS})
  target_link_libraries(${tname} PRIVATE y12m_legacy)
  add_test(NAME ${tname} COMMAND ${tname})

endforeach()

foreach(tname
    test_y12mb_errors_sp
    test_y12mb_errors_dp)

  add_executable(${tname} test_y12mb_errors.F90)
  target_compile_options(${tname} PRIVATE ${TEST_COMPILE_FLAGS})
  target_link_libraries(${tname} PRIVATE y12m_legacy)
  # Check if the target name ends with _sp
  if(tname MATCHES "_sp$")
    target_compile_definitions(${tname} PRIVATE TEST_SINGLE_PRECISION)
  endif()
  add_test(NAME ${tname} COMMAND ${tname})

endforeach()

foreach(tname
    test_y12mc_errors_sp
    test_y12mc_errors_dp)

  add_executable(${tname} test_y12mc_errors.F90)
  target_compile_options(${tname} PRIVATE ${TEST_COMPILE_FLAGS})
  target_link_libraries(${tname} PRIVATE y12m_legacy)
  if(tname MATCHES "_sp$")
    target_compile_definitions(${tname} PRIVATE TEST_SINGLE_PRECISION)
  endif()
  add_test(NAME ${tname} COMMAND ${tname})

endforeach()

foreach(tname
    test_y12md_errors_sp
    test_y12md_errors_dp)

  add_executable(${tname} test_y12md_errors.F90)
  target_compile_options(${tname} PRIVATE ${TEST_COMPILE_FLAGS})
  target_link_libraries(${tname} PRIVATE y12m_legacy)
  if(tname MATCHES "_sp$")
    target_compile_definitions(${tname} PRIVATE TEST_SINGLE_PRECISION)
  endif()
  add_test(NAME ${tname} COMMAND ${tname})

endforeach()

foreach(tname
    test_y12mc_z_intent_sp
    test_y12mc_z_intent_dp)

  add_executable(${tname} test_y12mc_z_intent.F90)
  target_compile_options(${tname} PRIVATE ${TEST_COMPILE_FLAGS})
  target_link_libraries(${tname} PRIVATE y12m_legacy)
  if(tname MATCHES "_sp$")
    target_compile_definitions(${tname} PRIVATE TEST_SINGLE_PRECISION)
  endif()
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

# ---- Random 9x9 test (single and double precision, requires LAPACK) ----
if(LAPACK_FOUND)
  foreach(tname
      test_y12ma_rand_sp
      test_y12ma_rand_dp)

    add_executable(${tname} test_y12ma_rand.F90)
    target_compile_options(${tname} PRIVATE ${TEST_COMPILE_FLAGS})
    target_link_libraries(${tname} PRIVATE y12m_legacy LAPACK::LAPACK)
    if(tname MATCHES "_sp$")
      target_compile_definitions(${tname} PRIVATE TEST_SINGLE_PRECISION)
    endif()
    add_test(NAME ${tname} COMMAND ${tname})

  endforeach()
else()
  message(STATUS
    "LAPACK not found -- skipping test_y12ma_rand_sp/dp (sparse-dense comparison tests)")
endif()
