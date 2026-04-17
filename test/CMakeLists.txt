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
    test_y12mg_mh)

  add_executable(${tname} ${tname}.f90)
  target_compile_options(${tname} PRIVATE ${TEST_COMPILE_FLAGS})
  target_link_libraries(${tname} PRIVATE y12m_legacy)
  add_test(NAME ${tname} COMMAND ${tname})

endforeach()
