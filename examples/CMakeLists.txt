
if(FC_IS_LFORTRAN)
    set(TEST_COMPILE_FLAGS --fixed-form --implicit-typing --implicit-interface)
else()
    set(TEST_COMPILE_FLAGS -Wall -pedantic)
endif()


# Step 1 of the examples refactoring: matrix generators and the timing helper
# are now compiled as a proper Fortran module (y12m_example_util.f90) instead
# of the four separate fixed-form files (matrd1.f, matre1.f, matrf2.f, time.f90).
add_library(y12m_mgen OBJECT y12m_example_util.f90)
target_link_libraries(y12m_mgen PUBLIC y12m_legacy)

set_target_properties(y12m_mgen PROPERTIES
    Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/include")
target_include_directories(y12m_mgen INTERFACE
    "${CMAKE_CURRENT_BINARY_DIR}/include")

if (NOT FC_IS_LFORTRAN)

foreach(tname
    test1
    test2)

  add_executable(${tname} ${tname}.f)
  target_compile_options(${tname} PRIVATE ${TEST_COMPILE_FLAGS})
  target_link_libraries(${tname} PRIVATE y12m_legacy)
#  add_test(NAME ${tname} COMMAND ${tname})

endforeach()

foreach(tname
    maine
    maind
    mainf
    mainf2)

  add_executable(${tname} ${tname}.f)
  target_compile_options(${tname} PRIVATE ${TEST_COMPILE_FLAGS})
  target_link_libraries(${tname} PRIVATE y12m_legacy y12m_mgen)
#  add_test(NAME ${tname} COMMAND ${tname})

endforeach()

# 1. Create the file as you did before
file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/input.txt" "1 600 600 600 4 44 204\n")

# 2. Call the shell directly to handle the redirection
add_test(NAME maind
    COMMAND /bin/sh -c "$<TARGET_FILE:maind> < input.txt"
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)

endif() # FC_IS_LFORTRAN

# Unit tests for the matrix generators
add_executable(test_matrix_generators test_matrix_generators.f90)
target_compile_options(test_matrix_generators PRIVATE )
target_link_libraries(test_matrix_generators PRIVATE y12m_legacy y12m_mgen)
add_test(NAME test_matrix_generators COMMAND test_matrix_generators)

