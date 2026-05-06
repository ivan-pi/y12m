
set(TEST_COMPILE_FLAGS -Wall -pedantic)

foreach(tname
    test1
    test2)

  add_executable(${tname} ${tname}.f)
  target_compile_options(${tname} PRIVATE ${TEST_COMPILE_FLAGS})
  target_link_libraries(${tname} PRIVATE y12m_legacy)
#  add_test(NAME ${tname} COMMAND ${tname})

endforeach()


# Step 1 of the examples refactoring: matrix generators and the timing helper
# are now compiled as a proper Fortran module (y12m_example_util.f90) instead
# of the four separate fixed-form files (matrd1.f, matre1.f, matrf2.f, time.f90).
add_library(y12m_mgen OBJECT y12m_example_util.f90)
target_link_libraries(y12m_mgen PUBLIC y12m_legacy)

set_target_properties(y12m_mgen PROPERTIES
    Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/include")
target_include_directories(y12m_mgen INTERFACE
    "${CMAKE_CURRENT_BINARY_DIR}/include")

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

# Unit tests for the matrix generators
add_executable(test_matrix_generators test_matrix_generators.f90)
target_compile_options(test_matrix_generators PRIVATE ${TEST_COMPILE_FLAGS})
target_link_libraries(test_matrix_generators PRIVATE y12m_legacy y12m_mgen)
add_test(NAME test_matrix_generators COMMAND test_matrix_generators)

# 1. Create the file as you did before
file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/input.txt" "1 600 600 600 4 44 204\n")

# 2. Call the shell directly to handle the redirection
add_test(NAME maind
    COMMAND /bin/sh -c "$<TARGET_FILE:maind> < input.txt"
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)

# Shared helpers for the modern example programs
add_library(y12m_example_helpers OBJECT y12m_example_helpers.f90)
set_target_properties(y12m_example_helpers PROPERTIES
    Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/include")
target_include_directories(y12m_example_helpers INTERFACE
    "${CMAKE_CURRENT_BINARY_DIR}/include")

# Double-precision timing driver for class D and E matrices
add_executable(timings_de timings_de.f90)
target_compile_options(timings_de PRIVATE ${TEST_COMPILE_FLAGS})
target_link_libraries(timings_de PRIVATE y12m_legacy y12m_example_helpers)

# Smoke tests: one for class D and one for class E
add_test(NAME timings_de_d
    COMMAND timings_de -M d -n 100 -c 10
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)

add_test(NAME timings_de_e
    COMMAND timings_de -M e -n 100 -c 10
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)

# Double-precision timing driver for class F2 matrices
add_executable(matf_bench matf_bench.f90)
target_compile_options(matf_bench PRIVATE ${TEST_COMPILE_FLAGS})
target_link_libraries(matf_bench PRIVATE y12m_legacy y12m_example_helpers)

# Smoke test: small F2 matrix (n=50, c=15, r=3)
add_test(NAME matf_bench
    COMMAND matf_bench -n 50 -c 15 -r 3
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)

# Simple function: Name of the executable and the exact source file
function(add_pde_test NAME SOURCE)
    add_executable(${NAME} ${SOURCE})
    target_compile_options(${NAME} PRIVATE ${TEST_COMPILE_FLAGS})
    target_link_libraries(${NAME} PRIVATE y12m_legacy)

    add_test(NAME ${NAME}
             COMMAND ${NAME}
             WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
endfunction()


add_pde_test(poisson_2d poisson_2d.F90)
add_pde_test(poisson_3d poisson_3d.f90)
add_pde_test(biharmonic_13pt biharmonic_13pt.f90)
target_link_libraries(biharmonic_13pt PRIVATE y12m_example_helpers)
add_pde_test(conv_diff_cds conv_diff_cds.f90)
target_link_libraries(conv_diff_cds PRIVATE y12m_example_helpers)
add_pde_test(fem_anisotropic fem_anisotropic.f90)
target_link_libraries(fem_anisotropic PRIVATE y12m_example_helpers)
add_pde_test(darcy_flow darcy_flow.f90)
add_pde_test(elasticity elasticity.f90)
add_pde_test(euler_bernoulli_beam euler_bernoulli_beam.f90)

# --------------------------------------------------------------------------
# LU-reuse examples: demonstrate the low-level y12mb / y12mc / y12md API
# --------------------------------------------------------------------------

# Case (ii): implicit heat equation -- factorize once, solve many times
add_pde_test(heat_implicit heat_implicit.f90)

# Case (ii): Poisson with multiple load cases -- one factorization, nrhs solves
add_pde_test(multiple_loads multiple_loads.f90)

# Case (iii): Newton iteration for Bratu -- reuse Markowitz ordering
add_pde_test(newton_bratu newton_bratu.f90)