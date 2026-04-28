# CMake generated Testfile for 
# Source directory: /home/runner/work/y12m/y12m/examples
# Build directory: /home/runner/work/y12m/y12m/build2/examples
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(test_matrix_generators "/home/runner/work/y12m/y12m/build2/examples/test_matrix_generators")
set_tests_properties(test_matrix_generators PROPERTIES  _BACKTRACE_TRIPLES "/home/runner/work/y12m/y12m/examples/CMakeLists.txt;44;add_test;/home/runner/work/y12m/y12m/examples/CMakeLists.txt;0;")
add_test(maind "/bin/sh" "-c" "/home/runner/work/y12m/y12m/build2/examples/maind < input.txt")
set_tests_properties(maind PROPERTIES  WORKING_DIRECTORY "/home/runner/work/y12m/y12m/build2/examples" _BACKTRACE_TRIPLES "/home/runner/work/y12m/y12m/examples/CMakeLists.txt;50;add_test;/home/runner/work/y12m/y12m/examples/CMakeLists.txt;0;")
add_test(timings_de_d "/home/runner/work/y12m/y12m/build2/examples/timings_de" "-M" "d" "-n" "100" "-c" "10")
set_tests_properties(timings_de_d PROPERTIES  WORKING_DIRECTORY "/home/runner/work/y12m/y12m/build2/examples" _BACKTRACE_TRIPLES "/home/runner/work/y12m/y12m/examples/CMakeLists.txt;61;add_test;/home/runner/work/y12m/y12m/examples/CMakeLists.txt;0;")
add_test(timings_de_e "/home/runner/work/y12m/y12m/build2/examples/timings_de" "-M" "e" "-n" "100" "-c" "10")
set_tests_properties(timings_de_e PROPERTIES  WORKING_DIRECTORY "/home/runner/work/y12m/y12m/build2/examples" _BACKTRACE_TRIPLES "/home/runner/work/y12m/y12m/examples/CMakeLists.txt;66;add_test;/home/runner/work/y12m/y12m/examples/CMakeLists.txt;0;")
add_test(matf_bench "/home/runner/work/y12m/y12m/build2/examples/matf_bench" "-n" "50" "-c" "15" "-r" "3")
set_tests_properties(matf_bench PROPERTIES  WORKING_DIRECTORY "/home/runner/work/y12m/y12m/build2/examples" _BACKTRACE_TRIPLES "/home/runner/work/y12m/y12m/examples/CMakeLists.txt;77;add_test;/home/runner/work/y12m/y12m/examples/CMakeLists.txt;0;")
add_test(poisson_9pt "/home/runner/work/y12m/y12m/build2/examples/poisson_9pt")
set_tests_properties(poisson_9pt PROPERTIES  WORKING_DIRECTORY "/home/runner/work/y12m/y12m/build2/examples" _BACKTRACE_TRIPLES "/home/runner/work/y12m/y12m/examples/CMakeLists.txt;87;add_test;/home/runner/work/y12m/y12m/examples/CMakeLists.txt;0;")
