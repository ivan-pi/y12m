# CMake generated Testfile for 
# Source directory: /home/runner/work/y12m/y12m/examples
# Build directory: /home/runner/work/y12m/y12m/build/examples
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(test_matrix_generators "/home/runner/work/y12m/y12m/build/examples/test_matrix_generators")
set_tests_properties(test_matrix_generators PROPERTIES  _BACKTRACE_TRIPLES "/home/runner/work/y12m/y12m/examples/CMakeLists.txt;44;add_test;/home/runner/work/y12m/y12m/examples/CMakeLists.txt;0;")
add_test(maind "/bin/sh" "-c" "/home/runner/work/y12m/y12m/build/examples/maind < input.txt")
set_tests_properties(maind PROPERTIES  WORKING_DIRECTORY "/home/runner/work/y12m/y12m/build/examples" _BACKTRACE_TRIPLES "/home/runner/work/y12m/y12m/examples/CMakeLists.txt;50;add_test;/home/runner/work/y12m/y12m/examples/CMakeLists.txt;0;")
