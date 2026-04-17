
set(TEST_COMPILE_FLAGS -Wall -pedantic)

foreach(tname
    test1
    test2)

  add_executable(${tname} ${tname}.f)
  target_compile_options(${tname} PRIVATE ${TEST_COMPILE_FLAGS})
  target_link_libraries(${tname} PRIVATE y12m_legacy)
#  add_test(NAME ${tname} COMMAND ${tname})

endforeach()


add_library(y12m_mgen OBJECT matrd1.f matre1.f matrf2.f time.f90)

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

