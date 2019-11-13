 
# @author  Barthélémy von Haller


# Compile C++ tests (here only one, but it will grab whatever cxx is in the test directory)
file(GLOB_RECURSE TEST_FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)


set(MAIN_LIBRARY wannierlib)

# Iterate over all tests found. For each, declare an executable and add it to the tests list.
set(ALL_TESTS)
foreach(FILE ${TEST_FILES})
    string(REGEX REPLACE "[.cpp/]" "" NAME ${FILE}) # create exe name
    message( "FILE = ${NAME}")
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/test/${NAME}") #create a directory for the test
    add_executable(${NAME} ${FILE})
    target_link_libraries(${NAME} ${MAIN_LIBRARY} ) # link against our lib and boost test
    add_test(NAME ${NAME} COMMAND ${NAME}  )                              # this is how to add tests to CMake
endforeach(FILE ${TEST_FILES})

include(CTest)
