 
# @author  Barthélémy von Haller


# Compile C++ tests (grabs every .cpp in the test directory)
file(GLOB_RECURSE TEST_FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)
# The WannierBerri independent cross-check is off by default (needs the package +
# a generated reference). Enable with -DW2SP_WBERRI_CHECK=ON (label `wberri`).
list(REMOVE_ITEM TEST_FILES wberri_crosscheck.cpp)

set( TEST_ROOT ${CMAKE_BINARY_DIR}/test)
configure_file(graphene_hr.dat "${TEST_ROOT}/model_from_files/graphene_hr.dat" COPYONLY)

configure_file(spin_graphene_hr.dat "${TEST_ROOT}/density_and_current_from_files/spin_graphene_hr.dat" COPYONLY)
configure_file(spin_graphene.uc "${TEST_ROOT}/density_and_current_from_files/spin_graphene.uc" COPYONLY)
configure_file(spin_graphene.xyz "${TEST_ROOT}/density_and_current_from_files/spin_graphene.xyz" COPYONLY)
configure_file(spin_graphene_hr.dat "${TEST_ROOT}/csr_count_equivalence/spin_graphene_hr.dat" COPYONLY)
configure_file(spin_graphene.uc "${TEST_ROOT}/csr_count_equivalence/spin_graphene.uc" COPYONLY)
configure_file(spin_graphene.xyz "${TEST_ROOT}/csr_count_equivalence/spin_graphene.xyz" COPYONLY)


set(MAIN_LIBRARY wannierlib)

# Iterate over all tests found. For each, declare an executable and add it to the tests list.
set(ALL_TESTS)
foreach(FILE ${TEST_FILES})
    string(REGEX REPLACE ".cpp" "" NAME ${FILE}) # create exe name
    set( TEST_DIR "${TEST_ROOT}/${NAME}")
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${TEST_DIR}" ) #create a directory for the test
    add_executable(${NAME} ${FILE})
    target_link_libraries(${NAME} ${MAIN_LIBRARY} ) # link against our lib and boost test
    add_test(NAME ${NAME} COMMAND ${NAME} WORKING_DIRECTORY "${TEST_DIR}"  )                              # this is how to add tests to CMake
endforeach(FILE ${TEST_FILES})

# Optional WannierBerri independent cross-check (Plan 10C).
if(W2SP_WBERRI_CHECK)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${TEST_ROOT}/wberri_crosscheck")
    add_executable(wberri_crosscheck wberri_crosscheck.cpp)
    target_link_libraries(wberri_crosscheck ${MAIN_LIBRARY})
    add_test(NAME wberri_crosscheck COMMAND wberri_crosscheck WORKING_DIRECTORY "${TEST_ROOT}/wberri_crosscheck")
    set_tests_properties(wberri_crosscheck PROPERTIES LABELS wberri)
endif()

include(CTest)
