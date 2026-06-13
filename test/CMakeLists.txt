 
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

# Plan 11 -- regenerate the WannierBerri committed goldens (manual; needs the
# package + a fixture with <seed>.chk). NOT part of `ctest`; run explicitly:
#   W2SP_WS_CONVENTION=<conv> cmake --build build --target regen-wberri-golden
# Fixture prefixes default to the local /tmp paths; override on the cmake line.
set(W2SP_SPIN_FIXTURE   "/tmp/fix/Fe/Fe"        CACHE STRING "spin fixture <prefix> for the WannierBerri golden")
set(W2SP_ORBITAL_FIXTURE "/tmp/fix/copper/copper" CACHE STRING "orbital fixture <prefix> for the WannierBerri golden")
add_custom_target(regen-wberri-golden
    COMMENT "Regenerating WannierBerri goldens (needs: pip install wannierberri + fixtures with .chk + W2SP_WS_CONVENTION)"
    COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_SOURCE_DIR}/golden
    COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/fixtures/wberri_reference.py ${W2SP_SPIN_FIXTURE}    spin    ${CMAKE_CURRENT_SOURCE_DIR}/golden
    COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/fixtures/wberri_reference.py ${W2SP_ORBITAL_FIXTURE} orbital ${CMAKE_CURRENT_SOURCE_DIR}/golden
    VERBATIM)

include(CTest)
