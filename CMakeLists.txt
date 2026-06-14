cmake_minimum_required(VERSION 3.14)

project(wannier2sparse
        VERSION 1.0.0
        DESCRIPTION "Expand a Wannier90 tight-binding model into a supercell and export sparse (CSR) operators"
        LANGUAGES CXX)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

option(W2SP_BUILD_TESTS "Build the wannier2sparse test suite" ON)

# Keep the executable at the top of the build tree (where it has always been).
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")

# ---------------------------------------------------------------------------
# Eigen3 dependency
#
# Preferred: a real Eigen3Config.cmake (gives the Eigen3::Eigen imported target,
# which propagates downstream via find_dependency).  Fallback: a plain include
# directory (-DEIGEN_DIR=... or the usual system path) with a warning that an
# installed wannier2sparse package will not be able to re-find Eigen.
# ---------------------------------------------------------------------------
find_package(Eigen3 3.3 NO_MODULE QUIET)
if(TARGET Eigen3::Eigen)
    set(W2SP_EIGEN_TARGET Eigen3::Eigen)
    set(W2SP_EIGEN_FROM_CONFIG TRUE)
    message(STATUS "Found Eigen3 via config: ${Eigen3_DIR}")
else()
    set(W2SP_EIGEN_FROM_CONFIG FALSE)
    find_path(EIGEN_DIR
        NAMES signature_of_eigen3_matrix_library
        HINTS /usr/include/eigen3 /usr/local/include/eigen3
        DOC "Path to the Eigen3 headers (the directory containing the Eigen/ folder)")
    if(NOT EIGEN_DIR)
        message(FATAL_ERROR "Eigen3 not found. Install libeigen3-dev or pass -DEIGEN_DIR=/path/to/eigen3")
    endif()
    add_library(w2sp_eigen INTERFACE)
    target_include_directories(w2sp_eigen INTERFACE $<BUILD_INTERFACE:${EIGEN_DIR}>)
    install(TARGETS w2sp_eigen EXPORT wannier2sparseTargets)
    set(W2SP_EIGEN_TARGET w2sp_eigen)
    message(WARNING
        "Eigen3Config.cmake not found; using the headers at '${EIGEN_DIR}' for this "
        "build only. An installed wannier2sparse package will not be able to "
        "find_dependency(Eigen3); consumers must provide Eigen themselves.")
endif()

enable_testing()

add_subdirectory(src)
if(W2SP_BUILD_TESTS)
    add_subdirectory(test)
endif()

# ---------------------------------------------------------------------------
# CMake config-package so downstream projects can do
#   find_package(wannier2sparse CONFIG REQUIRED)
#   target_link_libraries(app PRIVATE wannier2sparse::wannierlib)
# ---------------------------------------------------------------------------
install(EXPORT wannier2sparseTargets
        FILE wannier2sparseTargets.cmake
        NAMESPACE wannier2sparse::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/wannier2sparse)

configure_package_config_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/wannier2sparseConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/wannier2sparseConfig.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/wannier2sparse)

write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/wannier2sparseConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion)

install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/wannier2sparseConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/wannier2sparseConfigVersion.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/wannier2sparse)
