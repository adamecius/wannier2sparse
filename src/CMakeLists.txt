# @author  Jose H. Garcia

set(LIBRARY_NAME wannierlib)

set(PUBLIC_HEADERS
    ${PROJECT_SOURCE_DIR}/include/wannier_parser.hpp
    ${PROJECT_SOURCE_DIR}/include/hopping_list.hpp
    ${PROJECT_SOURCE_DIR}/include/sparse_matrix.hpp
    ${PROJECT_SOURCE_DIR}/include/operator_algebra.hpp
    ${PROJECT_SOURCE_DIR}/include/tbmodel.hpp
    ${PROJECT_SOURCE_DIR}/include/w2sp_arguments.hpp
)

add_library(${LIBRARY_NAME} STATIC
    wannier_parser.cpp
    hopping_list.cpp
)
add_library(wannier2sparse::wannierlib ALIAS ${LIBRARY_NAME})

target_include_directories(${LIBRARY_NAME} PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_compile_features(${LIBRARY_NAME} PUBLIC cxx_std_11)
target_link_libraries(${LIBRARY_NAME} PUBLIC ${W2SP_EIGEN_TARGET})

set_target_properties(${LIBRARY_NAME} PROPERTIES VERSION ${PROJECT_VERSION})

# The command-line application.
add_executable(wannier2sparse main.cpp)
target_link_libraries(wannier2sparse PRIVATE ${LIBRARY_NAME})

# Install the library, the application and the public headers, registering the
# library in the export set consumed by the config-package.
install(TARGETS ${LIBRARY_NAME} wannier2sparse
    EXPORT wannier2sparseTargets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
install(FILES ${PUBLIC_HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
