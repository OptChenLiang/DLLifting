cmake_minimum_required(VERSION 3.14)
project(dllifting VERSION 1.1.0 LANGUAGES CXX C)

option(DLLIFTING_BUILD_SHARED "Build shared library" OFF)
option(DLLIFTING_BUILD_EXAMPLES "Build example programs" ON)
option(DLLIFTING_BUILD_TESTS "Build unit tests" ON)
option(DLLIFTING_REDUCTION "Enable capacity reduction during lifting" ON)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

set(DLLIFTING_SOURCES src/DLLifting.cpp src/dllifting_c.cpp)
set(DLLIFTING_HEADERS include/DLLifting.h include/dllifting_c.h)

if(DLLIFTING_BUILD_SHARED)
  add_library(dllifting SHARED ${DLLIFTING_SOURCES})
  target_compile_definitions(dllifting PRIVATE DLLIFTING_EXPORTS)
else()
  add_library(dllifting STATIC ${DLLIFTING_SOURCES})
endif()

target_include_directories(dllifting
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

if(DLLIFTING_REDUCTION)
  target_compile_definitions(dllifting PUBLIC DLLIFTING_REDUCTION)
endif()

target_compile_options(dllifting PRIVATE -Wall -Wextra)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(dllifting PRIVATE -O2)
endif()

find_library(M_LIB m)
if(M_LIB)
  target_link_libraries(dllifting PUBLIC ${M_LIB})
endif()

include(GNUInstallDirs)

install(TARGETS dllifting
  EXPORT dlliftingTargets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

install(FILES ${DLLIFTING_HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(EXPORT dlliftingTargets
  FILE dlliftingTargets.cmake
  NAMESPACE dllifting::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/dllifting
)

include(CMakePackageConfigHelpers)
configure_package_config_file(
  cmake/dlliftingConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/dlliftingConfig.cmake
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/dllifting
)
write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/dlliftingConfigVersion.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion
)
install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/dlliftingConfig.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/dlliftingConfigVersion.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/dllifting
)

if(DLLIFTING_BUILD_EXAMPLES)
  add_subdirectory(examples)
endif()

if(DLLIFTING_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()
