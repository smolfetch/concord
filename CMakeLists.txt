cmake_minimum_required(VERSION 3.14)
project(concord VERSION 0.1.0 LANGUAGES CXX)

# Find GeographicLib (provides LocalCartesian & UTMUPS)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "/usr/share/cmake/geographiclib/")
find_package(GeographicLib REQUIRED)

# Header‑only library
add_library(concord INTERFACE)

target_include_directories(concord INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)

# Pull in GeographicLib for consumers
target_link_libraries(concord INTERFACE GeographicLib::Geographic)

# Install headers and export config
include(GNUInstallDirs)

install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(TARGETS concord EXPORT concordTargets)
install(EXPORT concordTargets
  FILE concordTargets.cmake
  NAMESPACE concord::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/concord
)

# Allow use via FetchContent without installing
export(EXPORT concordTargets
  FILE "${CMAKE_CURRENT_BINARY_DIR}/concordTargets.cmake"
  NAMESPACE concord::
)

