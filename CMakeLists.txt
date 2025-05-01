cmake_minimum_required(VERSION 3.14)
project(concord VERSION 0.2.0 LANGUAGES CXX)

# -----------------------------------------------------------------------------
# 1) Create the header‐only INTERFACE library
# -----------------------------------------------------------------------------
add_library(concord INTERFACE)

# Allow users to link via `concord::concord`
add_library(concord::concord ALIAS concord)

# -----------------------------------------------------------------------------
# 2) Public include path
# -----------------------------------------------------------------------------
# This makes <concord/enu_converter.h> etc. resolve to include/concord/...
target_include_directories(concord
  INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

# -----------------------------------------------------------------------------
# 3) Installation & Export for find_package()
# -----------------------------------------------------------------------------
include(GNUInstallDirs)

# 3a) Install headers
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# 3b) Install and export the INTERFACE target
install(TARGETS concord EXPORT concordTargets)

# 3c) Install the export file so that `find_package(concord)` works
install(EXPORT concordTargets
  FILE concordTargets.cmake
  NAMESPACE concord::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/concord
)

# 3d) Also export it for in‐source FetchContent usage
export(EXPORT concordTargets
  FILE "${CMAKE_CURRENT_BINARY_DIR}/concordTargets.cmake"
  NAMESPACE concord::
)

