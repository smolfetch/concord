cmake_minimum_required(VERSION 3.14)
project(concord VERSION 0.1.0 LANGUAGES CXX)

# -----------------------------------------------------------------------------
# 1) Create the header‐only INTERFACE library
# -----------------------------------------------------------------------------
add_library(concord INTERFACE)

# Create an alias so users can write `concord::concord`
add_library(concord::concord ALIAS concord)

# -----------------------------------------------------------------------------
# 2) Public include path
# -----------------------------------------------------------------------------
# This makes <concord/enu_converter.h> etc. available
target_include_directories(concord
  INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

# -----------------------------------------------------------------------------
# 4) Installation and export
# -----------------------------------------------------------------------------
include(GNUInstallDirs)

# Install headers
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# Install and export the INTERFACE target
install(TARGETS concord EXPORT concordTargets)

# Export file for find_package()
install(EXPORT concordTargets
  FILE concordTargets.cmake
  NAMESPACE concord::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/concord
)

# Allow FetchContent users to import without installing
export(EXPORT concordTargets
  FILE "${CMAKE_CURRENT_BINARY_DIR}/concordTargets.cmake"
  NAMESPACE concord::
)
