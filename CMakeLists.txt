cmake_minimum_required(VERSION 3.14)
project(concord VERSION 0.1.0 LANGUAGES CXX)

# 1) Your header‑only target
add_library(concord INTERFACE)
target_include_directories(concord INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)

# 2) Install headers
include(GNUInstallDirs)
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# 3) Export for install
install(TARGETS concord EXPORT concordTargets)
install(EXPORT concordTargets
  FILE      concordTargets.cmake
  NAMESPACE concord::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/concord
)

# 4) Export for build‑tree / FetchContent
export(EXPORT concordTargets
  FILE      "${CMAKE_CURRENT_BINARY_DIR}/concordTargets.cmake"
  NAMESPACE concord::
)

# 5) Generate package‑config files
include(CMakePackageConfigHelpers)
set(CONFIG_INSTALL_DIR ${CMAKE_INSTALL_LIBDIR}/cmake/concord)

write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/concordConfigVersion.cmake"
  VERSION     ${PROJECT_VERSION}
  COMPATIBILITY AnyNewerVersion
)

# You'll need to add a file at cmake/concordConfig.cmake.in (see below)
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/concordConfig.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/concordConfig.cmake"
  @ONLY
)

# 6) Install the package‑config files
install(FILES
  "${CMAKE_CURRENT_BINARY_DIR}/concordConfig.cmake"
  "${CMAKE_CURRENT_BINARY_DIR}/concordConfigVersion.cmake"
  DESTINATION ${CONFIG_INSTALL_DIR}
)

