cmake_minimum_required(VERSION 3.15)

set(project_name concord)
project(${project_name} VERSION 0.5.5 LANGUAGES CXX)

add_compile_options(-Wall -Wextra -Wpedantic)
cmake_policy(SET CMP0074 NEW)
cmake_policy(SET CMP0135 NEW)
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED On)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(params
    -Wno-reorder
)

option(BUILD_EXAMPLES "Build examples" OFF)

# --------------------------------------------------------------------------------------------------
add_library(${project_name} INTERFACE)
add_library(concord::concord ALIAS concord)
target_include_directories(${project_name} INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

# --------------------------------------------------------------------------------------------------
include(GNUInstallDirs)

#Install headers
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

#Install and export the INTERFACE target
install(TARGETS ${project_name} EXPORT ${project_name}Targets)

#Install the export file so that `find_package(${project_name})` works
install(EXPORT ${project_name}Targets
  FILE ${project_name}Targets.cmake
  NAMESPACE ${project_name}::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${project_name}
)

#Also export it for in‐source FetchContent usage
export(EXPORT ${project_name}Targets
  FILE "${CMAKE_CURRENT_BINARY_DIR}/${project_name}Targets.cmake"
  NAMESPACE ${project_name}::
)


# --------------------------------------------------------------------------------------------------
if(BUILD_EXAMPLES)
  set(external_deps
  )
  set(exe
    examples/main.cpp
  )
  # ----------------------------------------------
  set(exec_names)
  foreach(src_file IN LISTS exe)
    get_filename_component(exec_name "${src_file}" NAME_WE)
    add_executable(${exec_name} "${src_file}")
      foreach(lib_file IN LISTS internal_deps)
        target_compile_options(${exec_name} PRIVATE ${params})
        target_sources(${exec_name} PRIVATE "${lib_file}")
      endforeach()
    target_link_libraries(${exec_name} ${external_deps})
    install(TARGETS ${exec_name} DESTINATION bin)
    list(APPEND exec_names ${exec_name})
  endforeach()
  # ----------------------------------------------
  foreach(exec IN LISTS exec_names)
    file(REMOVE "${CMAKE_CURRENT_LIST_DIR}/.execs")
    file(WRITE "${CMAKE_CURRENT_LIST_DIR}/.execs")
    get_filename_component(exec_name "${exec}" NAME)
    file(APPEND "${CMAKE_CURRENT_LIST_DIR}/.execs" "${CMAKE_CURRENT_LIST_DIR}/build/${exec_name}\n")
  endforeach()
  # ----------------------------------------------
  install(TARGETS
    ${exec_names}
    DESTINATION lib/${PROJECT_NAME}
  )
endif()
