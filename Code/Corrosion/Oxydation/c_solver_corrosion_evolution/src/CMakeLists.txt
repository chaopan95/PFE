cmake_minimum_required(VERSION 2.8)

set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ../bin)

set (C_SOURCE_FILES
    c_solver_corrosion_evolution_functions.c               equadiff.c
    c_solver_corrosion_evolution_optimisation_functions.c  main.c
    c_solver_corrosion_evolution_parameters.c)

IF(NOT DEFINED INIPARSER_INCLUDE_DIR OR NOT DEFINED INIPARSER_LIB_DIR)
  SET(INIPARSER_INCLUDE_DIR $ENV{INIPARSER_INSTALL_DIR}/include)
  SET(INIPARSER_LIB_DIR $ENV{INIPARSER_INSTALL_DIR}/lib)
ENDIF()

IF(NOT DEFINED MAPUTILS_INCLUDE_DIR OR NOT DEFINED MAPUTILS_LIB_DIR)
  SET(MAPUTILS_INCLUDE_DIR $ENV{MAPUTILS_INSTALL_DIR}/include)
  SET(MAPUTILS_LIB_DIR $ENV{MAPUTILS_INSTALL_DIR}/lib)
ENDIF()

link_directories(
  ${INIPARSER_LIB_DIR}
  ${MAPUTILS_LIB_DIR}
)

include_directories(
  ${INIPARSER_INCLUDE_DIR}
  ${MAPUTILS_INCLUDE_DIR}
)

add_executable(${PROJECT_NAME} ${C_SOURCE_FILES})
add_dependencies(${PROJECT_NAME} libiniparser libmaputils)

target_link_libraries(${PROJECT_NAME}
  m iniparser argtable2  maputils)

install(TARGETS ${PROJECT_NAME}
  DESTINATION bin)
