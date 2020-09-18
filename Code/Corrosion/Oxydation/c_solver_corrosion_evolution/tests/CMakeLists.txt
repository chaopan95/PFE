# The project name is the name of the map component
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

# Copy source in build dir
file(COPY . DESTINATION .
  PATTERN ".svn" EXCLUDE)

enable_testing()

add_test(test_${PROJECT_NAME} ${PYTHON_EXECUTABLE} non_regression_tests.py)

add_test(test_map_${PROJECT_NAME}
   ${PYTHON_EXECUTABLE} non_regression_tests.py TestMapengine)

set_property(TEST test_map_${PROJECT_NAME}
  PROPERTY LABELS "mapengine")

# install all tests
install(DIRECTORY .
  DESTINATION share/map/tests/components/${PROJECT_NAME}
  PATTERN ".svn" EXCLUDE
  PATTERN "CMakeLists.txt" EXCLUDE
  PATTERN "*.cmake" EXCLUDE
  )

install(FILES CTestTestfileInstall.cmake
  DESTINATION share/map/tests/components/${PROJECT_NAME}
  RENAME CTestTestfile.cmake)