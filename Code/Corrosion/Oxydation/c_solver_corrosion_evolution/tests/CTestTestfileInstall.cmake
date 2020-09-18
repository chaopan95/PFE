# Post install testing script for ctest
# The conventional filename for such a script is CTestTestfile.cmake
# The funny actual file name avoids conflicts when building in source


add_test(test_c_solver_corrosion_evolution python non_regression_tests.py)
add_test(test_map_c_solver_corrosion_evolution
  python non_regression_tests.py TestMapengine)
set_tests_properties(TEST test_map_c_solver_corrosion_evolution
  PROPERTIES LABELS "mapengine"
  DEPENDS "test_c_solver_corrosion_evolution")
