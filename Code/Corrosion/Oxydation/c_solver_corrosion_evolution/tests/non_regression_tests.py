#!/usr/bin/env python
"""Simple python script file for automatic non-regression tests.

Exit status is zero only if all tests pass.

"""
import os
import unittest
import filecmp

class Test_c_solver_corrosion_evolution(unittest.TestCase):
    """ Tests for for the c_solver_corrosion_evolution component."""

    component_name = 'c_solver_corrosion_evolution'
    def setUp(self):
        """Guess how to run.
        If the original component executable is nearby, this is probably
        development time, use it.
        Otherwise, it's  post install time, use environment """
        _component_name = 'c_solver_corrosion_evolution'
        _component_script = '../bin/%s' % _component_name
        if os.path.isfile(_component_script):
            self.run_script = "%s -v -f " % _component_script
        else :
            self.run_script = '%s -v -f ' % _component_name

    def run_tests(self, test_names):
        """Base method for testing."""
        output_dir = '.'
        ref_dir = "refs"
        for _testname in test_names:
            os.popen("%s %s.input" % (self.run_script,
                                      os.path.join(output_dir, _testname)))
            _ref = os.path.join(ref_dir, "%s.output" % _testname)
            _new = os.path.join(output_dir, "%s.output" % _testname)
            self.assertTrue(filecmp.cmp(_ref, _new, shallow=False),
                            "Files %s and %s differ" % (_ref, _new))
            os.remove(_new)

    def test_output_optimisation(self):
        """Test identity of optimisation output file."""
        test_names = ('c_solver_corrosion_evolution_optimisation_1',)
        output_dir = '.'
        ref_dir = "refs"
        os.popen("c_solver_corrosion_evolution")
        for _testname in test_names:
            _ref = os.path.join(ref_dir, "%s.output" % _testname)
            _new = os.path.join(output_dir, "%s.output" % _testname)
            self.assertTrue(filecmp.cmp(_ref, _new, shallow=False),
                            "Files %s and %s differ" % (_ref, _new))
            os.remove(_new)


    def test_output_2010(self):
        """Test identity of output files for Seyeux_2010 model."""
        self.run_tests(['c_solver_corrosion_evolution_test_1',
                        'c_solver_corrosion_evolution_test_2',
                        'c_solver_corrosion_evolution_test_3'])

    def test_output_2012_vO_and_dissol(self):
        """Test identity of output files for Leistner_2012 model with vO and
        dissol only."""
        self.run_tests(['c_solver_corrosion_evolution_test_4',
                        'c_solver_corrosion_evolution_test_5'])
    def test_output_2012_vO_and_dissol_and_MCr(self):
        """Test identity of output files for Leistner_2012 model with mCr."""
        self.run_tests(['c_solver_corrosion_evolution_test_6'])

    def test_output_2012_vO_and_dissol_and_MCr_and_ICr(self):
        """Test identity of output files for Leistner_2012 model with mCr and
        ICr."""
        self.run_tests(['c_solver_corrosion_evolution_test_7',])

    def test_output_2014_H(self):
        """Test identity of output files for Voyshnis_2014 model with H output."""
        self.run_tests(['c_solver_corrosion_evolution_test_8',])


class TestMapengine(Test_c_solver_corrosion_evolution):
    """Tests for the component with mapengine.

    Run the same tests as the parent class, through mapengine"""

    def setUp(self):
        """Define how to run.
        The component is run with mapengine"""
        self.run_script = 'map run -n %s -i ' % self.component_name

if __name__ == '__main__':
    unittest.main(defaultTest='Test_c_solver_corrosion_evolution')

