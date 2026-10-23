from anuga import Domain
from anuga import Dirichlet_boundary
from anuga import Operator

import numpy as num
import unittest

class Test_Operator(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_create_operator(self):
        points = num.array([[0.0,0.0],[1.0,0.0],[0.0,1.0]])

        elements = num.array([[0,1,2]])
        boundary_map = {}
        boundary_map[(0,0)] = 'edge0'
        boundary_map[(0,1)] = 'edge1'
        boundary_map[(0,2)] = 'edge2'

        domain = Domain(points, elements, boundary_map)

        operator = Operator(domain)

        message = operator.statistics()
        assert message == 'You need to implement operator statistics for your operator'

        message = operator.timestepping_statistics()
        assert message == 'You need to implement timestepping statistics for your operator'

        domain.timestep = 3.0

        assert operator.get_timestep() == domain.get_timestep()

        try:
            operator()
        except Exception:
            pass
        else:
            raise Exception('should have raised an exception')


class Test_Operator_extra(unittest.TestCase):
    """Additional tests for base_operator.py."""

    def _make_domain(self):
        points = num.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
        elements = num.array([[0, 1, 2]])
        boundary_map = {(0, 0): 'edge0', (0, 1): 'edge1', (0, 2): 'edge2'}
        from anuga import Domain
        return Domain(points, elements, boundary_map)

    def test_description_not_none(self):
        """Operator with description stores it (line 40)."""
        domain = self._make_domain()
        op = Operator(domain, description='test description')
        self.assertEqual(op.description, 'test description')

    def test_label_not_none(self):
        """Operator with label stores prefixed label (line 104)."""
        domain = self._make_domain()
        op = Operator(domain, label='myop')
        self.assertIn('myop', op.label)

    def test_print_statistics(self):
        """print_statistics calls print (line 86)."""
        domain = self._make_domain()
        op = Operator(domain)
        op.print_statistics()  # should not raise

    def test_print_timestepping_statistics(self):
        """print_timestepping_statistics calls print (line 90)."""
        domain = self._make_domain()
        op = Operator(domain)
        op.print_timestepping_statistics()  # should not raise

    def test_log_timestepping_statistics_when_logging(self):
        """log_timestepping_statistics logs when self.logging=True (line 97)."""
        import tempfile
        import os
        domain = self._make_domain()
        op = Operator(domain, logging=True)
        orig = os.getcwd()
        with tempfile.TemporaryDirectory() as tmpdir:
            os.chdir(tmpdir)
            try:
                op.activate_logging()   # creates log file
                op.log_timestepping_statistics()  # should write to log
            finally:
                os.chdir(orig)


################################################################################

if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(Test_Operator)
    runner = unittest.TextTestRunner()
    runner.run(suite)
