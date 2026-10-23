"""
Tests for anuga.validation_utilities

Covers the four public helpers exported from the package:

  run_validation_script  -- thin wrapper around anuga.run_anuga_script
  typeset_report         -- wraps pdflatex; silently no-ops when absent
  save_parameters_tex    -- writes a LaTeX verbatim block of domain parameters
  produce_report         -- orchestrates script + plot + typeset steps
"""

import os
import sys
import tempfile
import shutil
import unittest
import importlib
from unittest.mock import patch, MagicMock, call

# __init__.py imports the function `produce_report` into the package namespace,
# shadowing the submodule of the same name for getattr-based lookup. Use
# importlib.import_module to get the actual module via sys.modules.
_pr_mod = importlib.import_module('anuga.validation_utilities.produce_report')


class TestRunValidationScript(unittest.TestCase):
    """run_validation_script is a backwards-compatibility shim over run_anuga_script."""

    def test_delegates_to_run_anuga_script(self):
        with patch('anuga.run_anuga_script') as mock_run:
            from anuga.validation_utilities.run_validation import run_validation_script
            run_validation_script('myscript.py', np=2, cfl=0.9, alg='DE1',
                                  verbose=True)
            mock_run.assert_called_once_with('myscript.py', np=2, cfl=0.9,
                                             alg='DE1', verbose=True)

    def test_default_args(self):
        with patch('anuga.run_anuga_script') as mock_run:
            from anuga.validation_utilities.run_validation import run_validation_script
            run_validation_script('script.py')
            mock_run.assert_called_once_with('script.py', np=1, cfl=None,
                                             alg=None, verbose=False)

    def test_accessible_from_package(self):
        from anuga.validation_utilities import run_validation_script
        self.assertTrue(callable(run_validation_script))


class TestTypeset(unittest.TestCase):
    """typeset_report wraps pdflatex; must not raise when pdflatex is absent."""

    def test_does_not_raise_when_pdflatex_missing(self):
        """Simulate pdflatex absent by making check_output raise OSError."""
        import subprocess
        with patch('subprocess.check_output',
                   side_effect=OSError('pdflatex not found')):
            from anuga.validation_utilities.typeset_report import typeset_report
            # Should complete silently — CalledProcessError and OSError are caught
            typeset_report(report_name='dummy', verbose=False)

    def test_does_not_raise_on_called_process_error(self):
        import subprocess
        err = subprocess.CalledProcessError(1, 'pdflatex')
        with patch('subprocess.check_output', side_effect=err):
            from anuga.validation_utilities.typeset_report import typeset_report
            typeset_report(report_name='report', verbose=False)

    def test_accessible_from_package(self):
        from anuga.validation_utilities import typeset_report
        self.assertTrue(callable(typeset_report))


class TestSaveParametersTex(unittest.TestCase):
    """save_parameters_tex writes a parameters.tex file in the cwd."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        self.old_cwd = os.getcwd()
        os.chdir(self.tmpdir)

    def tearDown(self):
        os.chdir(self.old_cwd)
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def _mock_domain(self, params=None):
        if params is None:
            params = {'alg': 'DE1', 'cfl': 1.0}
        domain = MagicMock()
        domain.get_algorithm_parameters.return_value = params
        return domain

    def test_creates_parameters_tex(self):
        from anuga.validation_utilities.save_parameters_tex import save_parameters_tex
        domain = self._mock_domain()
        save_parameters_tex(domain)
        self.assertTrue(os.path.exists('parameters.tex'))

    def test_file_contains_verbatim_delimiters(self):
        from anuga.validation_utilities.save_parameters_tex import save_parameters_tex
        domain = self._mock_domain({'alg': 'DE0'})
        save_parameters_tex(domain)
        content = open('parameters.tex').read()
        self.assertIn(r'\begin{verbatim}', content)
        self.assertIn(r'\end{verbatim}', content)

    def test_file_contains_parameter_value(self):
        from anuga.validation_utilities.save_parameters_tex import save_parameters_tex
        domain = self._mock_domain({'flow_alg': 'DE1', 'cfl': 0.9})
        save_parameters_tex(domain)
        content = open('parameters.tex').read()
        self.assertIn('DE1', content)

    def test_calls_get_algorithm_parameters(self):
        from anuga.validation_utilities.save_parameters_tex import save_parameters_tex
        domain = self._mock_domain()
        save_parameters_tex(domain)
        domain.get_algorithm_parameters.assert_called_once()

    def test_accessible_from_package(self):
        from anuga.validation_utilities import save_parameters_tex
        self.assertTrue(callable(save_parameters_tex))


class TestProduceReport(unittest.TestCase):
    """produce_report orchestrates run_anuga_script × 2 + typeset_report."""

    def _make_args(self, verbose=False, np=1, no_run=False):
        args = MagicMock()
        args.verbose = verbose
        args.np = np
        # MagicMock hands back a truthy object for any attribute not set here,
        # so every flag produce_report reads has to be spelled out or the mock
        # silently takes the opposite branch.
        args.no_run = no_run
        return args

    def test_calls_run_script_then_plot_then_typeset(self):
        args = self._make_args()
        with patch.object(_pr_mod, 'run_anuga_script', return_value=0) as mock_run, \
             patch.object(_pr_mod, 'typeset_report') as mock_ts:
            _pr_mod.produce_report('myscript.py', args=args)

        self.assertEqual(mock_run.call_count, 2)
        self.assertEqual(mock_ts.call_count, 1)

    def test_first_run_uses_supplied_script(self):
        args = self._make_args()
        with patch.object(_pr_mod, 'run_anuga_script', return_value=0) as mock_run, \
             patch.object(_pr_mod, 'typeset_report'):
            _pr_mod.produce_report('sim.py', args=args)

        first_call_script = mock_run.call_args_list[0][0][0]
        self.assertEqual(first_call_script, 'sim.py')

    def test_second_run_uses_plot_results(self):
        args = self._make_args()
        with patch.object(_pr_mod, 'run_anuga_script', return_value=0) as mock_run, \
             patch.object(_pr_mod, 'typeset_report'):
            _pr_mod.produce_report('sim.py', args=args)

        second_call_script = mock_run.call_args_list[1][0][0]
        self.assertEqual(second_call_script, 'plot_results.py')

    def test_no_run_skips_the_simulation(self):
        """-nr reuses the existing sww: only plot_results.py runs."""
        args = self._make_args(no_run=True)
        with patch.object(_pr_mod, 'run_anuga_script', return_value=0) as mock_run, \
             patch.object(_pr_mod, 'typeset_report') as mock_ts:
            _pr_mod.produce_report('sim.py', args=args)

        self.assertEqual(mock_run.call_count, 1)
        self.assertEqual(mock_run.call_args_list[0][0][0], 'plot_results.py')
        self.assertEqual(mock_ts.call_count, 1)

    def test_np_set_to_1_for_plot_step(self):
        """produce_report forces np=1 for the plot_results.py step."""
        args = self._make_args(np=4)
        with patch.object(_pr_mod, 'run_anuga_script', return_value=0), \
             patch.object(_pr_mod, 'typeset_report'):
            _pr_mod.produce_report('sim.py', args=args)

        self.assertEqual(args.np, 1)

    def test_failed_simulation_exits_without_plotting(self):
        """A failing simulation is fatal, so the report driver sees it."""
        args = self._make_args()
        with patch.object(_pr_mod, 'run_anuga_script', return_value=1) as mock_run, \
             patch.object(_pr_mod, 'typeset_report') as mock_ts:
            with self.assertRaises(SystemExit) as cm:
                _pr_mod.produce_report('sim.py', args=args)

        self.assertIn('sim.py', str(cm.exception.code))
        self.assertEqual(mock_run.call_count, 1)
        self.assertEqual(mock_ts.call_count, 0)

    def test_failed_plot_exits_without_typesetting(self):
        args = self._make_args()
        with patch.object(_pr_mod, 'run_anuga_script', side_effect=[0, 2]), \
             patch.object(_pr_mod, 'typeset_report') as mock_ts:
            with self.assertRaises(SystemExit) as cm:
                _pr_mod.produce_report('sim.py', args=args)

        self.assertIn('plot_results.py', str(cm.exception.code))
        self.assertEqual(mock_ts.call_count, 0)

    def test_accessible_from_package(self):
        from anuga.validation_utilities import produce_report
        self.assertTrue(callable(produce_report))


if __name__ == '__main__':
    unittest.main()
