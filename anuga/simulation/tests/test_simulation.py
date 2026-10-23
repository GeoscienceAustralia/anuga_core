"""
Tests for anuga.simulation.simulation

The Simulation class itself is tightly coupled to the MPI/parallel
infrastructure (sequential_distribute_dump/load, barrier, finalize) and is
integration-test territory.  These unit tests cover the two standalone
argument-parsing helpers:

  parse_args                  -- requires a callable argument_adder
  parse_args_and_parameters   -- argument_adder is optional; silently skips a
                                 missing project.py
"""

import sys
import unittest

from anuga.simulation.simulation import parse_args, parse_args_and_parameters


class TestParseArgs(unittest.TestCase):
    """Tests for parse_args (requires a callable argument_adder)."""

    @staticmethod
    def _noop_adder(parser):
        """No-op argument_adder — adds nothing extra."""

    def test_default_alg(self):
        args = parse_args(self._noop_adder, from_commandline=False)
        self.assertEqual(args.alg, 'DE1')

    def test_default_np(self):
        args = parse_args(self._noop_adder, from_commandline=False)
        self.assertEqual(args.np, 1)

    def test_default_verbose_false(self):
        args = parse_args(self._noop_adder, from_commandline=False)
        self.assertFalse(args.verbose)

    def test_default_checkpointing_false(self):
        args = parse_args(self._noop_adder, from_commandline=False)
        self.assertFalse(args.checkpointing)

    def test_default_outname(self):
        args = parse_args(self._noop_adder, from_commandline=False)
        self.assertEqual(args.outname, 'domain')

    def test_default_partition_dir(self):
        args = parse_args(self._noop_adder, from_commandline=False)
        self.assertEqual(args.partition_dir, 'PARTITIONS')

    def test_default_checkpoint_dir(self):
        args = parse_args(self._noop_adder, from_commandline=False)
        self.assertEqual(args.checkpoint_dir, 'CHECKPOINTS')

    def test_default_checkpoint_time(self):
        args = parse_args(self._noop_adder, from_commandline=False)
        self.assertAlmostEqual(args.checkpoint_time, -1.0)

    def test_kwarg_overrides_alg(self):
        args = parse_args(self._noop_adder, from_commandline=False, alg='DE0')
        self.assertEqual(args.alg, 'DE0')

    def test_kwarg_overrides_outname(self):
        args = parse_args(self._noop_adder, from_commandline=False,
                          outname='my_run')
        self.assertEqual(args.outname, 'my_run')

    def test_kwarg_overrides_np(self):
        args = parse_args(self._noop_adder, from_commandline=False, np=4)
        self.assertEqual(args.np, 4)

    def test_extra_argument_from_adder(self):
        def adder(parser):
            parser.add_argument('--extra', type=int, default=42)

        args = parse_args(adder, from_commandline=False)
        self.assertEqual(args.extra, 42)

    def test_extra_argument_kwarg_override(self):
        def adder(parser):
            parser.add_argument('--extra', type=int, default=42)

        args = parse_args(adder, from_commandline=False, extra=99)
        self.assertEqual(args.extra, 99)

    def test_returns_namespace(self):
        import argparse
        args = parse_args(self._noop_adder, from_commandline=False)
        self.assertIsInstance(args, argparse.Namespace)


class TestParseArgsAndParameters(unittest.TestCase):
    """Tests for parse_args_and_parameters (argument_adder is optional)."""

    def test_no_argument_adder_returns_namespace(self):
        import argparse
        args = parse_args_and_parameters(argument_adder=None,
                                         from_commandline=False)
        self.assertIsInstance(args, argparse.Namespace)

    def test_default_alg(self):
        args = parse_args_and_parameters(from_commandline=False)
        self.assertEqual(args.alg, 'DE1')

    def test_default_np(self):
        args = parse_args_and_parameters(from_commandline=False)
        self.assertEqual(args.np, 1)

    def test_default_verbose_false(self):
        args = parse_args_and_parameters(from_commandline=False)
        self.assertFalse(args.verbose)

    def test_default_outname(self):
        args = parse_args_and_parameters(from_commandline=False)
        self.assertEqual(args.outname, 'domain')

    def test_kwarg_overrides_alg(self):
        args = parse_args_and_parameters(from_commandline=False, alg='DE2')
        self.assertEqual(args.alg, 'DE2')

    def test_kwarg_overrides_outname(self):
        args = parse_args_and_parameters(from_commandline=False,
                                         outname='tsunami_run')
        self.assertEqual(args.outname, 'tsunami_run')

    def test_kwarg_overrides_np(self):
        args = parse_args_and_parameters(from_commandline=False, np=8)
        self.assertEqual(args.np, 8)

    def test_project_py_absent_does_not_raise(self):
        """If no project.py is importable the function must still return defaults."""
        # Ensure project is not importable by temporarily hiding it
        had_project = 'project' in sys.modules
        sys.modules.pop('project', None)

        try:
            args = parse_args_and_parameters(from_commandline=False)
            self.assertEqual(args.outname, 'domain')
        finally:
            if had_project:
                pass  # don't restore a potentially broken cached module

    def test_with_argument_adder(self):
        def adder(parser):
            parser.add_argument('--mode', type=str, default='serial')

        args = parse_args_and_parameters(argument_adder=adder,
                                         from_commandline=False)
        self.assertEqual(args.mode, 'serial')

    def test_kwargs_override_adder_default(self):
        def adder(parser):
            parser.add_argument('--mode', type=str, default='serial')

        args = parse_args_and_parameters(argument_adder=adder,
                                         from_commandline=False,
                                         mode='parallel')
        self.assertEqual(args.mode, 'parallel')


if __name__ == '__main__':
    unittest.main()
