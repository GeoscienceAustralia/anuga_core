"""Tests for anuga.scenario.emit_script (the `anuga_toml_run --emit-script`
generator)."""

import ast
import inspect
import os
import pkgutil
import tempfile
import unittest
from importlib import import_module
from pathlib import Path

import anuga.scenario
from anuga.scenario.emit_script import build_run_script, write_run_script


def _domain_setup_helpers():
    """Names of the setup_<x>(domain, project) helpers in anuga.scenario."""
    names = []
    for mod in pkgutil.iter_modules(anuga.scenario.__path__):
        if not mod.name.startswith('setup_'):
            continue
        func = getattr(import_module('anuga.scenario.' + mod.name), mod.name, None)
        if func is None:
            continue
        if list(inspect.signature(func).parameters)[:2] == ['domain', 'project']:
            names.append(mod.name)
    return sorted(names)


def _runner_source():
    """Text of scripts/anuga_toml_run.py from a source checkout, or None."""
    for start in (Path.cwd().resolve(), Path(__file__).resolve()):
        for parent in start.parents:
            candidate = parent / 'scripts' / 'anuga_toml_run.py'
            if candidate.is_file():
                return candidate.read_text(encoding='utf-8')
    return None


class Test_emit_script(unittest.TestCase):

    def test_generated_script_is_valid_python(self):
        """The emitted script must at least parse as Python."""
        text = build_run_script('scenario.toml')
        ast.parse(text)   # raises SyntaxError on failure

    def test_config_name_is_substituted(self):
        text = build_run_script('my_flood.toml')
        # Passed to PrepareData and shown in the docstring/example.
        self.assertIn("PrepareData('my_flood.toml'", text)
        self.assertIn('anuga_toml_run my_flood.toml', text)

    def test_drives_the_standard_phases(self):
        """All the runner's phases appear in the generated script, in order."""
        text = build_run_script('c.toml')
        phases = [
            'setup_mesh.setup_mesh(project)',
            'setup_initial_conditions.setup_initial_conditions(domain, project)',
            'setup_riverwalls.setup_riverwalls(domain, project)',
            'setup_rainfall.setup_rainfall(domain, project)',
            'setup_inlets.setup_inlets(domain, project)',
            'setup_culverts.setup_culverts(domain, project)',
            'setup_weirs.setup_weirs(domain, project)',
            'setup_erosion.setup_erosion(domain, project)',
            'setup_sediment.setup_sediment(domain, project)',
            'setup_boundary_conditions.setup_boundary_conditions(domain, project)',
            'domain.evolve(',
        ]
        last = -1
        for p in phases:
            idx = text.find(p)
            self.assertNotEqual(idx, -1, 'missing phase: %s' % p)
            self.assertGreater(idx, last, 'phase out of order: %s' % p)
            last = idx

    def test_calls_every_setup_helper(self):
        """A TOML section is ignored unless its setup_* helper is called, as
        [[culverts]] and [[weirs]] were until 2026-10. Check the emitted
        script and the runner against the helpers the package defines."""
        helpers = _domain_setup_helpers()
        self.assertIn('setup_culverts', helpers)
        sources = {'emitted script': build_run_script('c.toml')}
        runner = _runner_source()
        if runner is not None:
            sources['scripts/anuga_toml_run.py'] = runner
        for where, text in sources.items():
            for name in helpers:
                call = '%s.%s(domain, project)' % (name, name)
                self.assertTrue(call in text,
                                '%s never calls %s' % (where, name))

    def test_evolve_uses_config_timestepping(self):
        text = build_run_script('c.toml')
        self.assertIn('yieldstep=project.yieldstep', text)
        self.assertIn('finaltime=project.finaltime', text)

    def test_parallel_aware(self):
        """MPI-safe: barrier + rank-guarded output + parallel sww merge."""
        text = build_run_script('c.toml')
        self.assertIn('finalize()', text)
        self.assertIn('if numprocs > 1:', text)
        self.assertIn('domain.sww_merge(delete_old=True)', text)

    def test_script_example_uses_given_path(self):
        text = build_run_script('c.toml', script_path='/some/dir/my_run.py')
        self.assertIn('python my_run.py', text)

    def test_write_run_script_writes_file(self):
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, 'run.py')
            returned = write_run_script('c.toml', path)
            self.assertEqual(returned, path)
            self.assertTrue(os.path.exists(path))
            with open(path, encoding='utf-8') as fh:
                ast.parse(fh.read())


if __name__ == '__main__':
    unittest.main()
