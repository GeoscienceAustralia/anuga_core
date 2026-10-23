"""
Unit tests for anuga.scenario.setup_boundary_conditions

Tests cover:
  - Reflective boundary applied to all tags
  - Stage boundary with a CSV timeseries file
  - Flather_Stage boundary with a CSV timeseries file
  - Unknown boundary type raises NotImplementedError
  - Unset boundary tag raises Exception
"""
import os
import tempfile
import shutil
import unittest
from unittest.mock import MagicMock

try:
    import anuga
    from anuga import rectangular_cross_domain
    from anuga.scenario.setup_boundary_conditions import setup_boundary_conditions
    HAS_MODULE = True
    SKIP_REASON = ''
except ImportError as _e:
    HAS_MODULE = False
    SKIP_REASON = str(_e)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

ALL_TAGS = ['left', 'right', 'bottom', 'top']


def _make_domain():
    """Return a tiny 3×3 rectangular cross domain with named boundary tags."""
    return rectangular_cross_domain(3, 3)


def _make_project(boundary_data, boundary_tags=None):
    """Return a simple mock for the project object.

    If *boundary_tags* is omitted it defaults to the tags present in
    *boundary_data* (first element of each entry).
    """
    if boundary_tags is None:
        boundary_tags = [bd[0] for bd in boundary_data]
    p = MagicMock()
    p.boundary_data = boundary_data
    p.boundary_tags = boundary_tags
    return p


def _all_reflective(domain):
    """Return boundary_data that maps every domain tag to Reflective."""
    return [[t, 'Reflective'] for t in ALL_TAGS]


def _write_stage_csv(path, times, stages):
    """Write a two-column CSV (time, stage) with a header row."""
    with open(path, 'w') as f:
        f.write('time,stage\n')
        for t, s in zip(times, stages):
            f.write('%s,%s\n' % (t, s))


def _boundary_types(domain):
    """Return a set of boundary object types attached to the domain."""
    return {type(B) for _, B in domain.boundary_objects}


def _boundary_for_tag(domain, tag):
    """Return the first boundary object whose edge has the given tag."""
    for (vol_id, edge_id), B in domain.boundary_objects:
        if domain.boundary[(vol_id, edge_id)] == tag:
            return B
    return None


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

@unittest.skipUnless(HAS_MODULE, SKIP_REASON)
class TestSetupBoundaryConditions(unittest.TestCase):

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        self.domain = _make_domain()

    def tearDown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def _stage_csv(self, name='stage.csv'):
        return os.path.join(self.tmpdir, name)

    # --- Reflective ---

    def test_reflective_all_tags_applied(self):
        project = _make_project(_all_reflective(self.domain), ALL_TAGS)
        setup_boundary_conditions(self.domain, project)
        types = _boundary_types(self.domain)
        self.assertIn(anuga.Reflective_boundary, types)

    def test_reflective_every_edge_has_correct_type(self):
        project = _make_project(_all_reflective(self.domain), ALL_TAGS)
        setup_boundary_conditions(self.domain, project)
        for tag in ALL_TAGS:
            B = _boundary_for_tag(self.domain, tag)
            self.assertIsNotNone(B, msg='No boundary for tag ' + tag)
            self.assertIsInstance(B, anuga.Reflective_boundary)

    def test_reflective_sets_boundary_objects(self):
        project = _make_project(_all_reflective(self.domain), ALL_TAGS)
        setup_boundary_conditions(self.domain, project)
        self.assertGreater(len(self.domain.boundary_objects), 0)

    # --- Stage ---

    def test_stage_boundary_left_is_not_reflective(self):
        csv = self._stage_csv()
        _write_stage_csv(csv, times=[0, 3600, 7200], stages=[0.0, 1.0, 0.5])
        bd = [
            ['left', 'Stage', csv, 0.0],
            ['right', 'Reflective'],
            ['bottom', 'Reflective'],
            ['top', 'Reflective'],
        ]
        project = _make_project(bd, ALL_TAGS)
        setup_boundary_conditions(self.domain, project)
        B_left = _boundary_for_tag(self.domain, 'left')
        self.assertIsNotNone(B_left)
        self.assertNotIsInstance(B_left, anuga.Reflective_boundary)

    def test_stage_boundary_other_tags_are_reflective(self):
        csv = self._stage_csv()
        _write_stage_csv(csv, times=[0, 3600], stages=[0.0, 1.0])
        bd = [
            ['left', 'Stage', csv, 0.0],
            ['right', 'Reflective'],
            ['bottom', 'Reflective'],
            ['top', 'Reflective'],
        ]
        project = _make_project(bd, ALL_TAGS)
        setup_boundary_conditions(self.domain, project)
        for tag in ['right', 'bottom', 'top']:
            self.assertIsInstance(_boundary_for_tag(self.domain, tag),
                                  anuga.Reflective_boundary)

    def test_stage_boundary_start_time_offset(self):
        # Timeseries starts at t=1000; domain starts at t=0
        csv = self._stage_csv('offset.csv')
        _write_stage_csv(csv, times=[1000, 2000, 3000], stages=[0.5, 1.0, 1.5])
        bd = [['left', 'Stage', csv, 1000.0],
              ['right', 'Reflective'],
              ['bottom', 'Reflective'],
              ['top', 'Reflective']]
        project = _make_project(bd, ALL_TAGS)
        # Should not raise (time offset converts 1000→0, 2000→1000, 3000→2000)
        setup_boundary_conditions(self.domain, project)

    # --- Flather_Stage ---

    def test_flather_stage_boundary_applied(self):
        csv = self._stage_csv('flather.csv')
        _write_stage_csv(csv, times=[0, 3600], stages=[0.0, 0.5])
        bd = [
            ['left', 'Flather_Stage', csv, 0.0],
            ['right', 'Reflective'],
            ['bottom', 'Reflective'],
            ['top', 'Reflective'],
        ]
        project = _make_project(bd, ALL_TAGS)
        setup_boundary_conditions(self.domain, project)
        B_left = _boundary_for_tag(self.domain, 'left')
        self.assertIsNotNone(B_left)
        self.assertNotIsInstance(B_left, anuga.Reflective_boundary)

    # --- Error cases ---

    def test_unknown_boundary_type_raises(self):
        # Can only provide 'left' boundary since set_boundary needs all 4 —
        # unknown type raises before set_boundary is called
        project = _make_project(
            [['left', 'FancyNewBoundary'],
             ['right', 'Reflective'],
             ['bottom', 'Reflective'],
             ['top', 'Reflective']],
            ALL_TAGS,
        )
        with self.assertRaises(NotImplementedError):
            setup_boundary_conditions(self.domain, project)

    def test_unset_boundary_tag_raises(self):
        # boundary_data only covers 3 of 4 tags; set_boundary should raise
        project = _make_project(
            [['left', 'Reflective'],
             ['right', 'Reflective'],
             ['bottom', 'Reflective']],
            boundary_tags=['left', 'right', 'bottom', 'top'],  # 'top' missing
        )
        with self.assertRaises(Exception):
            setup_boundary_conditions(self.domain, project)

    def test_mixed_reflective_and_stage(self):
        csv = self._stage_csv('mix.csv')
        _write_stage_csv(csv, times=[0, 1000], stages=[0.1, 0.2])
        bd = [
            ['left',   'Stage',      csv, 0.0],
            ['right',  'Reflective'],
            ['bottom', 'Reflective'],
            ['top',    'Reflective'],
        ]
        project = _make_project(bd, ALL_TAGS)
        setup_boundary_conditions(self.domain, project)
        self.assertNotIsInstance(_boundary_for_tag(self.domain, 'left'),
                                 anuga.Reflective_boundary)
        for tag in ['right', 'bottom', 'top']:
            self.assertIsInstance(_boundary_for_tag(self.domain, tag),
                                  anuga.Reflective_boundary)


if __name__ == '__main__':
    unittest.main()
