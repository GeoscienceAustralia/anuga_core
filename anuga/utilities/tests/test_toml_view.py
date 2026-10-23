"""Unit tests for anuga.utilities.toml_view (the anuga_toml_view backend)."""
import textwrap
import unittest

from anuga.utilities.toml_view import collapse_toml, highlight, render


def _make(n_friction):
    lines = ['[project]', 'scenario = "t"', '',
             '[mesh]', 'default_res = 1000.0', '',
             '[[mesh.boundary_tags]]', 'tag = "south"', 'edges = [0]',
             '[[mesh.boundary_tags]]', 'tag = "east"', 'edges = [1]', '',
             '[initial_conditions]',
             '# Friction section']
    for i in range(n_friction):
        lines += ['[[initial_conditions.friction]]',
                  f'polygon = "p{i}.csv"', 'value = 0.15']
    return '\n'.join(lines) + '\n'


class TestCollapse(unittest.TestCase):
    def test_collapses_long_run(self):
        col = collapse_toml(_make(20), threshold=6)
        # first friction block kept, rest collapsed to a marker
        self.assertEqual(col.count('[[initial_conditions.friction]]'), 2)  # 1 block + marker text
        self.assertIn('19 more [[initial_conditions.friction]] blocks collapsed', col)
        self.assertLess(col.count('\n'), _make(20).count('\n'))

    def test_short_run_not_collapsed(self):
        # boundary_tags has only 2 entries -> untouched at default threshold
        col = collapse_toml(_make(3), threshold=6)
        self.assertEqual(col.count('[[mesh.boundary_tags]]'), 2)  # both kept
        # 3 friction < threshold 6 -> not collapsed, no marker
        self.assertNotIn('collapsed', col)

    def test_threshold_respected(self):
        col = collapse_toml(_make(4), threshold=3)
        self.assertIn('collapsed', col)             # 4 > 3
        col2 = collapse_toml(_make(4), threshold=4)
        self.assertNotIn('collapsed', col2)         # 4 not > 4

    def test_section_comment_travels_with_block(self):
        col = collapse_toml(_make(20), threshold=6)
        # the '# Friction section' comment must still precede the kept block
        idx_c = col.index('# Friction section')
        idx_b = col.index('[[initial_conditions.friction]]')
        self.assertLess(idx_c, idx_b)

    def test_marker_is_a_toml_comment(self):
        # markers start with '#', so they stay valid-ish and highlight as comments
        for ln in collapse_toml(_make(20), threshold=6).splitlines():
            if 'collapsed' in ln:
                self.assertTrue(ln.lstrip().startswith('#'))

    def test_no_headers_passthrough(self):
        txt = 'a = 1\nb = 2\n'
        self.assertEqual(collapse_toml(txt).strip(), txt.strip())

    def test_singular_plural_marker(self):
        # exactly threshold+? -> a run of 8 collapses to "7 more ... blocks"
        col = collapse_toml(_make(8), threshold=6)
        self.assertIn('7 more', col)
        self.assertIn('blocks collapsed', col)


class TestHighlightRender(unittest.TestCase):
    def test_highlight_off_is_identity(self):
        txt = _make(3)
        self.assertEqual(highlight(txt, color=False), txt)

    def test_highlight_on_returns_text(self):
        # pygments may or may not be installed; either way we get a str back
        out = highlight('[project]\nx = 1\n', color=True)
        self.assertIsInstance(out, str)
        self.assertIn('project', out)

    def test_render_collapses_and_is_string(self):
        out = render(_make(20), collapse=True, color=False, threshold=6)
        self.assertIn('collapsed', out)
        out_full = render(_make(20), collapse=False, color=False)
        self.assertNotIn('collapsed', out_full)


if __name__ == '__main__':
    unittest.main()
