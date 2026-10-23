"""The Domain reference page lists its methods by hand -- keep it honest.

``docs/source/reference/anuga.Domain.rst`` is a *committed* autosummary stub:
unlike the stubs under ``docs/source/reference/generated/`` (git-ignored and
rebuilt every time), its ``~Domain.<name>`` list only changes when somebody
edits it. Adding a public method to Domain therefore leaves the reference page
silently incomplete -- Sphinx emits no warning for a method that was never
listed, so a clean docs build proves nothing.

That is how the page came to advertise ``set_multiprocessor_mode`` while
omitting its replacement ``set_compute_mode``, along with the whole tracer and
sediment API (anuga-community/anuga_core#319).

The docs live outside the installed package, so this skips unless run from a
source checkout.
"""

import re
import unittest
from pathlib import Path

import anuga

REFERENCE_RST = (Path(anuga.__file__).resolve().parents[1]
                 / 'docs' / 'source' / 'reference' / 'anuga.Domain.rst')


def _public_methods():
    import inspect
    return {name for name, _ in inspect.getmembers(anuga.Domain, callable)
            if not name.startswith('_')}


def _listed_methods():
    """The names under the "Methods" rubric, ignoring the Attributes one."""
    text = REFERENCE_RST.read_text(encoding='utf-8-sig')
    body = text.split('.. rubric:: Methods', 1)[1].split('.. rubric::', 1)[0]
    return set(re.findall(r'~Domain\.(\w+)', body))


@unittest.skipUnless(REFERENCE_RST.is_file(),
                     'needs a source checkout: %s' % REFERENCE_RST)
class TestDomainReferenceDocs(unittest.TestCase):

    def test_every_public_method_is_listed(self):
        missing = sorted(_public_methods() - _listed_methods())
        self.assertEqual(
            missing, [],
            'Domain gained public methods that the reference page does not '
            'list, so they are absent from the published docs. Add them to '
            '%s under "rubric:: Methods":\n    %s'
            % (REFERENCE_RST.name, '\n    '.join(missing)))

    def test_no_listed_method_has_been_removed(self):
        # __init__ is listed deliberately; it is documented via automethod.
        stale = sorted(_listed_methods() - _public_methods() - {'__init__'})
        self.assertEqual(
            stale, [],
            'The reference page lists methods that no longer exist on Domain. '
            'Remove them from %s:\n    %s'
            % (REFERENCE_RST.name, '\n    '.join(stale)))

    def test_the_current_compute_mode_api_is_documented(self):
        # The regression that prompted the guard: the superseded call was
        # listed and its replacement was not.
        listed = _listed_methods()
        absent = [n for n in ('set_compute_mode', 'get_compute_mode')
                  if n not in listed]
        self.assertEqual(absent, [],
                         'the current compute-mode API is missing from %s '
                         'while the superseded set_multiprocessor_mode is '
                         'listed' % REFERENCE_RST.name)


if __name__ == '__main__':
    unittest.main()
