"""A ``:spec:`` role linking the bracketed physics labels to their definitions.

Labels like ``[G-4]`` name a term in the sediment physics rather than a
reference, and they appear across the task page, the appendix and the source.
Written as plain literals they gave the reader a name without saying where to
find it, so ``:spec:`G-4``` renders ``[G-4]`` in the same monospace as before
and links it to the equation that defines it in the physics specification
(``appendices/physics_spec.md``), where every label has a ``(spec-x-n)=``
target on the block that defines it.

The target for label ``X-N`` is the RST label ``spec-x-n``. Resolution is done
here rather than through ``:ref:`` because the std domain rebuilds the link
text as a plain inline, which would drop the monospace and leave linked and
unlinked labels looking like different things. A label with no anchor is a
build warning, not silently dead text.
"""

from docutils import nodes
from sphinx.transforms.post_transforms import SphinxPostTransform
from sphinx.util import logging
from sphinx.util.docutils import SphinxRole

logger = logging.getLogger(__name__)


class spec_ref(nodes.Inline, nodes.Element):
    """Placeholder resolved once every document's labels are known."""


class SpecLitRole(SphinxRole):
    """``:speclit:`T-6``` -- the label, styled, but NOT a link.

    For the cell that DEFINES a label. Linking it would point at itself, but
    it still needs the ``spec-label`` class: without it the theme's
    ``code.literal { white-space: normal }`` applies and the browser breaks
    the label at its hyphen, rendering "[D-" on one line and "2]" on the next
    in a narrow column.
    """

    def run(self):
        label = self.text.strip()
        return [nodes.literal('', '[%s]' % label, classes=['spec-label'])], []


class SpecRole(SphinxRole):
    def run(self):
        node = spec_ref()
        node['spec'] = self.text.strip()
        node['docname'] = self.env.docname
        node.line = self.lineno
        return [node], []


class ResolveSpecRefs(SphinxPostTransform):
    default_priority = 5

    def run(self, **kwargs):
        std = self.env.get_domain('std')
        for node in list(self.document.findall(spec_ref)):
            label = node['spec']
            literal = nodes.literal('', '[%s]' % label, classes=['spec-label'])
            entry = std.anonlabels.get('spec-%s' % label.lower())
            if entry is None:
                logger.warning('no anchor for spec label [%s]; expected a '
                               '`(spec-%s)=` target in physics_spec.md',
                               label, label.lower(),
                               location=(node['docname'], node.line))
                node.replace_self(literal)
                continue
            docname, labelid = entry
            uri = self.app.builder.get_relative_uri(self.env.docname, docname)
            if labelid:
                uri += '#' + labelid
            ref = nodes.reference('', '', internal=True, refuri=uri,
                                  classes=['spec-ref'])
            ref += literal
            node.replace_self(ref)


def setup(app):
    app.add_role('spec', SpecRole())
    app.add_role('speclit', SpecLitRole())
    app.add_post_transform(ResolveSpecRefs)
    return {'version': '1.0',
            'parallel_read_safe': True,
            'parallel_write_safe': True}
