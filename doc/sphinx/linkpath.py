import os

from docutils import nodes
from docutils.parsers.rst import Directive
from docutils.parsers.rst import directives

import sphinx
from sphinx.util.docutils import SphinxDirective

from findpath import FindPath, Assert

class LinkPath(SphinxDirective):
    has_content = False
    required_arguments = 1

    def run(self):
        location = (self.env.docname, self.lineno)
        args = self.arguments
        relpath = args[0]
        abspath = FindPath(relpath, location)
        Assert(abspath, "Target not found '{:}'".format(relpath), location)
        text = abspath
        p = nodes.problematic(text=text)
        return [p]


def setup(app):
    app.add_directive("linkpath", LinkPath)

    return {
        'version': '0.1',
        'parallel_read_safe': True,
        'parallel_write_safe': True,
    }
