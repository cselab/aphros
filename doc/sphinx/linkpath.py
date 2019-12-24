import os

from docutils import nodes
from docutils.parsers.rst import Directive

from findpath import FindPath


class LinkPath(Directive):
    has_content = True

    def run(self):
        assert len(self.content) == 1, \
            "Expected relative path, got '{:}'".format(self.content)
        relpath = self.content[0]
        abspath = FindPath(relpath)
        assert abspath, "Target not found '{:}'".format(relpath)
        return [nodes.line(text=abspath)]


def setup(app):
    app.add_directive("linkpath", LinkPath)

    return {
        'version': '0.1',
        'parallel_read_safe': True,
        'parallel_write_safe': True,
    }
