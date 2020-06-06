import os

from docutils import nodes
from docutils.parsers.rst import Directive
from docutils.parsers.rst import directives

import sphinx
from sphinx.util.docutils import SphinxDirective

from findpath import FindPath, Assert


def PathToGithub(relpath, root):
    if not relpath:
        return relpath
    return os.path.join(root, relpath)


class LinkPath(SphinxDirective):
    has_content = False
    required_arguments = 1

    def run(self):
        location = (self.env.docname, self.lineno)
        args = self.arguments
        relpath = args[0]
        abspath, relpath = FindPath(relpath, location)
        Assert(abspath, "Target not found '{:}'".format(relpath), location)
        config = self.env.config
        if config.linkpath_link_github:
            uri = PathToGithub(relpath, config.linkpath_github_root)
        else:
            uri = abspath
        return [nodes.problematic(text=relpath, refuri=uri)]


def linkpath_role(name,
                  rawtext,
                  text,
                  lineno,
                  inliner,
                  options={},
                  content=[]):
    """
    name: ?
    rawtext: ?
    text: arguments for supplied to role, relative path to file
    inliner: document instance
    lineno: line number
    """
    config = inliner.document.settings.env.app.config
    docname = os.path.splitext(inliner.document.attributes["source"])[0]
    location = (docname, lineno)
    relpath = text
    abspath, relpath = FindPath(relpath, location)
    Assert(abspath, "Target not found '{:}'".format(relpath), location)
    if config.linkpath_link_github:
        uri = PathToGithub(relpath, config.linkpath_github_root)
    else:
        uri = abspath
    return [nodes.reference(text=relpath, refuri=uri)], []


def setup(app):
    global LINK_GITHUB, GITHUB_ROOT
    app.add_config_value('linkpath_link_github', False, 'html')
    app.add_config_value('linkpath_github_root', None, 'html')
    #LINK_GITHUB = app.config.linkpath_link_github
    #GITHUB_ROOT = app.config.linkpath_github_root
    app.add_directive("linkpath", LinkPath)
    app.add_role("linkpath", linkpath_role)

    return {
        'version': '0.1',
        'parallel_read_safe': True,
        'parallel_write_safe': True,
    }
