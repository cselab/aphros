"""
    sphinx.ext.includecode
    ~~~~~~~~~~~~~~~~~~~

    Provides the ``includecode`` directive that allows to write documentation
    that is included depending on configuration variables.

    Usage::

        .. includecode:: releaselevel in ('alpha', 'beta', 'rc')

           This stuff is only included in the built docs for unstable versions.

    The argument for ``includecode`` is a plain Python expression, evaluated in the
    namespace of the project configuration (that is, all variables from
    ``conf.py`` are available.)

    :copyright: Copyright 2007-2019 by the Sphinx team, see AUTHORS.
    :license: BSD, see LICENSE for details.
"""
import sys
import warnings

from docutils import nodes
from docutils.parsers.rst import directives
from docutils.statemachine import StringList

import sphinx
from sphinx.util.docutils import SphinxDirective
from sphinx.util.nodes import nested_parse_with_titles
from sphinx.directives.code import container_wrapper

from sphinx import addnodes
from sphinx.deprecation import RemovedInSphinx40Warning
from sphinx.locale import __
from sphinx.util import logging
from sphinx.util import parselinenos
from sphinx.util.docutils import SphinxDirective

if False:
    # For type annotation
    from typing import Any, Dict, List  # NOQA
    from sphinx.application import Sphinx  # NOQA


class includecode(nodes.Element):
    pass

class LiteralIncludeReader:
    INVALID_OPTIONS_PAIR = []

    def __init__(self, filename, options, config):
        # type: (str, Dict, Config) -> None
        self.filename = filename
        self.options = options
        self.encoding = options.get('encoding', config.source_encoding)

        self.parse_options()

    def parse_options(self):
        # type: () -> None
        for option1, option2 in self.INVALID_OPTIONS_PAIR:
            if option1 in self.options and option2 in self.options:
                raise ValueError(__('Cannot use both "%s" and "%s" options') %
                                 (option1, option2))

    def read_file(self, filename, location=None):
        # type: (str, Tuple[str, int]) -> List[str]
        try:
            with open(filename, encoding=self.encoding, errors='strict') as f:
                text = f.read()
                if 'tab-width' in self.options:
                    text = text.expandtabs(self.options['tab-width'])

                return text.splitlines(True)
        except OSError:
            raise OSError(__('Include file %r not found or reading it failed') % filename)
        except UnicodeError:
            raise UnicodeError(__('Encoding %r used for reading included file %r seems to '
                                  'be wrong, try giving an :encoding: option') %
                               (self.encoding, filename))

    def read(self, location=None):
        # type: (Tuple[str, int]) -> Tuple[str, int]
        filters = [self.pyobject_filter,
                   self.lines_filter,
                   self.prepend_filter,
                   self.append_filter,
                   self.dedent_filter]
        lines = self.read_file(self.filename, location=location)
        for func in filters:
            lines = func(lines, location=location)

        return ''.join(lines), len(lines)


    def pyobject_filter(self, lines, location=None):
        # type: (List[str], Tuple[str, int]) -> List[str]
        pyobject = self.options.get('pyobject')
        if pyobject:
            from sphinx.pycode import ModuleAnalyzer
            analyzer = ModuleAnalyzer.for_file(self.filename, '')
            tags = analyzer.find_tags()
            if pyobject not in tags:
                raise ValueError(__('Object named %r not found in include file %r') %
                                 (pyobject, self.filename))
            else:
                start = tags[pyobject][1]
                end = tags[pyobject][2]
                lines = lines[start - 1:end]

        return lines

    def lines_filter(self, lines, location=None):
        # type: (List[str], Tuple[str, int]) -> List[str]
        linespec = self.options.get('lines')
        if linespec:
            linelist = parselinenos(linespec, len(lines))
            if any(i >= len(lines) for i in linelist):
                logger.warning(__('line number spec is out of range(1-%d): %r') %
                               (len(lines), linespec), location=location)

            lines = [lines[n] for n in linelist if n < len(lines)]
            if lines == []:
                raise ValueError(__('Line spec %r: no lines pulled from include file %r') %
                                 (linespec, self.filename))

        return lines

    def prepend_filter(self, lines, location=None):
        # type: (List[str], Tuple[str, int]) -> List[str]
        prepend = self.options.get('prepend')
        if prepend:
            lines.insert(0, prepend + '\n')

        return lines

    def append_filter(self, lines, location=None):
        # type: (List[str], Tuple[str, int]) -> List[str]
        append = self.options.get('append')
        if append:
            lines.append(append + '\n')

        return lines

    def dedent_filter(self, lines, location=None):
        # type: (List[str], Tuple[str, int]) -> List[str]
        if 'dedent' in self.options:
            return dedent_lines(lines, self.options.get('dedent'), location=location)
        else:
            return lines



class IncludeCode(SphinxDirective):
    has_content = True
    required_arguments = 1
    optional_arguments = 0
    final_argument_whitespace = True
    option_spec = {
        'dedent': int,
        'linenos': directives.flag,
        'tab-width': int,
        'language': directives.unchanged_required,
        'encoding': directives.encoding,
        'pyobject': directives.unchanged_required,
        'lines': directives.unchanged_required,
        'start-after': directives.unchanged_required,
        'end-before': directives.unchanged_required,
        'start-at': directives.unchanged_required,
        'end-at': directives.unchanged_required,
        'prepend': directives.unchanged_required,
        'append': directives.unchanged_required,
        'emphasize-lines': directives.unchanged_required,
        'caption': directives.unchanged,
        'class': directives.class_option,
        'name': directives.unchanged,
    }


    def run(self):
        # type: () -> List[nodes.Node]
        document = self.state.document
        if not document.settings.file_insertion_enabled:
            return [document.reporter.warning('File insertion disabled',
                                              line=self.lineno)]

        try:
            location = self.state_machine.get_source_and_line(self.lineno)
            rel_filename, filename = self.env.relfn2path(self.arguments[0])
            self.env.note_dependency(rel_filename)

            reader = LiteralIncludeReader(filename, self.options, self.config)
            text, lines = reader.read(location=location)

            retnode = nodes.literal_block(text, text, source=filename)  # type: nodes.Element
            #self.set_source_info(retnode)
            if 'language' in self.options:
                retnode['language'] = self.options['language']
            if ('linenos' in self.options or 'lineno-start' in self.options or
                    'lineno-match' in self.options):
                retnode['linenos'] = True
            retnode['classes'] += self.options.get('class', [])
            extra_args = retnode['highlight_args'] = {}
            if 'emphasize-lines' in self.options:
                hl_lines = parselinenos(self.options['emphasize-lines'], lines)
                if any(i >= lines for i in hl_lines):
                    logger.warning(__('line number spec is out of range(1-%d): %r') %
                                   (lines, self.options['emphasize-lines']),
                                   location=location)
                extra_args['hl_lines'] = [x + 1 for x in hl_lines if x < lines]

            if 'caption' in self.options:
                caption = self.options['caption'] or self.arguments[0]
                retnode = container_wrapper(self, retnode, caption)

            # retnode will be note_implicit_target that is linked from caption and numref.
            # when options['name'] is provided, it should be primary ID.
            self.add_name(retnode)

            return [retnode]
        except Exception as exc:
            return [document.reporter.warning(exc, line=self.lineno)]


def process_includecode_nodes(app, doctree, docname):
    # type: (Sphinx, nodes.document, str) -> None
    ns = {confval.name: confval.value for confval in app.config}
    ns.update(app.config.__dict__.copy())
    ns['builder'] = app.builder.name
    for node in doctree.traverse(includecode):
        try:
            res = eval(node['expr'], ns)
        except Exception as err:
            # handle exceptions in a clean fashion
            from traceback import format_exception_only
            msg = ''.join(format_exception_only(err.__class__, err))
            newnode = doctree.reporter.error('Exception occured in '
                                             'includecode expression: \n%s' %
                                             msg, base_node=node)
            node.replace_self(newnode)
        else:
            if not res:
                node.replace_self([])
            else:
                node.replace_self(node.children)


def setup(app):
    # type: (Sphinx) -> Dict[str, Any]
    app.add_node(includecode)
    app.add_directive('includecode', IncludeCode)
    app.connect('doctree-resolved', process_includecode_nodes)
    return {'version': sphinx.__display_version__, 'parallel_read_safe': True}
