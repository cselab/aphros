import sys
import warnings
import re

from docutils import nodes
from docutils.parsers.rst import directives
from docutils.statemachine import StringList

import sphinx
from sphinx.util.docutils import SphinxDirective
from sphinx.util.nodes import nested_parse_with_titles
from sphinx.directives.code import container_wrapper, dedent_lines

from sphinx import addnodes
from sphinx.locale import __
from sphinx.util import logging
from sphinx.util import parselinenos
from sphinx.util.docutils import SphinxDirective

from findpath import FindPath

if False:
    # For type annotation
    from typing import Any, Dict, List  # NOQA
    from sphinx.application import Sphinx  # NOQA


class includecode(nodes.Element):
    pass


def ExtractBraces(text):
    pos = 0
    s = 0
    r = []

    def peek():
        return text[pos] if pos < len(text) else ''

    while pos < len(text):
        c = text[pos]
        pos += 1
        if s == 0:  # inside outer brace
            if c == '{':
                cnt = 1
                s = 1
        elif s == 1:  # inside inner brace
            if c == '{':
                cnt += 1
            if c == '}':
                cnt -= 1
                if cnt == 0:
                    break
            if c == '"':
                s = 2
            if c == '/' and peek() == '/':
                s = 3
            if c == '/' and peek() == '*':
                s = 4
            if c == "'":
                pos += 1
        elif s == 2:  # inside string
            if c == '\\':
                pos += 1
            if c == '"':
                s = 1
        elif s == 3:  # inside // comment
            if c == '\\':
                pos += 1
            if c == '\n':
                s = 1
        elif s == 4:  # inside /* comment
            if c == '*' and peek() == '/':
                s = 1
    if pos + 1 < len(text) and text[pos] == ';':
        pos += 1
    return text[:pos]


# Extracts function prototype from text
# text: text
# func: function name
# comment: if True, include comments above the prototype
# impl: if True, adds function implementation
def GetFunc(text, func, filename, comment=True, impl=False):
    p = "(^.*?\s+" + func + "\([^\)]*?\)[^:\{;=\$]*)([:\{;=\$]*[^\}]+\})?"
    m = re.search(p, text, re.MULTILINE)
    if not func or not m:
        res = "      // Error: function '{:}' not found in '{:}'".format(
            func, filename)
    else:
        if comment:
            s = m.start(0)
            ll = text[:s].splitlines(True)
            i = 0
            while i > -len(ll) and re.search(" *//", ll[i - 1]):
                i -= 1
            res = "".join(ll[i:]) if i != 0 else ""
        else:
            res = ""
        res += m.group(1)
        if impl:
            if m.lastindex > 1:
                res += ExtractBraces(text[m.end(1):])
            else:
                res += "\n      // Error: impl=True but implementation not found'".\
                        format(func, filename)
    return res.rstrip()


# Extracts struct definition.
# text: text
# f: struct name
def GetStruct(text, f):
    p = "(^.*struct " + f + " {[\w\W]*?};.*?$)"
    m = re.search(p, text, re.MULTILINE)
    if not f or not m:
        text = "// Error: struct '{:}' not found in '{:}'".format(f, filename)
    else:
        text = ExtractBraces(text[m.start(1):])
    return text


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
                raise ValueError(
                    __('Cannot use both "%s" and "%s" options') %
                    (option1, option2))

    def read_file(self, filename, location=None):
        # type: (str, Tuple[str, int]) -> List[str]
        try:
            with open(filename, encoding=self.encoding, errors='strict') as f:
                text = f.read()
                if 'func' in self.options:
                    func = self.options['func'].strip()
                    comment = 'comment' in self.options
                    impl = 'impl' in self.options
                    text = GetFunc(text, func, filename, comment, impl)
                if 'struct' in self.options:
                    f = self.options['struct'].strip()
                    text = GetStruct(text, f)
                if 'tab-width' in self.options:
                    text = text.expandtabs(self.options['tab-width'])

                return text.splitlines(True)
        except OSError:
            raise OSError(
                __('Include file %r not found or reading it failed') %
                filename)
        except UnicodeError:
            raise UnicodeError(
                __('Encoding %r used for reading included file %r seems to '
                   'be wrong, try giving an :encoding: option') %
                (self.encoding, filename))

    def read(self, location=None):
        # type: (Tuple[str, int]) -> Tuple[str, int]
        filters = [self.prepend_filter, self.append_filter, self.dedent_filter]
        lines = self.read_file(self.filename, location=location)
        for func in filters:
            lines = func(lines, location=location)

        return ''.join(lines), len(lines)

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
            return dedent_lines(lines,
                                self.options.get('dedent'),
                                location=location)
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
        'comment': directives.flag,
        'impl': directives.flag,
        'tab-width': int,
        'language': directives.unchanged_required,
        'encoding': directives.encoding,
        'start-after': directives.unchanged_required,
        'end-before': directives.unchanged_required,
        'start-at': directives.unchanged_required,
        'end-at': directives.unchanged_required,
        'prepend': directives.unchanged_required,
        'append': directives.unchanged_required,
        'caption': directives.unchanged,
        'class': directives.class_option,
        'name': directives.unchanged,
        'func': directives.unchanged,
        'struct': directives.unchanged,
    }

    def run(self):
        # type: () -> List[nodes.Node]
        document = self.state.document
        if not document.settings.file_insertion_enabled:
            return [
                document.reporter.warning('File insertion disabled',
                                          line=self.lineno)
            ]

        try:
            location = self.state_machine.get_source_and_line(self.lineno)
            relpath = self.arguments[0]
            abspath, relpath = FindPath(relpath,
                                        (self.env.docname, self.lineno))
            self.env.note_dependency(abspath)

            reader = LiteralIncludeReader(abspath, self.options, self.config)
            text, lines = reader.read(location=location)

            retnode = nodes.literal_block(
                text, text, source=abspath)  # type: nodes.Element
            #self.set_source_info(retnode)
            if 'language' not in self.options:
                self.options['language'] = 'cpp'
            if 'language' in self.options:
                retnode['language'] = self.options['language']
            if ('linenos' in self.options):
                retnode['linenos'] = True
            retnode['classes'] += self.options.get('class', [])
            extra_args = retnode['highlight_args'] = {}

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
                                             msg,
                                             base_node=node)
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
