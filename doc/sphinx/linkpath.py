from docutils import nodes
from docutils.parsers.rst import Directive
import os


class LinkPath(Directive):
    has_content = True

    def run(self):
        def Try(abspath, *path):
            s = os.path.abspath(os.path.join(*path))
            if os.path.exists(s):
                assert not abspath, \
                    "Ambiguous path: target found in '{:}' and '{:}'".format(
                            abspath, s)
                return s
            return abspath

        assert len(self.content) == 1, \
            "Expected relative path, got '{:}'".format(self.content)
        r = self.content[0] # relative path

        repo = os.path.abspath(os.path.join(os.getcwd(), "../../"))
        a = "" # absolute path
        a = Try(a, repo, r)
        a = Try(a, repo, "deploy", r)
        a = Try(a, repo, "deploy/lib", r)
        a = Try(a, repo, "deploy/tool", r)
        s = Try(a, repo, "deploy/wrap", r)
        a = Try(a, repo, "src", r)

        assert a, "Target not found '{:}'".format(r)
        paragraph_node = nodes.paragraph(text=a)
        return [paragraph_node]


def setup(app):
    app.add_directive("linkpath", LinkPath)

    return {
        'version': '0.1',
        'parallel_read_safe': True,
        'parallel_write_safe': True,
    }
