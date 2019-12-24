from docutils import nodes
from docutils.parsers.rst import Directive


class HelloWorld(Directive):
    has_content = True

    def run(self):
        paragraph_node = nodes.paragraph(text='Hello World!' + str(self.content))
        return [paragraph_node]


def setup(app):
    app.add_directive("helloworld", HelloWorld)

    return {
        'version': '0.1',
        'parallel_read_safe': True,
        'parallel_write_safe': True,
    }
