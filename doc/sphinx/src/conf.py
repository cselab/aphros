import glob
import re
import os

extensions = [
              'sphinx.ext.githubpages',
              'sphinx.ext.imgmath',
              'sphinx.ext.viewcode',
              'includecode',
              ]

templates_path = []

source_suffix = '.rst'

master_doc = 'index'


project = 'Aphros'
#copyright = '2019, ETH Zurich'
#author = 'ETH Zurich'
author = ''

# short X.Y version
#version = "0.1"
# full version
#release = "0.1"

language = None

pygments_style = 'sphinx'

numfig = True

# html
html_theme = 'classic'

html_static_path = ['_static']

html_context = {'css_files': ['_static/center.css']}

html_theme_options = {
      'nosidebar' : True,
      'body_max_width' : '800px',
      }

html_sidebars = {
   #'**' : ['localtoc.html']
}

htmlhelp_basename = 'hydrodoc'

html_math_renderer = 'imgmath'
#html_math_renderer = ''

# latex
latex_elements = {}
latex_documents = [
    (master_doc, 'hydro.tex', u'Aphros documentation', author, 'manual'),
]

# man
man_pages = [
    (master_doc, 'hydro', u'Aphros documentation', '', 1)
]
lib = "lib/"
for p in glob.glob(lib + "*.rst"):
    if p in [lib + "index.rst"]:
        continue
    path = os.path.splitext(p)[0]
    name = os.path.basename(path)
    desc = ''
    with open(p) as f:
        line = f.readline()
        if line.startswith('..'):
            desc = line[3:]
    man_pages.append((path, name, desc, '', 1))

# texinfo
texinfo_documents = [
    (master_doc, 'hydro', u'Aphros documentation',
     author, 'hydro', 'One line description of project.',
     'Miscellaneous'),
]
