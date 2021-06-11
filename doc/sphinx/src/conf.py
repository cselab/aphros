import glob
import re
import os
import sys

sys.path.append(os.path.abspath("./"))

extensions = [
              'sphinx.ext.githubpages',
              'sphinx.ext.imgmath',
              'sphinx.ext.viewcode',
              'includecode',
              'linkpath',
              ]

templates_path = []
today_fmt = ' '

source_suffix = '.rst'

master_doc = 'index'


project = 'Aphros'
copyright = '2020 ETH Zurich'
#author = 'ETH Zurich'
author = ''

linkpath_github_root = "https://github.com/cselab/aphros/blob/master/"
linkpath_link_github = True

# short X.Y version
#version = "0.1"
# full version
#release = "0.1"

language = None

pygments_style = 'sphinx'

numfig = True

# html
html_theme = 'classic'

html_theme_options = {
      'body_max_width' : '800px',
      }

html_sidebars = {
   #'**' : ['localtoc.html']
}

htmlhelp_basename = 'hydrodoc'

html_math_renderer = 'imgmath'
imgmath_image_format = 'svg'
imgmath_font_size = 15
imgmath_use_preview = True

# latex
latex_elements = {}
latex_documents = [
    (master_doc, 'hydro.tex', u'Aphros documentation', author, 'manual'),
]

# man
man_pages = [
    (master_doc, 'hydro', u'Aphros documentation', '', 1),
    ('log/index', 'ap.log', u'development log', '', 7),
    ('styleguide', 'ap.styleguide', u'coding style guide', '', 7),
]


for lib in glob.glob("lib?"):
    sec = int(lib[3])
    for p in glob.glob(lib + "/*.rst"):
        if p in [lib + "/index.rst"]:
            continue
        path = os.path.splitext(p)[0]
        name = os.path.basename(path)
        desc = ''
        with open(p) as f:
            line = f.readline()
            if line.startswith('..'):
                desc = line[3:]
        man_pages.append((path, name, desc, '', sec))

# texinfo
texinfo_documents = [
    (master_doc, 'hydro', u'Aphros documentation',
     author, 'hydro', 'One line description of project.',
     'Miscellaneous'),
]


todo_include_todos = True
