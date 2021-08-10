import glob
import re
import os
import sys

sys.path.append(os.path.abspath('..'))
extensions = [
              'sphinx.ext.githubpages',
              'sphinx.ext.imgmath',
              'sphinx.ext.viewcode',
              'includecode',
              'linkpath',
              'sphinxcontrib.bibtex',
              ]

templates_path = []
today_fmt = ' '

source_suffix = '.rst'

master_doc = 'index'

project = 'Aphros'
copyright = '2021, ETH Zurich'
author = 'Petr Karnakov, Sergey Litvinov'

linkpath_github_root = 'https://github.com/cselab/aphros/blob/master/'
linkpath_link_github = True

bibtex_bibfiles = ['main.bib']
# short X.Y version
version = '0.1.2'
# full version
release = '0.1.2'

language = 'en'

pygments_style = 'sphinx'

numfig = True

# html
html_theme = 'theme'
html_theme_path = ['.']

html_theme_options = {
      'body_max_width' : '800px',
      }

html_sidebars = {
   #'**' : ['localtoc.html']
}

htmlhelp_basename = 'aphrosdoc'

html_math_renderer = 'imgmath'
imgmath_image_format = 'svg'
imgmath_font_size = 15
imgmath_use_preview = True

# latex
latex_elements = {
    'maketitle':
    r'''
\date{Jul 27, 2021}
\sphinxmaketitle
        ''',
}
latex_documents = [
    (master_doc, 'aphros.tex', u'Aphros documentation',
     author.replace(', ', '\\and ').replace(' and ', '\\and and '), 'manual'),
]

# man
man_pages = [
    (master_doc, 'aphros', u'Aphros documentation', '', 1),
    ('log/index', 'ap.log', u'development log', '', 7),
    ('styleguide', 'ap.styleguide', u'coding style guide', '', 7),
]


for lib in glob.glob('lib?'):
    sec = int(lib[3])
    for p in glob.glob(lib + '/*.rst'):
        if p in [lib + '/index.rst']:
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
    (master_doc, 'aphros', u'Aphros documentation',
     author, 'aphros', 'One line description of project.',
     'Miscellaneous'),
]

