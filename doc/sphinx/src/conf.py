extensions = ['sphinx.ext.mathjax',
              'sphinx.ext.githubpages',
              'sphinx.ext.viewcode',
              'includecode',
              ]

templates_path = []

source_suffix = '.rst'

master_doc = 'index'


project = 'hydro'
copyright = '2019, ETH Zurich'
author = 'ETH Zurich'

# short X.Y version
version = "0.1"
# full version
release = "0.1"

language = None

exclude_patterns = []

pygments_style = 'sphinx'

numfig = True


# html
html_theme = 'classic'

html_theme_options = {}

html_sidebars = {
   #'**' : ['localtoc.html']
   '**' : []
}

htmlhelp_basename = 'hydrodoc'

# latex
latex_elements = {}
latex_documents = [
    (master_doc, 'hydro.tex', u'cubism-hydro Documentation',
     u'ETH Zurich', 'manual'),
]


# man
man_pages = [
    (master_doc, 'hydro', u'hydro Documentation',
     [author], 1)
]

# texinfo
texinfo_documents = [
    (master_doc, 'hydro', u'hydro Documentation',
     author, 'hydro', 'One line description of project.',
     'Miscellaneous'),
]
