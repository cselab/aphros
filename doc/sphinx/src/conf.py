extensions = [
              'sphinx.ext.githubpages',
              'sphinx.ext.imgmath',
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

html_static_path = ['_static']

html_context = {'css_files': ['_static/center.css']}

html_theme_options = {
      'nosidebar' : True,
      'body_max_width' : '80%',
      }

html_sidebars = {
   '**' : ['localtoc.html']
   #'**' : []
}

htmlhelp_basename = 'hydrodoc'

html_math_renderer = 'imgmath'
#html_math_renderer = ''

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
