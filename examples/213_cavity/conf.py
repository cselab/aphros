import os
import sys

sys.path.append(os.path.abspath('..'))
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

linkpath_github_root = 'https://github.com/cselab/aphros/blob/master/'
linkpath_link_github = True

project = 'Lid-driven cavity'

language = 'en'

pygments_style = 'sphinx'

numfig = True

# html
html_theme = 'theme'
html_theme_path = ['.']
htmlhelp_basename = 'aphrosdoc'

html_math_renderer = 'imgmath'
imgmath_image_format = 'svg'
imgmath_font_size = 15
imgmath_use_preview = True
