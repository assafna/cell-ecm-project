import os

import plotly


def div_to_website(_div, _filepath):
    _web_page = '''
    <html>
        <head>
            <script src="https://cdn.plot.ly/plotly-latest.min.js"></script> 
        </head>
        <body>
    ''' + _div + '''
        </body>
    </html>
    '''

    try:
        with open(_filepath + '.html', 'w') as _html:
            _html.write(_web_page)
    finally:
        _html.close()


def to_html(_fig, _path, _filename):
    os.makedirs(_path, exist_ok=True)

    _div = plotly.offline.plot(
        figure_or_data=_fig,
        output_type='div',
        include_plotlyjs=False
    )

    div_to_website(
        _div=_div,
        _filepath=os.path.join(_path, _filename)
    )
