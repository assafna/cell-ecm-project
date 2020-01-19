import os

import inspect
import time

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


def get_module_name():
    return os.path.splitext(os.path.basename(inspect.stack()[1].filename))[0]


def to_html(_fig, _path, _filename):
    _div = plotly.offline.plot(
        figure_or_data=_fig,
        output_type='div',
        include_plotlyjs=False
    )

    # regular save
    os.makedirs(_path, exist_ok=True)
    div_to_website(
        _div=_div,
        _filepath=os.path.join(_path, _filename)
    )

    # backup save
    _time = time.strftime('%Y_%m_%d-%H_%M_%S')
    _backup_filename = _filename + '-' + _time
    _backup_path = os.path.join(_path, 'History')
    os.makedirs(_backup_path, exist_ok=True)
    div_to_website(
        _div=_div,
        _filepath=os.path.join(_backup_path, _filename)
    )

    print('Saved plot:', _filename)
