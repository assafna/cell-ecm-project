import inspect
import os
import time

import plotly


def replace_image_download(_div):
    return _div.replace(
        '{"responsive": true}',
        '{modeBarButtonsToAdd: [{ name: "SVG", '
        'icon: Plotly.Icons.camera, click: function(gd) {Plotly.downloadImage(gd, {format: "svg", '
        'height: window.innerHeight, width: window.innerWidth})}}]},{"responsive": true} '
    )


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


def to_html(_fig, _path, _filename, _and_to_image=True):
    # update theme
    # options: plotly_white, presentation, none
    _fig.update_layout(
        template='presentation',
        xaxis={
            'showgrid': False,
            'zerolinewidth': 2,
            'zerolinecolor': 'black'
        },
        yaxis={
            'showgrid': False,
            'zerolinewidth': 2,
            'zerolinecolor': 'black',
            'tickangle': -90
        },
        font={
            'family': 'Arial',
            'size': 25,
            'color': 'black'
        },
        paper_bgcolor='rgba(0,0,0,0)',
        plot_bgcolor='rgba(0,0,0,0)'
    )

    # create div
    _div = plotly.offline.plot(
        figure_or_data=_fig,
        output_type='div',
        include_plotlyjs=False
    )

    # update download button
    _div = replace_image_download(_div)

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
        _filepath=os.path.join(_backup_path, _backup_filename)
    )

    print('Saved plot:', _filename)

    if _and_to_image:
        to_image(_fig, _path, _filename)


def to_image(_fig, _path, _filename, _format='svg', _width=500, _height=500):
    os.makedirs(_path, exist_ok=True)
    _fig.write_image(
        os.path.join(_path, _filename + '.' + _format),
        width=_width,
        height=_height
    )

    print('Saved image:', _filename)
