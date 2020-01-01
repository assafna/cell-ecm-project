import numpy as np
import plotly.graph_objs as go


def create_error_bars_plot(_x_array, _y_array, _names_array, _mode_array, _dash_array, _x_axis_title, _y_axis_title, _title):
    _fig = go.Figure()

    for _x, _y, _name, _mode, _dash in zip(_x_array, _y_array, _names_array, _mode_array, _dash_array):
        _averages = []
        _std = []
        for _array in _y:
            _averages.append(np.mean(_array))
            _std.append(np.std(_array))
        _fig.add_trace(
            go.Scatter(
                y=_averages,
                x=_x,
                name=_name,
                error_y=dict(
                    type='data',
                    array=_std,
                    thickness=1
                ),
                mode=_mode,
                line=dict(
                    dash=_dash
                )
            )
        )

    _fig.update_layout(
        title=_title,
        xaxis_title=_x_axis_title,
        yaxis_title=_y_axis_title
    )

    return _fig


def create_plot(_x_array, _y_array, _names_array, _modes_array, _x_axis_title, _y_axis_title, _title):
    _data = []
    for _x, _y, _name, _mode in zip(_x_array, _y_array, _names_array, _modes_array):
        _data.append(
            go.Scatter(
                x=_x,
                y=_y,
                name=_name,
                mode=_mode
            )
        )

    _fig = go.Figure(_data)

    _fig.update_layout(
        title=_title,
        xaxis_title=_x_axis_title,
        yaxis_title=_y_axis_title
    )

    return _fig


def update_y_axis(_fig, _color, _width, _format, _range):
    _fig.update_layout(
        yaxis=dict(
            zerolinecolor=_color,
            zerolinewidth=_width,
            tickformat=_format,
            range=_range
        )
    )

    return _fig
