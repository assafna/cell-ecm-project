import numpy as np
import plotly.graph_objs as go


def create_error_bars_plot(_x_array, _y_array, _names_array, _modes_array, _dashes_array, _x_axis_title, _y_axis_title,
                           _title=None):
    _fig = go.Figure()

    for _x, _y, _name, _mode, _dash in zip(_x_array, _y_array, _names_array, _modes_array, _dashes_array):
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


def create_plot(_x_array, _y_array, _names_array, _modes_array, _show_legend_array, _x_axis_title, _y_axis_title,
                _title=None):
    _data = []
    for _x, _y, _name, _mode, _show_legend in zip(_x_array, _y_array, _names_array, _modes_array, _show_legend_array):
        _data.append(
            go.Scatter(
                x=_x,
                y=_y,
                name=_name,
                mode=_mode,
                showlegend=_show_legend
            )
        )

    _fig = go.Figure(_data)

    _fig.update_layout(
        title=_title,
        xaxis_title=_x_axis_title,
        yaxis_title=_y_axis_title
    )

    return _fig


def add_line(_fig, _x1, _y1, _x2, _y2, _name, _color, _showlegend):
    _fig.add_trace(
        go.Scatter(
            x=[_x1, _x2],
            y=[_y1, _y2],
            name=_name,
            marker=dict(
                color=_color
            ),
            mode='lines',
            showlegend=_showlegend
        )
    )

    return _fig
