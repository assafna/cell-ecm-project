import plotly.graph_objs as go


def create_plot(_y_array, _names_array, _x_axis_title, _y_axis_title, _title):
    _data = []
    for _y, _name in zip(_y_array, _names_array):
        _data.append(
            go.Box(
                y=_y,
                name=_name,
                boxpoints='all',
                jitter=0.3,
                pointpos=-1.8
            )
        )

    _fig = go.Figure(_data)

    _fig.update_layout(
        title=_title,
        xaxis_title=_x_axis_title,
        yaxis_title=_y_axis_title
    )

    return _fig


def create_group_plot(_x_array, _y_array, _names_array, _x_axis_title, _y_axis_title, _title):
    _data = []
    _x_ticks = []
    for _x, _y, _name in zip(_x_array, _y_array, _names_array):
        _data.append(
            go.Box(
                x=_x,
                y=_y,
                name=_name
            )
        )
        for _x_value in _x:
            _x_ticks.append(_x_value) if _x_value not in _x_ticks else None

    _fig = go.Figure(_data)

    _fig.update_layout(
        title=_title,
        xaxis_title=_x_axis_title,
        xaxis=dict(
            tickmode='array',
            tickvals=sorted(_x_ticks)
        ),
        yaxis_title=_y_axis_title,
        boxmode='group'
    )

    return _fig
