def update_y_axis(_fig, _range, _color=None, _width=None):
    _fig.update_layout(
        yaxis=dict(
            zerolinecolor=_color,
            zerolinewidth=_width,
            range=_range
        )
    )

    return _fig


def update_x_axis(_fig, _range, _color=None, _width=None):
    _fig.update_layout(
        xaxis=dict(
            zerolinecolor=_color,
            zerolinewidth=_width,
            range=_range
        )
    )

    return _fig
