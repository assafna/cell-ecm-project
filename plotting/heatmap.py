import plotly.graph_objs as go


def create_plot(_x_labels, _y_labels, _z_array, _x_axis_title, _y_axis_title, _color_scale=None, _show_scale=True,
                _zmin=None, _zmax=None, _zsmooth=None, _title=None):
    _fig = go.Figure(
        data=go.Heatmap(
            x=_x_labels,
            y=_y_labels,
            z=_z_array,
            colorscale=_color_scale,
            showscale=_show_scale,
            zmin=_zmin,
            zmax=_zmax,
            zsmooth=_zsmooth
        )
    )

    _fig.update_layout(
        title=_title,
        xaxis_title=_x_axis_title,
        yaxis_title=_y_axis_title
    )

    return _fig
