import plotly.graph_objs as go


def create_plot(_x_labels, _y_labels, _z_array, _x_axis_title, _y_axis_title, _color_scale, _show_scale, _title):
    _fig = go.Figure(
        data=go.Contour(
            x=_x_labels,
            y=_y_labels,
            z=_z_array,
            colorscale=_color_scale,
            showscale=_show_scale,
            zmin=0,
            zmax=1
        )
    )

    _fig.update_layout(
        title=_title,
        xaxis_title=_x_axis_title,
        yaxis_title=_y_axis_title
    )

    return _fig
