import plotly.graph_objs as go


def create_plot(_x_labels, _y_labels, _z_array, _x_axis_title, _y_axis_title, _color_scale, _show_scale, _title):
    _fig = go.Figure(
        data=go.Heatmap(
            x=_x_labels,
            y=_y_labels,
            z=_z_array,
            colorscale=_color_scale,
            # colorbar=dict(
            #     tickformat=',.0%'
            # ),
            showscale=_show_scale,
            zmin=0,
            zmax=20
        )
    )

    _fig.update_layout(
        title=_title,
        xaxis_title=_x_axis_title,
        yaxis_title=_y_axis_title
    )

    return _fig
