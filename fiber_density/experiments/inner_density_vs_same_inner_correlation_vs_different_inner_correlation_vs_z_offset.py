import os

import numpy as np
import plotly.graph_objs as go

from fiber_density.experiments import inner_density_vs_offsets_in_axes, \
    same_inner_correlation_vs_different_inner_correlation_offsets_in_axes
from fiber_density.experiments.inner_density_vs_offsets_in_axes import OFFSET_Y_START, OFFSET_Y_END, OFFSET_Y_STEP
from libs import compute_lib
from libs.experiments import paths
from plotting import save

OFFSET_Z = 0


def main():
    print('Computing fiber density')
    _fiber_densities = inner_density_vs_offsets_in_axes.compute_z_array(
        _offset_z_start=OFFSET_Z, _offset_z_end=OFFSET_Z, _offset_z_step=OFFSET_Z + 1).flatten()

    print('Computing "same vs. different"')
    _communications = same_inner_correlation_vs_different_inner_correlation_offsets_in_axes.compute_z_array(
        _offset_z_start=OFFSET_Z, _offset_z_end=OFFSET_Z, _offset_z_step=OFFSET_Z + 1).flatten()

    print('Correlation:', compute_lib.correlation(_fiber_densities, _communications, _with_p_value=True))

    # plot
    _x = np.arange(start=OFFSET_Y_START, stop=OFFSET_Y_END + OFFSET_Y_STEP, step=OFFSET_Y_STEP)
    _fig = go.Figure(
        data=[
            go.Scatter(
                x=_x,
                y=_fiber_densities,
                name='Fiber density',
                mode='markers+lines',
                marker={
                    'size': 10,
                    'color': '#844b00'
                },
                line={
                    'width': 1
                },
                showlegend=True
            ),
            go.Scatter(
                x=_x,
                y=_communications,
                yaxis='y2',
                name='Communication',
                mode='markers+lines',
                marker={
                    'size': 10,
                    'color': '#edbc80'
                },
                line={
                    'width': 1
                },
                showlegend=True
            )
        ],
        layout={
            'xaxis': {
                'title': 'Offset in Z axis (cell diameter)',
                'zeroline': False
            },
            'yaxis': {
                'title': 'Fiber density (z-score)',
                'titlefont': {
                    'color': '#844b00'
                },
                'tickfont': {
                    'color': '#844b00'
                },
                'zeroline': False
            },
            'yaxis2': {
                'title': '"same" > "different" (%)',
                'titlefont': {
                    'color': '#edbc80'
                },
                'tickfont': {
                    'color': '#edbc80'
                },
                'overlaying': 'y',
                'side': 'right',
                'showgrid': False,
                'zeroline': False,
                'tickangle': -90
            },
            'legend': {
                'xanchor': 'left',
                'x': 0.1,
                'yanchor': 'top',
                'bordercolor': 'black',
                'borderwidth': 2,
                'bgcolor': 'white'
            }
        }
    )

    save.to_html(
        _fig=_fig,
        _path=os.path.join(paths.PLOTS, save.get_module_name()),
        _filename='plot'
    )


if __name__ == '__main__':
    main()
