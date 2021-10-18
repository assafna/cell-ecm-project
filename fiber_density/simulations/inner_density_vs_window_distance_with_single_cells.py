import os

import numpy as np
import plotly.graph_objs as go

from fiber_density import inner_density_vs_window_distance_single_cells
from fiber_density.simulations import inner_density_vs_window_distance
from libs.simulations import config, paths
from plotting import save


def main(_low_connectivity=False):
    _x_array = []
    _y_array = []
    _names_array = []
    _max_offsets_x = []

    print('Single cells')
    _single_cells_fiber_densities = inner_density_vs_window_distance_single_cells.compute_simulations_data(_low_connectivity)

    print('Cell pairs')
    _names_array, _x_array, _y_array = inner_density_vs_window_distance.compute_cell_pairs(_low_connectivity)

    # plot
    _colors_array = config.colors(4)
    _fig = go.Figure(
        data=[
            go.Scatter(
                x=inner_density_vs_window_distance_single_cells.OFFSETS_X,
                y=[np.mean(_array) for _array in _single_cells_fiber_densities],
                name='Single cells',
                error_y={
                    'type': 'data',
                    'array': [np.std(_array) for _array in _single_cells_fiber_densities],
                    'thickness': 1,
                    'color': _colors_array[0]
                },
                mode='markers',
                marker={
                    'size': 15,
                    'color': _colors_array[0]
                },
                opacity=0.7
            )] + [
            go.Scatter(
                x=_x,
                y=[np.mean(_array) for _array in _y],
                name=_name,
                error_y={
                    'type': 'data',
                    'array': [np.std(_array) for _array in _y],
                    'thickness': 1,
                    'color': _color
                },
                mode='markers',
                marker={
                    'size': 15,
                    'color': _color
                },
                opacity=0.7
            ) for _x, _y, _name, _color in zip(_x_array, _y_array, _names_array, _colors_array[1:])
        ],
        layout={
            'xaxis': {
                'title': 'Window distance (cell diameter)',
                'zeroline': False
            },
            'yaxis': {
                'title': 'Fiber density (z-score)',
                'range': [-1.7, 13],
                'zeroline': False,
                'tickmode': 'array',
                'tickvals': [0, 4, 8, 12]
            },
            'legend': {
                'xanchor': 'right',
                'yanchor': 'top',
                'bordercolor': 'black',
                'borderwidth': 2
            },
            'shapes': [
                {
                    'type': 'line',
                    'x0': -0.2,
                    'y0': -1.5,
                    'x1': 3.4,
                    'y1': -1.5,
                    'line': {
                        'color': 'black',
                        'width': 2
                    }
                },
                {
                    'type': 'line',
                    'x0': -0.2,
                    'y0': -1.5,
                    'x1': -0.2,
                    'y1': 13,
                    'line': {
                        'color': 'black',
                        'width': 2
                    }
                }
            ]
        }
    )

    save.to_html(
        _fig=_fig,
        _path=os.path.join(paths.PLOTS, save.get_module_name()),
        _filename='plot_low_con_' + str(_low_connectivity)
    )


if __name__ == '__main__':
    main()
