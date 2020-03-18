import os

import numpy as np
import plotly.graph_objs as go

from fibers_density import fibers_vs_distance_pairs, fibers_vs_distance_single_cells
from fibers_density.fibers_vs_distance_pairs import OFFSETS_X
from libs import paths_lib
from plotting import scatter, save

OFFSET_X_STEP = 0.2


def main(_low_connectivity=False):
    print('Experiments')
    _experiments_pairs_fibers_densities = fibers_vs_distance_pairs.compute_experiments_data()
    _experiments_single_cells_fibers_densities = fibers_vs_distance_single_cells.compute_experiments_data()

    _experiments_pairs_fibers_densities_averages = \
        np.array([np.mean(_array) for _array in _experiments_pairs_fibers_densities])
    _experiments_single_cells_fibers_densities_averages = \
        np.array([np.mean(_array) for _array in _experiments_single_cells_fibers_densities])

    _experiments_fibers_densities_differences = \
        _experiments_pairs_fibers_densities_averages - _experiments_single_cells_fibers_densities_averages

    print('Simulations')
    _simulations_pairs_fibers_densities = \
        fibers_vs_distance_pairs.compute_simulations_data(_low_connectivity)
    _simulations_single_cells_fibers_densities = \
        fibers_vs_distance_single_cells.compute_simulations_data(_low_connectivity)
    _simulations_fibers_densities_differences = \
        np.mean(_simulations_pairs_fibers_densities, axis=1) - \
        np.mean(_simulations_single_cells_fibers_densities, axis=1)

    # plot
    _fig = go.Figure(
        data=[
            go.Scatter(
                x=OFFSETS_X,
                y=_simulations_fibers_densities_differences,
                name='Simulations',
                mode='markers',
                marker={
                    'size': 15
                },
                opacity=0.7
            ),
            go.Scatter(
                x=OFFSETS_X,
                y=_experiments_fibers_densities_differences,
                name='Experiments',
                mode='markers',
                marker={
                    'size': 15
                },
                opacity=0.7
            )
        ],
        layout={
            'xaxis': {
                'title': 'Distance from Cell (cell size)',
                'zeroline': False
            },
            'yaxis': {
                'title': 'Fibers Density Z-score Difference',
                'range': [-0.2, 3.5],
                'zeroline': False
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
                    'x0': -OFFSET_X_STEP,
                    'y0': 0,
                    'x1': OFFSETS_X[-1] + OFFSET_X_STEP,
                    'y1': 0,
                    'line': {
                        'color': 'black',
                        'width': 2
                    }
                },
                {
                    'type': 'line',
                    'x0': -OFFSET_X_STEP,
                    'y0': 0,
                    'x1': -OFFSET_X_STEP,
                    'y1': 3.5,
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
        _path=os.path.join(paths_lib.PLOTS, save.get_module_name()),
        _filename='plot_low_con_' + str(_low_connectivity)
    )


if __name__ == '__main__':
    main()
