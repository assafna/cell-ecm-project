import os

import numpy as np
import plotly.graph_objs as go

from fibers_density import fibers_vs_distance_pairs, fibers_vs_distance_single_cells
from fibers_density.fibers_vs_distance_single_cells import OFFSETS_X, OFFSET_X_STEP
from libs import paths_lib
from plotting import save


def main(_low_connectivity=False):
    print('Simulations')
    _simulations_pairs_fibers_densities = \
        fibers_vs_distance_pairs.compute_simulations_data(_low_connectivity)
    _simulations_single_cells_fibers_densities = \
        fibers_vs_distance_single_cells.compute_simulations_data(_low_connectivity)

    _min_simulations_fibers_densities_length = \
        min(len(_simulations_pairs_fibers_densities), len(_simulations_single_cells_fibers_densities))

    _simulations_fibers_densities_differences = \
        np.mean(_simulations_pairs_fibers_densities[:_min_simulations_fibers_densities_length], axis=1) - \
        np.mean(_simulations_single_cells_fibers_densities[:_min_simulations_fibers_densities_length], axis=1)

    print('Experiments')
    _experiments_pairs_fibers_densities, _ = fibers_vs_distance_pairs.compute_experiments_data()
    _experiments_single_cells_fibers_densities = fibers_vs_distance_single_cells.compute_experiments_data()

    _experiments_pairs_fibers_densities_averages = \
        np.array([np.mean(_array) for _array in _experiments_pairs_fibers_densities])
    _experiments_single_cells_fibers_densities_averages = \
        np.array([np.mean(_array) for _array in _experiments_single_cells_fibers_densities])

    _min_experiments_fibers_densities_length = \
        min(len(_experiments_pairs_fibers_densities_averages), len(_experiments_single_cells_fibers_densities_averages))
    _experiments_fibers_densities_differences = \
        _experiments_pairs_fibers_densities_averages[:_min_experiments_fibers_densities_length] - \
        _experiments_single_cells_fibers_densities_averages[:_min_experiments_fibers_densities_length]

    # plot
    _fig = go.Figure(
        data=[
            go.Scatter(
                x=OFFSETS_X,
                y=_simulations_fibers_densities_differences,
                name='Simulations',
                mode='markers',
                marker={
                    'size': 15,
                    'color': '#005b96'
                },
                opacity=0.7
            ),
            go.Scatter(
                x=OFFSETS_X,
                y=_experiments_fibers_densities_differences,
                name='Experiments',
                mode='markers',
                marker={
                    'size': 15,
                    'color': '#ea8500'
                },
                opacity=0.7
            )
        ],
        layout={
            'xaxis': {
                'title': 'Distance from cell (cell diameter)',
                'zeroline': False
            },
            'yaxis': {
                'title': 'Fibers density z-score difference',
                'range': [-0.2, 6],
                'zeroline': False,
                'tickvals': [0, 2, 4]
            },
            'legend': {
                'xanchor': 'right',
                'yanchor': 'top',
                'bordercolor': 'black',
                'borderwidth': 2,
                'bgcolor': 'white'
            },
            'shapes': [
                {
                    'type': 'line',
                    'x0': -0.2,
                    'y0': 0,
                    'x1': 2.6,
                    'y1': 0,
                    'line': {
                        'color': 'black',
                        'width': 2
                    }
                },
                {
                    'type': 'line',
                    'x0': -0.2,
                    'y0': 0,
                    'x1': -0.2,
                    'y1': 6,
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
