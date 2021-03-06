import os

import numpy as np
import plotly.graph_objs as go

from fiber_density import inner_density_vs_window_distance_cell_pairs, inner_density_vs_window_distance_single_cells
from fiber_density.inner_density_vs_window_distance_single_cells import OFFSETS_X
from libs import paths_lib
from plotting import save


def main(_low_connectivity=False):
    print('Simulations')
    _simulations_pairs_fiber_densities = \
        inner_density_vs_window_distance_cell_pairs.compute_simulations_data(_low_connectivity)
    _simulations_single_cells_fiber_densities = \
        inner_density_vs_window_distance_single_cells.compute_simulations_data(_low_connectivity)

    _min_simulations_fiber_densities_length = \
        min(len(_simulations_pairs_fiber_densities), len(_simulations_single_cells_fiber_densities))

    _simulations_fiber_densities_differences = \
        np.mean(_simulations_pairs_fiber_densities[:_min_simulations_fiber_densities_length], axis=1) - \
        np.mean(_simulations_single_cells_fiber_densities[:_min_simulations_fiber_densities_length], axis=1)

    print('Experiments')
    _experiments_pairs_fiber_densities, _ = inner_density_vs_window_distance_cell_pairs.compute_experiments_data()
    _experiments_single_cells_fiber_densities = inner_density_vs_window_distance_single_cells.compute_experiments_data()

    _experiments_pairs_fiber_densities_averages = \
        np.array([np.mean(_array) for _array in _experiments_pairs_fiber_densities])
    _experiments_single_cells_fiber_densities_averages = \
        np.array([np.mean(_array) for _array in _experiments_single_cells_fiber_densities])

    _min_experiments_fiber_densities_length = \
        min(len(_experiments_pairs_fiber_densities_averages), len(_experiments_single_cells_fiber_densities_averages))
    _experiments_fiber_densities_differences = \
        _experiments_pairs_fiber_densities_averages[:_min_experiments_fiber_densities_length] - \
        _experiments_single_cells_fiber_densities_averages[:_min_experiments_fiber_densities_length]

    # plot
    _fig = go.Figure(
        data=[
            go.Scatter(
                x=OFFSETS_X,
                y=_simulations_fiber_densities_differences,
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
                y=_experiments_fiber_densities_differences,
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
                'title': 'Window distance (cell diameter)',
                'zeroline': False
            },
            'yaxis': {
                'title': 'Fiber density difference (z-score)',
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
