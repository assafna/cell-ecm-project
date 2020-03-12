import os

import numpy as np

from fibers_density import fibers_vs_distance_pairs, fibers_vs_distance_single_cells
from libs import paths_lib
from plotting import scatter, save

OFFSET_X_STEP = 0.2
OFFSETS_X = np.arange(start=0, stop=100, step=OFFSET_X_STEP)


def main():
    print('Experiments')
    _experiments_pairs_fibers_densities = fibers_vs_distance_pairs.compute_experiments_data()
    _experiments_single_cells_fibers_densities = fibers_vs_distance_single_cells.compute_experiments_data()
    _minimum_length = min(len(_experiments_pairs_fibers_densities), len(_experiments_single_cells_fibers_densities))
    _experiments_fibers_densities_differences = \
        np.mean(_experiments_pairs_fibers_densities, axis=1)[:_minimum_length] - \
        np.mean(_experiments_single_cells_fibers_densities, axis=1)[:_minimum_length]

    print('Simulations')
    _simulations_pairs_fibers_densities = fibers_vs_distance_pairs.compute_simulations_data()
    _simulations_single_cells_fibers_densities = fibers_vs_distance_single_cells.compute_simulations_data()
    _minimum_length = min(len(_simulations_pairs_fibers_densities), len(_simulations_single_cells_fibers_densities))
    _simulations_fibers_densities_differences = \
        np.mean(_simulations_pairs_fibers_densities, axis=1)[:_minimum_length] - \
        np.mean(_simulations_single_cells_fibers_densities, axis=1)[:_minimum_length]

    # plot
    _fig = scatter.create_plot(
        _x_array=[OFFSETS_X] * 2,
        _y_array=[_experiments_fibers_densities_differences, _simulations_fibers_densities_differences],
        _names_array=['Experiments', 'Simulations'],
        _modes_array=['lines+markers'] * 2,
        _show_legend_array=[True] * 2,
        _x_axis_title='Distance from Left Cell (cell size)',
        _y_axis_title='Fibers Density Z-score Difference'
    )

    save.to_html(
        _fig=_fig,
        _path=os.path.join(paths_lib.PLOTS, save.get_module_name()),
        _filename='plot'
    )


if __name__ == '__main__':
    main()
