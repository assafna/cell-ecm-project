import os

import numpy as np
import plotly.graph_objs as go
from scipy.stats import ranksums

from fiber_density import inner_density_vs_time_cell_pairs, inner_density_vs_time_single_cells
from libs.experiments import config, paths
from plotting import save


def main():
    print('Single cells')
    _single_cells_fiber_densities = inner_density_vs_time_single_cells.compute_experiments_fiber_densities()

    print('Cell pairs')
    _cell_pairs_fiber_densities = inner_density_vs_time_cell_pairs.compute_experiments_data()

    # wilcoxon
    print('Wilcoxon rank sums test between single cells and cell pairs for every time point:')
    for _index, _time in enumerate(np.array(range(inner_density_vs_time_cell_pairs.EXPERIMENTS_TIME_FRAMES)) * 15):
        if _index == len(_single_cells_fiber_densities) or _index == len(_cell_pairs_fiber_densities):
            break
        _result = ranksums(_single_cells_fiber_densities[_index], _cell_pairs_fiber_densities[_index])
        print('Time: ', _time)
        print(_result)

    # plot
    _colors = config.colors(2)
    _fig = go.Figure(
        data=[
            go.Scatter(
                x=np.array(range(inner_density_vs_time_single_cells.EXPERIMENTS_TIME_FRAMES)) *
                inner_density_vs_time_single_cells.EXPERIMENTS_TEMPORAL_RESOLUTION,
                y=[np.mean(_array) for _array in _single_cells_fiber_densities],
                name='Single cells',
                error_y={
                    'type': 'data',
                    'array': [np.std(_array) for _array in _single_cells_fiber_densities],
                    'thickness': 1,
                    'color': _colors[0]
                },
                mode='markers',
                marker={
                    'size': 15,
                    'color': _colors[0]
                },
                opacity=0.7
            ),
            go.Scatter(
                x=np.array(range(inner_density_vs_time_cell_pairs.EXPERIMENTS_TIME_FRAMES)) * 15,
                y=[np.mean(_array) for _array in _cell_pairs_fiber_densities],
                name='Cell pairs',
                error_y={
                    'type': 'data',
                    'array': [np.std(_array) for _array in _cell_pairs_fiber_densities],
                    'thickness': 1,
                    'color': _colors[1]
                },
                mode='markers',
                marker={
                    'size': 15,
                    'color': _colors[1]
                },
                opacity=0.7
            )
        ],
        layout={
            'xaxis': {
                'title': 'Time (minutes)',
                'zeroline': False
            },
            'yaxis': {
                'title': 'Fiber density (z-score)',
                'range': [-0.2, 14],
                'zeroline': False,
                'tickmode': 'array',
                'tickvals': [0, 4, 8, 12]
            },
            'legend': {
                'xanchor': 'left',
                'x': 0.1,
                'yanchor': 'top',
                'bordercolor': 'black',
                'borderwidth': 2,
                'bgcolor': 'white'
            },
            'shapes': [
                {
                    'type': 'line',
                    'x0': -15,
                    'y0': 0,
                    'x1': 270,
                    'y1': 0,
                    'line': {
                        'color': 'black',
                        'width': 2
                    }
                },
                {
                    'type': 'line',
                    'x0': -15,
                    'y0': 0,
                    'x1': -15,
                    'y1': 14,
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
        _filename='plot'
    )


if __name__ == '__main__':
    main()
