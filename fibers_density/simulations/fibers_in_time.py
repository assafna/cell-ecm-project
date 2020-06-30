import os

import numpy as np
import plotly.graph_objs as go

from libs import compute_lib
from libs.simulations import load, compute, filtering, paths
from libs.simulations.config import CELL_DIAMETER, ROI_WIDTH, ROI_HEIGHT
from plotting import scatter, save

TIME_POINTS = 50
CELLS_DISTANCE = 5.0
OFFSET_X = 0
OFFSET_Y = 0


def run(_simulations, _cell_id, _directions):
    _fibers_densities = list([None] * TIME_POINTS)
    for _simulation in _simulations:
        print(_simulation)
        _time_point_index = 0
        _normalization = load.normalization(_simulation)
        for _time_point in range(TIME_POINTS):
            _direction_fibers_densities = []
            for _direction in _directions:
                _fibers_density = compute.roi_fibers_density_time_point({
                    'simulation': _simulation,
                    'length_x': ROI_HEIGHT if _direction in ['up', 'down'] else ROI_WIDTH,
                    'length_y': ROI_WIDTH if _direction in ['up', 'down'] else ROI_HEIGHT,
                    'offset_x': OFFSET_Y if _direction in ['up', 'down'] else OFFSET_X,
                    'offset_y': OFFSET_X if _direction in ['up', 'down'] else OFFSET_Y,
                    'cell_id': _cell_id,
                    'direction': _direction,
                    'time_point': _time_point
                })[1]
                _normalized_fibers_density = compute_lib.z_score(
                    _fibers_density,
                    _normalization['average'],
                    _normalization['std']
                )
                _direction_fibers_densities.append(_normalized_fibers_density)

            if _fibers_densities[_time_point_index] is None:
                _fibers_densities[_time_point_index] = [np.mean(_direction_fibers_densities)]
            else:
                _fibers_densities[_time_point_index].append(np.mean(_direction_fibers_densities))
            _time_point_index += 1

    return _fibers_densities


def main():
    _simulations = load.structured()
    _simulations = filtering.by_time_points_amount(_simulations, _time_points=TIME_POINTS)

    # single cell
    _simulations_single_cells = filtering.by_categories(
        _simulations,
        _is_single_cell=True,
        _is_heterogeneity=False,
        _is_low_connectivity=False,
        _is_causality=False,
        _is_dominant_passive=False
    )
    _fibers_densities_single_cells = run(
        _simulations=_simulations_single_cells,
        _cell_id='cell',
        _directions=['left', 'right', 'up', 'down']
    )

    # pairs
    _simulations = filtering.by_categories(
        _simulations,
        _is_single_cell=False,
        _is_heterogeneity=False,
        _is_low_connectivity=False,
        _is_causality=False,
        _is_dominant_passive=False
    )
    _simulations_pairs = filtering.by_distance(_simulations, _distance=CELLS_DISTANCE)
    _fibers_densities_pairs_left = run(
        _simulations=_simulations_pairs,
        _cell_id='left_cell',
        _directions=['inside']
    )
    _fibers_densities_pairs_right = run(
        _simulations=_simulations_pairs,
        _cell_id='right_cell',
        _directions=['inside']
    )

    # combine left and right cells
    _fibers_densities_pairs = []
    for _left, _right in zip(_fibers_densities_pairs_left, _fibers_densities_pairs_right):
        _fibers_densities_pairs.append(_left + _right)

    # plot
    _fig = go.Figure(
        data=[
            go.Scatter(
                x=list(range(TIME_POINTS)),
                y=[np.mean(_array) for _array in _fibers_densities_single_cells],
                name='Single Cells',
                error_y={
                    'type': 'data',
                    'array': [np.std(_array) for _array in _fibers_densities_single_cells],
                    'thickness': 1
                },
                mode='lines+markers',
                line={'dash': 'dash'}
            ),
            go.Scatter(
                x=list(range(TIME_POINTS)),
                y=[np.mean(_array) for _array in _fibers_densities_pairs],
                name='Pairs',
                error_y={
                    'type': 'data',
                    'array': [np.std(_array) for _array in _fibers_densities_pairs],
                    'thickness': 1
                },
                mode='lines+markers',
                line={'dash': 'solid'}
            )
        ],
        layout={
            'xaxis': {
                'title': 'Time (cell contraction percentages)'
            },
            'yaxis': {
                'title': 'Fiber density (z-score)'
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
