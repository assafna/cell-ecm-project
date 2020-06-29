import os

import numpy as np
import plotly.graph_objs as go

from libs import compute_lib
from libs.simulations import load, filtering, compute, paths
from libs.simulations.config import CELL_DIAMETER, ROI_WIDTH, ROI_HEIGHT
from plotting import scatter, save

TIME_POINT = 50
CELLS_DISTANCE = 7.0
OFFSET_X_END = CELLS_DISTANCE * CELL_DIAMETER - (CELL_DIAMETER / 2) - ROI_WIDTH
OFFSET_X_STEP = CELL_DIAMETER / 8
OFFSET_Y = 0


def main():
    _simulations = load.structured()
    _simulations = filtering.by_categories(
        _simulations,
        _is_single_cell=False,
        _is_heterogeneity=False,
        _is_low_connectivity=False,
        _is_causality=False,
        _is_dominant_passive=False
    )
    _simulations = filtering.by_distance(_simulations, _distance=CELLS_DISTANCE)
    _simulations = filtering.by_time_points_amount(_simulations, _time_points=TIME_POINT)
    _fibers_densities = list([None] * len(np.arange(start=0, stop=OFFSET_X_END, step=OFFSET_X_STEP)))
    for _simulation in _simulations:
        print(_simulation)
        _offset_index = 0
        _normalization = load.normalization(_simulation)
        for _offset_x in np.arange(start=0, stop=OFFSET_X_END, step=OFFSET_X_STEP):
            _fibers_density = compute.roi_fibers_density_time_point({
                'simulation': _simulation,
                'length_x': ROI_WIDTH,
                'length_y': ROI_HEIGHT,
                'offset_x': _offset_x,
                'offset_y': OFFSET_Y,
                'cell_id': 'left_cell',
                'direction': 'right',
                'time_point': TIME_POINT
            })[1]
            _normalized_fibers_density = compute_lib.z_score(
                _fibers_density,
                _normalization['average'],
                _normalization['std']
            )

            if _fibers_densities[_offset_index] is None:
                _fibers_densities[_offset_index] = [_normalized_fibers_density]
            else:
                _fibers_densities[_offset_index].append(_normalized_fibers_density)
            _offset_index += 1

    # plot
    _fig = go.Figure(
        data=[
            go.Scatter(
                x=np.arange(start=0, stop=OFFSET_X_END, step=OFFSET_X_STEP) / CELL_DIAMETER,
                y=[np.mean(_array) for _array in _fibers_densities],
                name='Time-Point Last',
                error_y={
                    'type': 'data',
                    'array': [np.std(_array) for _array in _fibers_densities],
                    'thickness': 1
                },
                mode='lines+markers',
                line={'dash': 'solid'}
            )
        ],
        layout={
            'xaxis': {
                'title': 'Distance from left cell (cell size)'
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
