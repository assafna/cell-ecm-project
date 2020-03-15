import os
import time

import numpy as np
import plotly.graph_objs as go

from libs import compute_lib
from libs.simulations import load, filtering, compute, paths
from libs.simulations.config import ROI_WIDTH, ROI_HEIGHT, CELL_DIAMETER
from plotting import scatter, save

TIME_POINT = 50
OFFSET_X_START = 0
OFFSET_X_END = CELL_DIAMETER * 4
OFFSET_X_STEP = CELL_DIAMETER / 8
OFFSET_Y = 0


def main():
    _simulations = load.structured()
    _simulations = filtering.by_categories(
        _simulations,
        _is_single_cell=True,
        _is_heterogeneity=False,
        _is_low_connectivity=False,
        _is_causality=False,
        _is_dominant_passive=False
    )
    _simulations = filtering.by_time_points_amount(_simulations, _time_points=TIME_POINT)
    _fibers_densities = list([None] * len(np.arange(start=OFFSET_X_START, stop=OFFSET_X_END, step=OFFSET_X_STEP)))
    for _simulation in _simulations:
        print(_simulation)
        _offset_index = 0
        _normalization = load.normalization(_simulation)
        for _offset_x in np.arange(start=OFFSET_X_START, stop=OFFSET_X_END, step=OFFSET_X_STEP):
            _fibers_density_right = compute.roi_fibers_density_time_point({
                'simulation': _simulation,
                'length_x': ROI_WIDTH,
                'length_y': ROI_HEIGHT,
                'offset_x': _offset_x,
                'offset_y': OFFSET_Y,
                'cell_id': 'cell',
                'direction': 'right',
                'time_point': TIME_POINT
            })[1]
            _normalized_fibers_density_right = compute_lib.z_score(
                _fibers_density_right,
                _normalization['average'],
                _normalization['std']
            )
            _fibers_density_left = compute.roi_fibers_density_time_point({
                'simulation': _simulation,
                'length_x': ROI_WIDTH,
                'length_y': ROI_HEIGHT,
                'offset_x': _offset_x,
                'offset_y': OFFSET_Y,
                'cell_id': 'cell',
                'direction': 'left',
                'time_point': TIME_POINT
            })[1]
            _normalized_fibers_density_left = compute_lib.z_score(
                _fibers_density_left,
                _normalization['average'],
                _normalization['std']
            )
            _fibers_density_up = compute.roi_fibers_density_time_point({
                'simulation': _simulation,
                'length_x': ROI_HEIGHT,
                'length_y': ROI_WIDTH,
                'offset_x': OFFSET_Y,
                'offset_y': _offset_x,
                'cell_id': 'cell',
                'direction': 'up',
                'time_point': TIME_POINT
            })[1]
            _normalized_fibers_density_up = compute_lib.z_score(
                _fibers_density_up,
                _normalization['average'],
                _normalization['std']
            )
            _fibers_density_down = compute.roi_fibers_density_time_point({
                'simulation': _simulation,
                'length_x': ROI_HEIGHT,
                'length_y': ROI_WIDTH,
                'offset_x': OFFSET_Y,
                'offset_y': _offset_x,
                'cell_id': 'cell',
                'direction': 'down',
                'time_point': TIME_POINT
            })[1]
            _normalized_fibers_density_down = compute_lib.z_score(
                _fibers_density_down,
                _normalization['average'],
                _normalization['std']
            )
            _fibers_densities_average = \
                (_normalized_fibers_density_right +
                 _normalized_fibers_density_left +
                 _normalized_fibers_density_up +
                 _normalized_fibers_density_down) / 4
            if _fibers_densities[_offset_index] is None:
                _fibers_densities[_offset_index] = [_fibers_densities_average]
            else:
                _fibers_densities[_offset_index].append(_fibers_densities_average)
            _offset_index += 1

    # plot
    _fig = go.Figure(
        data=[
            go.Scatter(
                x=np.arange(start=OFFSET_X_START, stop=OFFSET_X_END, step=OFFSET_X_STEP) / CELL_DIAMETER,
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
                'title': 'Distance from Left Cell (cell size)'
            },
            'yaxis': {
                'title': 'Fibers Density Z-score'
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
