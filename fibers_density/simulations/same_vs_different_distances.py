import os
from itertools import product
from multiprocessing.pool import Pool

import plotly.graph_objs as go
from scipy.stats import wilcoxon
from tqdm import tqdm

from libs import compute_lib
from libs.config_lib import CPUS_TO_USE
from libs.simulations import compute, filtering, load, paths
from libs.simulations.config import ROI_WIDTH, ROI_HEIGHT
from plotting import save

# by low connectivity
TIME_POINT = {
    False: 50,
    True: 35
}
OFFSET_X = 0
OFFSET_Y = 0
DERIVATIVE = 2
STD = 0.5
CELLS_DISTANCES = [4, 5, 7, 9]


def compute_fibers_densities(_simulations, _low_connectivity):
    _arguments = []
    for _simulation in _simulations:
        for _cell_id in ['left_cell', 'right_cell']:
            _arguments.append({
                'simulation': _simulation,
                'length_x': ROI_WIDTH,
                'length_y': ROI_HEIGHT,
                'offset_x': OFFSET_X,
                'offset_y': OFFSET_Y,
                'cell_id': _cell_id,
                'direction': 'inside',
                'time_points': TIME_POINT[_low_connectivity]
            })

    _fibers_densities = {}
    with Pool(CPUS_TO_USE) as _p:
        for _keys, _value in tqdm(
                _p.imap_unordered(compute.roi_fibers_density_by_time, _arguments),
                total=len(_arguments), desc='Computing Rois & Fibers Densities'):
            _fibers_densities[
                (_keys['simulation'], _keys['cell_id'])] = _value
        _p.close()
        _p.join()

    return _fibers_densities


def main(_low_connectivity=False):
    _simulations = load.structured()
    _simulations = filtering.by_time_points_amount(_simulations, TIME_POINT[_low_connectivity])
    _simulations = filtering.by_categories(
        _simulations,
        _is_single_cell=False,
        _is_heterogeneity=True,
        _is_low_connectivity=_low_connectivity,
        _is_causality=False,
        _is_dominant_passive=False
    )
    _simulations = filtering.by_heterogeneity(_simulations, _std=STD)
    _simulations = filtering.by_distances(_simulations, _distances=CELLS_DISTANCES)
    print('Total simulations:', len(_simulations))

    _fibers_densities = compute_fibers_densities(_simulations, _low_connectivity)

    _y_arrays = [[] for _i in CELLS_DISTANCES]
    for _distance_index, _distance in enumerate(CELLS_DISTANCES):
        _distance_simulations = filtering.by_distance(_simulations, _distance=_distance)
        print('Distance:', _distance, 'Total simulations:', len(_distance_simulations))
        _higher_same_counter = 0
        for _same_index in tqdm(range(len(_distance_simulations)), desc='Main loop'):
            _same_simulation = _distance_simulations[_same_index]
            _same_left_cell_fibers_densities = _fibers_densities[(_same_simulation, 'left_cell')]
            _same_right_cell_fibers_densities = _fibers_densities[(_same_simulation, 'right_cell')]
            _same_correlation = compute_lib.correlation(
                compute_lib.derivative(_same_left_cell_fibers_densities, _n=DERIVATIVE),
                compute_lib.derivative(_same_right_cell_fibers_densities, _n=DERIVATIVE)
            )
            for _different_index in range(len(_distance_simulations)):
                if _same_index != _different_index:
                    _different_simulation = _distance_simulations[_different_index]
                    for _same_cell_id, _different_cell_id in product(['left_cell', 'right_cell'],
                                                                     ['left_cell', 'right_cell']):
                        _same_fibers_densities = \
                            _fibers_densities[(_same_simulation, _same_cell_id)]
                        _different_fibers_densities = \
                            _fibers_densities[(_different_simulation, _different_cell_id)]
                        _different_correlation = compute_lib.correlation(
                            compute_lib.derivative(_same_fibers_densities, _n=DERIVATIVE),
                            compute_lib.derivative(_different_fibers_densities, _n=DERIVATIVE)
                        )
                        _point_distance = compute_lib.distance_from_a_point_to_a_line(
                            _line=[-1, -1, 1, 1],
                            _point=[_same_correlation, _different_correlation]
                        )
                        if _same_correlation > _different_correlation:
                            _y_arrays[_distance_index].append(_point_distance)
                            _higher_same_counter += 1
                        else:
                            _y_arrays[_distance_index].append(-_point_distance)
        print('Total points:', len(_y_arrays[_distance_index]))
        print('Wilcoxon around the zero:')
        print(wilcoxon(_y_arrays[_distance_index]))
        print('Higher same amount:', _higher_same_counter / len(_y_arrays[_distance_index]))

    # plot
    _colors_array = ['#011f4b', '#00417c', '#2e82bf', '#56caed']
    _fig = go.Figure(
        data=[
            go.Box(
                y=_y,
                name=str(_distance),
                boxpoints=False,
                line={
                    'width': 1
                },
                marker={
                    'size': 10,
                    'color': _color
                },
                showlegend=False
            ) for _y, _distance, _color in zip(_y_arrays, CELLS_DISTANCES, _colors_array)
        ],
        layout={
            'xaxis': {
                'title': 'Cells distance (cell diameter)',
                'zeroline': False,
                'tickmode': 'array',
                'tickvals': CELLS_DISTANCES,
                'type': 'category'
            },
            'yaxis': {
                'title': 'Distance from y = x',
                'range': [-1, 1],
                'zeroline': False,
                'tickmode': 'array',
                'tickvals': [-1, -0.5, 0, 0.5, 1]
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
