import os
from multiprocessing.pool import Pool

import plotly.graph_objs as go
from scipy.stats import pearsonr
from tqdm import tqdm

from libs import compute_lib
from libs.config_lib import CPUS_TO_USE
from libs.simulations import compute, filtering, load, paths
from libs.simulations.config import ROI_WIDTH, ROI_HEIGHT
from plotting import save

TIME_POINT = 50
OFFSET_X = 0
OFFSET_Y = 0
DERIVATIVE = 2
STD = 0.5


def compute_fibers_densities(_simulations):
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
                'time_points': TIME_POINT
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


def main():
    _simulations = load.structured()
    _simulations = filtering.by_time_points_amount(_simulations, TIME_POINT)
    _simulations = filtering.by_categories(
        _simulations,
        _is_single_cell=False,
        _is_heterogeneity=True,
        _is_low_connectivity=False,
        _is_causality=False,
        _is_dominant_passive=False
    )
    _simulations = filtering.by_heterogeneity(_simulations, _std=STD)
    print('Total simulations:', len(_simulations))

    _fibers_densities = compute_fibers_densities(_simulations)

    _correlations = []
    _pair_distances = []
    for _simulation in tqdm(_simulations, desc='Simulations loop'):
        _properties = load.properties(_simulation)

        _left_cell_fibers_densities = _fibers_densities[(_simulation, 'left_cell')]
        _right_cell_fibers_densities = _fibers_densities[(_simulation, 'right_cell')]

        _correlations.append(compute_lib.correlation(
            compute_lib.derivative(_left_cell_fibers_densities, _n=DERIVATIVE),
            compute_lib.derivative(_right_cell_fibers_densities, _n=DERIVATIVE)
        ))
        _pair_distances.append(compute.cells_distance(_properties))

    print('Total pairs:', len(_correlations))
    print('Pearson:', pearsonr(_correlations, _pair_distances))

    # plot
    _fig = go.Figure(
        data=go.Scatter(
            x=_pair_distances,
            y=_correlations,
            mode='markers',
            marker={
                'size': 10,
                'color': 'black'
            },
            showlegend=False
        ),
        layout={
            'xaxis': {
                'title': 'Pair distance (cell diameter)',
                'zeroline': False,
                # 'tickmode': 'array',
                # 'tickvals': [4, 6, 8, 10]
            },
            'yaxis': {
                'title': 'Correlation',
                'zeroline': False,
                'range': [-1.1, 1.2],
                'tickmode': 'array',
                'tickvals': [-1, -0.5, 0, 0.5, 1]
            },
            # 'shapes': [
            #     {
            #         'type': 'line',
            #         'x0': -1,
            #         'y0': -1,
            #         'x1': -1,
            #         'y1': 1,
            #         'line': {
            #             'color': 'black',
            #             'width': 2
            #         }
            #     },
            #     {
            #         'type': 'line',
            #         'x0': -1,
            #         'y0': -1,
            #         'x1': 1,
            #         'y1': -1,
            #         'line': {
            #             'color': 'black',
            #             'width': 2
            #         }
            #     }
            # ]
        }
    )

    save.to_html(
        _fig=_fig,
        _path=os.path.join(paths.PLOTS, save.get_module_name()),
        _filename='plot'
    )


if __name__ == '__main__':
    main()
