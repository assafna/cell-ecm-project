import os
import sys

from scipy import stats

from libs import compute_lib
from libs.simulations import filtering, load, organize, compute, paths
from libs.simulations.config import ROI_WIDTH, ROI_HEIGHT
from plotting import scatter, save

TIME_POINTS = sys.maxsize
OFFSET_X = 0
OFFSET_Y = 0
DERIVATIVE = 1
DISTANCES = [3.0, 7.0, 12.0]
DIRECTION = 'inside'
HETEROGENEITY = False


def main():
    _simulations = load.structured()
    _simulations = filtering.by_categories(
        _simulations,
        _is_single_cell=False,
        _is_heterogeneity=HETEROGENEITY,
        _is_low_connectivity=False,
        _is_causality=False,
        _is_dominant_passive=False
    )
    # _simulations = filtering.by_distances(_simulations, DISTANCES)
    _simulations_by_distance = organize.by_distances(_simulations)
    _fibers_densities_by_distance = {}
    _change_in_fibers_densities_by_distance = {}
    for _distance in _simulations_by_distance:
        _distance_simulations = _simulations_by_distance[_distance]
        _fibers_densities_array = []
        _change_in_fibers_densities_array = []
        for _simulation in _distance_simulations:
            print(_simulation)
            _simulation_normalization = load.normalization(_simulation)
            _simulation_normalization = [_simulation_normalization['average'], _simulation_normalization['std']]
            for _cell_id in ['left_cell', 'right_cell']:
                _cell_fibers_densities = compute.roi_fibers_density_by_time({
                    'simulation': _simulation,
                    'length_x': ROI_WIDTH,
                    'length_y': ROI_HEIGHT,
                    'offset_x': OFFSET_X,
                    'offset_y': OFFSET_Y,
                    'cell_id': _cell_id,
                    'direction': DIRECTION,
                    'time_points': TIME_POINTS
                })
                _z_score_fibers_density = compute_lib.z_score_fibers_densities_array(
                    _cell_fibers_densities, _simulation_normalization
                )
                _fibers_densities_array += _z_score_fibers_density
                _change_in_fibers_densities_array += [0] + compute_lib.derivative(
                    _z_score_fibers_density, _n=DERIVATIVE
                )
        _fibers_densities_by_distance[_distance] = _fibers_densities_array
        _change_in_fibers_densities_by_distance[_distance] = _change_in_fibers_densities_array

    _fig = scatter.create_plot(
        _x_array=_fibers_densities_by_distance.values(),
        _y_array=_change_in_fibers_densities_by_distance.values(),
        _names_array=['Distance ' + str(_distance) for _distance in _fibers_densities_by_distance.keys()],
        _modes_array=['markers'] * len(_fibers_densities_by_distance.keys()),
        _show_legend_array=[True] * len(_fibers_densities_by_distance.keys()),
        _x_axis_title='Fibers Densities Z-Score',
        _y_axis_title='Change in Fibers Densities Z-Score',
        _title='Fibers Densities vs. Change in Fibers Densities - ' + DIRECTION.capitalize()
    )

    save.to_html(
        _fig=_fig,
        _path=os.path.join(paths.PLOTS, save.get_module_name()),
        _filename=DIRECTION + '_heterogeneity_' + str(HETEROGENEITY) + '_points'
    )

    # line of best fit
    _best_fit_lines_x_array = []
    _best_fit_lines_y_array = []
    for _distance in _simulations_by_distance:
        _x_array = _fibers_densities_by_distance[_distance]
        _y_array = _change_in_fibers_densities_by_distance[_distance]
        _slope, _intercept, _r_value, _p_value, _std_err = stats.linregress(_x_array, _y_array)
        _x1, _x2 = min(_x_array), max(_x_array)
        _y1, _y2 = _slope * _x1 + _intercept, _slope * _x2 + _intercept
        _best_fit_lines_x_array.append([_x1, _x2])
        _best_fit_lines_y_array.append([_y1, _y2])

    _fig = scatter.create_plot(
        _x_array=_best_fit_lines_x_array,
        _y_array=_best_fit_lines_y_array,
        _names_array=['Distance ' + str(_distance) for _distance in _simulations_by_distance],
        _modes_array=['lines'] * len(_simulations_by_distance),
        _show_legend_array=[True] * len(_simulations_by_distance),
        _x_axis_title='Fibers Densities Z-Score',
        _y_axis_title='Change in Fibers Densities Z-Score',
        _title='Fibers Densities vs. Change in Fibers Densities - ' + DIRECTION.capitalize() + ' - Line of Best Fit'
    )

    save.to_html(
        _fig=_fig,
        _path=os.path.join(paths.PLOTS, save.get_module_name()),
        _filename=DIRECTION + '_heterogeneity_' + str(HETEROGENEITY) + '_best_fit'
    )


if __name__ == '__main__':
    main()
