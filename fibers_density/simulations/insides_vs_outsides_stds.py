import os

from scipy.stats import wilcoxon

from libs import compute_lib
from libs.simulations import load, filtering, organize, compute, paths
from libs.simulations.config import ROI_WIDTH, ROI_HEIGHT
from plotting import box, save

TIME_POINTS = 50
OFFSET_X = 0
OFFSET_Y = 0
DERIVATIVE = 2
CELLS_DISTANCE = 5.0


def run(_simulation, _direction):
    _normalization = load.normalization(_simulation)
    _left_cell_fibers_densities = compute.roi_fibers_density_by_time({
        'simulation': _simulation,
        'length_x': ROI_WIDTH,
        'length_y': ROI_HEIGHT,
        'offset_x': OFFSET_X,
        'offset_y': OFFSET_Y,
        'cell_id': 'left_cell',
        'direction': _direction,
        'time_points': TIME_POINTS
    })
    _left_cell_fibers_densities_normalized = compute_lib.z_score(
        _left_cell_fibers_densities,
        _normalization['average'],
        _normalization['std']
    )
    _right_cell_fibers_densities = compute.roi_fibers_density_by_time({
        'simulation': _simulation,
        'length_x': ROI_WIDTH,
        'length_y': ROI_HEIGHT,
        'offset_x': OFFSET_X,
        'offset_y': OFFSET_Y,
        'cell_id': 'right_cell',
        'direction': _direction,
        'time_points': TIME_POINTS
    })
    _right_cell_fibers_densities_normalized = compute_lib.z_score(
        _right_cell_fibers_densities,
        _normalization['average'],
        _normalization['std']
    )

    return compute_lib.correlation(
        compute_lib.derivative(_left_cell_fibers_densities_normalized, _n=DERIVATIVE),
        compute_lib.derivative(_right_cell_fibers_densities_normalized, _n=DERIVATIVE)
    )


def main():
    _simulations = load.structured()
    _simulations = filtering.by_time_points_amount(_simulations, _time_points=TIME_POINTS)
    _simulations = filtering.by_categories(
        _simulations,
        _is_single_cell=False,
        _is_heterogeneity=None,
        _is_low_connectivity=False,
        _is_causality=False,
        _is_dominant_passive=False
    )
    _simulations = filtering.by_distance(_simulations, _distance=CELLS_DISTANCE)
    _simulations_by_std = organize.by_heterogeneity(_simulations)
    _insides_correlations = []
    _outsides_correlations = []
    _x_array = []
    for _std in _simulations_by_std:
        _std_insides_correlations = []
        _std_outsides_correlations = []
        for _simulation in _simulations_by_std[_std]:
            print(_simulation)
            _insides_correlation = run(_simulation, _direction='inside')
            _outsides_correlation = run(_simulation, _direction='outside')
            _insides_correlations.append(_insides_correlation)
            _outsides_correlations.append(_outsides_correlation)
            _std_insides_correlations.append(_insides_correlation)
            _std_outsides_correlations.append(_outsides_correlation)
            _x_array.append(_std)

        # wilcoxon
        print('Wilcoxon STD', _std, 'insides', wilcoxon(_std_insides_correlations), sep='\t')
        print('Wilcoxon STD', _std, 'outsides', wilcoxon(_std_outsides_correlations), sep='\t')


    # plot
    _fig = box.create_group_plot(
        _x_array=[_x_array] * 2,
        _y_array=[_insides_correlations, _outsides_correlations],
        _names_array=['Insides', 'Outsides'],
        _x_axis_title='STDs',
        _y_axis_title='Correlation',
        _title='Insides & Outsides Correlation by Cell Contractions STDs'
    )

    save.to_html(
        _fig=_fig,
        _path=os.path.join(paths.PLOTS, save.get_module_name()),
        _filename='plot'
    )


if __name__ == '__main__':
    main()
