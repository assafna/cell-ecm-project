import os

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
    _left_cell_fibers_densities = compute.roi_fibers_density_by_time(
        _simulation=_simulation,
        _length_x=ROI_WIDTH,
        _length_y=ROI_HEIGHT,
        _offset_x=OFFSET_X,
        _offset_y=OFFSET_Y,
        _cell_id='left_cell',
        _direction=_direction,
        _time_points=TIME_POINTS
    )
    _left_cell_fibers_densities_normalized = compute_lib.z_score(
        _left_cell_fibers_densities,
        _normalization['average'],
        _normalization['std']
    )
    _right_cell_fibers_densities = compute.roi_fibers_density_by_time(
        _simulation=_simulation,
        _length_x=ROI_WIDTH,
        _length_y=ROI_HEIGHT,
        _offset_x=OFFSET_X,
        _offset_y=OFFSET_Y,
        _cell_id='right_cell',
        _direction=_direction,
        _time_points=TIME_POINTS
    )
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
    _simulations = filtering.by_distances(_simulations, _distances=[CELLS_DISTANCE])
    _simulations_by_std = organize.by_heterogeneity(_simulations)
    _insides_correlations = []
    _outsides_correlations = []
    _x_array = []
    for _std in _simulations_by_std:
        for _simulation in _simulations_by_std[_std]:
            print(_simulation)
            _insides_correlations.append(run(_simulation, _direction='inside'))
            _outsides_correlations.append(run(_simulation, _direction='outside'))
            _x_array.append(_std)

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
