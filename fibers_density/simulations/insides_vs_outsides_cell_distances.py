import os

from libs import compute_lib
from libs.simulations import load, filtering, organize, compute, paths
from libs.simulations.config import ROI_WIDTH, ROI_HEIGHT, CELL_DIAMETER
from plotting import box, save, edit

TIME_POINTS = 50
OFFSET_X = CELL_DIAMETER * 1
OFFSET_Y = 0
DERIVATIVE = 2
STD = 0.5
CELLS_DISTANCES = [5.0, 7.0, 9.0, 12.0]


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
        _is_heterogeneity=True,
        _is_low_connectivity=False,
        _is_causality=False,
        _is_dominant_passive=False
    )
    _simulations = filtering.by_distances(_simulations, _distances=CELLS_DISTANCES)
    _simulations = filtering.by_heterogeneity(_simulations, _std=STD)
    _simulations_by_distances = organize.by_distances(_simulations)
    _insides_correlations = []
    _outsides_correlations = []
    _x_array = []
    for _distance in _simulations_by_distances:
        for _simulation in _simulations_by_distances[_distance]:
            print(_simulation)
            _insides_correlations.append(run(_simulation, _direction='inside'))
            _outsides_correlations.append(run(_simulation, _direction='outside'))
            _x_array.append(_distance)

    # plot
    _fig = box.create_group_plot(
        _x_array=[_x_array] * 2,
        _y_array=[_insides_correlations, _outsides_correlations],
        _names_array=['Insides', 'Outsides'],
        _x_axis_title='Cells Distances',
        _y_axis_title='Correlation',
        _title='Insides & Outsides Correlation by Cell Distance'
    )

    _fig = edit.update_y_axis(
        _fig=_fig,
        _range=[-0.5, 1.0]
    )

    save.to_html(
        _fig=_fig,
        _path=os.path.join(paths.PLOTS, save.get_module_name()),
        _filename='plot_offset_x_1_cell'
    )


if __name__ == '__main__':
    main()
