from libs import compute_lib
from libs.simulations import filtering, load, organize, compute
from libs.simulations.config import ROI_WIDTH, ROI_HEIGHT


TIME_POINTS = 51
OFFSET_X = 0
OFFSET_Y = 0
DERIVATIVE = 1
DIRECTION = 'inside'


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
    _simulations_by_distance = organize.by_distances(_simulations)
    _fibers_densities_by_distance = {}
    _change_in_fibers_densities_by_distance = {}
    for _distance in _simulations_by_distance:
        _distance_simulations = _simulations_by_distance[_distance]
        _fibers_densities_array = []
        _change_in_fibers_densities_array = []
        for _simulation in _distance_simulations:
            print(_simulation)
            for _cell_id in ['left_cell', 'right_cell']:
                _cell_fibers_densities = compute.roi_fibers_density_by_time(
                    _simulation=_simulation,
                    _length_x=ROI_WIDTH,
                    _length_y=ROI_HEIGHT,
                    _offset_x=OFFSET_X,
                    _offset_y=OFFSET_Y,
                    _cell_id=_cell_id,
                    _direction=DIRECTION,
                    _time_points=TIME_POINTS
                )
                _fibers_densities_array += _cell_fibers_densities
                _change_in_fibers_densities_array += [0] + compute_lib.derivative(_cell_fibers_densities, _n=DERIVATIVE)
        _fibers_densities_by_distance[_distance] = _fibers_densities_array
        _change_in_fibers_densities_by_distance[_distance] = _change_in_fibers_densities_array

    # TODO: create the plots


if __name__ == '__main__':
    main()
