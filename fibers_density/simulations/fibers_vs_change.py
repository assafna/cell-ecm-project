from libs import compute
from libs.simulations import compute, filtering
from libs.simulations.config import ROI_WIDTH, ROI_HEIGHT


TIME_POINTS = 51
OFFSET_X = 0
OFFSET_Y = 0
DERIVATIVE = 1
CELLS_DISTANCE = 3.0
DIRECTION = 'inside'


def main():
    _simulations = ['3D_1', '3D_2', '3D_3']
    _simulations = filtering.by_distance(_simulations, CELLS_DISTANCE)
    _fibers_densities_array = []
    _change_in_fibers_densities_array = []
    for _simulation in _simulations:
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
            _change_in_fibers_densities_array += [0] + compute.derivative(_cell_fibers_densities, _n=DERIVATIVE)

    # TODO: create the plots
    # TODO: create one plot for all distances


if __name__ == '__main__':
    main()
