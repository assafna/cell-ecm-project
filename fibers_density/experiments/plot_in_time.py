import sys

from matplotlib import pyplot as plt

from libs import compute_lib
from libs.experiments import compute

ROI_LENGTH = 1
ROI_HEIGHT = 1
ROI_WIDTH = 1
OFFSET_X = 0
OFFSET_Y = 0.5
OFFSET_Z = 0
DIRECTION = 'inside'
TIME_POINTS = sys.maxsize
DERIVATIVE = 0
OUT_OF_BOUNDARIES = True


def process_group(_experiment, _series_id, _group):
    _fibers_densities = compute.roi_fibers_density_by_time_pairs(
        _experiment=_experiment,
        _series_id=_series_id,
        _group=_group,
        _length_x=ROI_LENGTH,
        _length_y=ROI_HEIGHT,
        _length_z=ROI_WIDTH,
        _offset_x=OFFSET_X,
        _offset_y=OFFSET_Y,
        _offset_z=OFFSET_Z,
        _direction=DIRECTION,
        _time_points=TIME_POINTS,
        _out_of_borders=OUT_OF_BOUNDARIES
    )

    # _fibers_densities = compute.longest_same_indices_shared_in_borders_sub_array(
    #     _fibers_densities['left_cell'],
    #     _fibers_densities['right_cell']
    # )
    # _fibers_densities = {
    #     'left_cell': _fibers_densities[0],
    #     'right_cell': _fibers_densities[1]
    # }

    # remove out of boundaries
    if OUT_OF_BOUNDARIES:
        _fibers_densities = {
            'left_cell': [_fibers_density[0] for _fibers_density in _fibers_densities['left_cell']],
            'right_cell': [_fibers_density[0] for _fibers_density in _fibers_densities['right_cell']]
        }

    # delete problematic points
    # del(_fibers_densities['left_cell'][80])
    # del (_fibers_densities['right_cell'][80])
    # del(_fibers_densities['left_cell'][138])
    # del (_fibers_densities['right_cell'][138])

    # derivative
    if DERIVATIVE > 0:
        _fibers_densities = {
            'left_cell': compute_lib.derivative(_fibers_densities['left_cell'], _n=DERIVATIVE),
            'right_cell': compute_lib.derivative(_fibers_densities['right_cell'], _n=DERIVATIVE)
        }

    # correlation
    _correlation = compute_lib.correlation(_fibers_densities['left_cell'], _fibers_densities['right_cell'])

    # plot
    plt.plot(_fibers_densities['left_cell'])
    plt.plot(_fibers_densities['right_cell'])
    plt.legend(['left cell', 'right cell'])
    plt.title('Correlation: ' + str(_correlation))
    plt.show()


if __name__ == '__main__':
    process_group(_experiment='SN41', _series_id=2, _group='cells_0_1')
