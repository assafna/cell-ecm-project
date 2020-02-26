import os
import sys
from multiprocessing.pool import Pool

from matplotlib import pyplot as plt

from libs import compute_lib
from libs.config_lib import CPUS_TO_USE
from libs.experiments import compute, paths, load

ROI_LENGTH = 1
ROI_HEIGHT = 1
ROI_WIDTH = 1
OFFSET_X = 0
OFFSET_Y = 0
OFFSET_Z = 0
DIRECTION = 'inside'
TIME_POINTS = sys.maxsize
DERIVATIVE = 2


def process_group(_experiment, _series_id, _group):
    _arguments = {
        'experiment': _experiment,
        'series_id': _series_id,
        'group': _group,
        'length_x': ROI_LENGTH,
        'length_y': ROI_HEIGHT,
        'length_z': ROI_WIDTH,
        'offset_x': OFFSET_X,
        'offset_y': OFFSET_Y,
        'offset_z': OFFSET_Z,
        'direction': DIRECTION,
        'time_points': TIME_POINTS
    }

    _fibers_densities = compute.roi_fibers_density_by_time_pairs(_arguments)

    # remove blacklist
    _fibers_densities = {
        'left_cell': compute.remove_blacklist(_experiment, _series_id, 'left_cell', _fibers_densities['left_cell']),
        'right_cell': compute.remove_blacklist(_experiment, _series_id, 'left_cell', _fibers_densities['right_cell'])
    }

    # get longest sequence
    _fibers_densities = compute.longest_same_indices_shared_in_borders_sub_array(
        _fibers_densities['left_cell'],
        _fibers_densities['right_cell']
    )
    _fibers_densities = {
        'left_cell': _fibers_densities[0],
        'right_cell': _fibers_densities[1]
    }

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

    # save
    _path = paths.images(_experiment + ' - In Time', 'Series ' + str(_series_id))
    os.makedirs(_path, exist_ok=True)
    plt.savefig(os.path.join(_path, _group + '.png'))
    print('Saved:', _experiment, _series_id, _group, sep='\t')


def process_experiment(_experiment):
    _arguments = []
    for _tuple in load.experiment_groups_as_tuples(_experiment):
        _experiment, _series_id, _group = _tuple
        _arguments.append((_experiment, _series_id, _group))

    _p = Pool(CPUS_TO_USE)
    _p.starmap(process_group, _arguments)
    _p.close()


if __name__ == '__main__':
    process_experiment('SN41')
    # process_group('SN41', 8, 'cells_1_2')
