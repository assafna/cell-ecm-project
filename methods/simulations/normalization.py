import os
import time
from multiprocessing.pool import Pool

import numpy as np
from tqdm import tqdm

from libs.config_lib import CPUS_TO_USE
from libs.simulations import paths, compute, load, save
from libs.simulations.config import QUANTIFICATION_WINDOW_WIDTH_IN_CELL_DIAMETER, \
    QUANTIFICATION_WINDOW_HEIGHT_IN_CELL_DIAMETER, CELL_DIAMETER

BORDER = 0.2
OFFSET_X_START = -BORDER
OFFSET_X_END = BORDER
OFFSET_Y_START = -BORDER
OFFSET_Y_END = BORDER
STEP = 0.02
CELLS_ZONE_Y_START = -CELL_DIAMETER / 2
CELLS_ZONE_Y_END = CELL_DIAMETER / 2


def process_simulation(_simulation, _overwrite=False):
    _start_time = time.time()

    _normalization_pickle_path = os.path.join(paths.NORMALIZATION, _simulation + '.pkl')
    if not _overwrite and os.path.isfile(_normalization_pickle_path):
        return

    _fiber_density = []
    _offsets_x = list(np.arange(
        OFFSET_X_START, OFFSET_X_END - QUANTIFICATION_WINDOW_WIDTH_IN_CELL_DIAMETER * CELL_DIAMETER, STEP))
    _offsets_y = list(np.arange(
        OFFSET_Y_START, CELLS_ZONE_Y_START - QUANTIFICATION_WINDOW_HEIGHT_IN_CELL_DIAMETER * CELL_DIAMETER, STEP)) + \
        list(np.arange(
            CELLS_ZONE_Y_END, OFFSET_Y_END - QUANTIFICATION_WINDOW_HEIGHT_IN_CELL_DIAMETER * CELL_DIAMETER, STEP))

    # prepare windows
    _windows = []
    for _offset_x in _offsets_x:
        for _offset_y in _offsets_y:
            _windows.append((_offset_x, _offset_y, _offset_x + QUANTIFICATION_WINDOW_WIDTH_IN_CELL_DIAMETER *
                             CELL_DIAMETER, _offset_y + QUANTIFICATION_WINDOW_HEIGHT_IN_CELL_DIAMETER * CELL_DIAMETER))

    # compute
    _windows_sums = compute.quantification_windows_fiber_densities(_simulation, _time_point=0, _windows=_windows)

    # save
    _sums = list(_windows_sums.values())
    _normalization = {
        'average': np.mean(_sums),
        'std': np.std(_sums)
    }
    save.normalization(_simulation, _normalization)


def process_simulations(_simulations, _overwrite=False):
    _arguments = [_simulation for _simulation in _simulations]
    with Pool(CPUS_TO_USE) as _p:
        for _ in tqdm(_p.imap_unordered(process_simulation, _arguments), total=len(_arguments), desc='Simulations'):
            continue
        _p.close()
        _p.join()


def process_all_simulations(_overwrite=False):
    process_simulations(load.structured(), _overwrite)


if __name__ == '__main__':
    process_all_simulations()
