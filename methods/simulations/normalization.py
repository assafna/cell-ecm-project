import os
import time

import numpy as np

from libs.simulations import paths, compute, load, save
from libs.simulations.config import ROI_WIDTH, ROI_HEIGHT, CELL_DIAMETER

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

    _simulation_structured_path = paths.structured(_simulation)
    os.makedirs(_simulation_structured_path) if not os.path.isdir(_simulation_structured_path) else None

    _fibers_densities = []
    _offsets_x = list(np.arange(OFFSET_X_START, OFFSET_X_END - ROI_WIDTH, STEP))
    _offsets_y = list(np.arange(OFFSET_Y_START, CELLS_ZONE_Y_START - ROI_HEIGHT, STEP)) +\
                 list(np.arange(CELLS_ZONE_Y_END, OFFSET_Y_END - ROI_HEIGHT, STEP))
    _time_point_fibers_densities = load.fibers_densities(_simulation, _time_point=0)

    _rois = []
    for _offset_x in _offsets_x:
        for _offset_y in _offsets_y:
            _rois.append((_offset_x, _offset_y, _offset_x + ROI_WIDTH, _offset_y + ROI_HEIGHT))

    _rois_sums = compute.rois_fibers_densities(_simulation, _time_point=0, _rois=_rois)

    if len(_time_point_fibers_densities) != len(_rois_sums):
        _time_point_fibers_densities.update(_rois_sums)
        save.fibers_densities(_simulation, _time_point=0, _fibers_densities=_time_point_fibers_densities)


if __name__ == '__main__':
    process_simulation('3D_1')
