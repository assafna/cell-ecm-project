import gzip
import os
import pickle

import paths


def to_pickle(_object, _pickle_path):
    os.makedirs(os.path.dirname(_pickle_path)) if not os.path.isdir(os.path.dirname(_pickle_path)) else None
    try:
        with gzip.open(_pickle_path, 'wb') as _pickle:
            pickle.dump(_object, _pickle, -1)
    finally:
        _pickle.close()


def fibers_densities(_simulation, _time_point, _fibers_densities):
    _simulation_fibers_densities_path = paths.fibers_densities(_simulation)
    _time_point_fibers_densities_path = os.path.join(_simulation_fibers_densities_path, str(_time_point) + '.pkl')
    to_pickle(_fibers_densities, _time_point_fibers_densities_path)
