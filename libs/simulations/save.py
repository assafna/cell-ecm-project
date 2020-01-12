import os

from libs.save_lib import to_pickle
from libs.simulations import paths


def fibers_densities(_simulation, _time_point, _fibers_densities):
    _simulation_fibers_densities_path = paths.fibers_densities(_simulation)
    _time_point_fibers_densities_path = os.path.join(_simulation_fibers_densities_path, str(_time_point) + '.pkl')
    to_pickle(_fibers_densities, _time_point_fibers_densities_path)


def normalization(_simulation, _normalization):
    _normalization_path = os.path.join(paths.NORMALIZATION, _simulation + '.pkl')
    to_pickle(_normalization, _normalization_path)
