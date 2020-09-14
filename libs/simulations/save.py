import os

from libs.save_lib import to_pickle
from libs.simulations import paths


def normalization(_simulation, _normalization):
    _normalization_path = os.path.join(paths.NORMALIZATION, _simulation + '.pkl')
    to_pickle(_normalization, _normalization_path)
