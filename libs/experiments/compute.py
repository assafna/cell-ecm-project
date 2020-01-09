import math

from libs.experiments import load
from libs.experiments.config import CELL_DIAMETER_IN_MICRONS


def z_score(_x, _average, _std):
    return (_x - _average) / _std


def z_score_fibers_density_array(_fibers_density, _normalization):
    _average, _std = _normalization
    return {k: [z_score(_x, _average, _std) for _x in _fibers_density[k]] for k in _fibers_density.keys()}


def fibers_density_cut_edges(_fibers_density, _cut_amount=4):
    return {k: list(_fibers_density[k][_cut_amount:-_cut_amount]) for k in _fibers_density.keys()}


def fibers_density_cut_left_edge(_fibers_density, _cut_amount=4):
    return {k: list(_fibers_density[k][_cut_amount:]) for k in _fibers_density.keys()}


def cells_distance_in_cell_size(_experiment, _series, _cell_1_coordinates, _cell_2_coordinates):
    _series_id = _series.split()[1]
    _image_properties = load.image_properties(_experiment, _series)
    _image_resolutions = _image_properties['resolutions']
    _x1, _y1, _z1 = [float(_value) for _value in _cell_1_coordinates[0]]
    _x2, _y2, _z2 = [float(_value) for _value in _cell_2_coordinates[0]]
    _x1, _y1, _z1 = _x1 * _image_resolutions['x'], _y1 * _image_resolutions['y'], _z1 * _image_resolutions['z']
    _x2, _y2, _z2 = _x2 * _image_resolutions['x'], _y2 * _image_resolutions['y'], _z2 * _image_resolutions['z']

    return math.sqrt((_x1 - _x2) ** 2 + (_y1 - _y2) ** 2 + (_z1 - _z2) ** 2) / CELL_DIAMETER_IN_MICRONS
