import numpy as np


def z_score(_x, _average, _std):
    return (_x - _average) / _std


def z_score_fibers_density_array(_fibers_density, _normalization):
    _average, _std = _normalization
    return {k: [z_score(_x, _average, _std) for _x in _fibers_density[k]] for k in _fibers_density.keys()}


def fibers_density_cut_edges(_fibers_density, _cut_amount=4):
    return {k: list(_fibers_density[k][_cut_amount:-_cut_amount]) for k in _fibers_density.keys()}


def fibers_density_cut_left_edge(_fibers_density, _cut_amount=4):
    return {k: list(_fibers_density[k][_cut_amount:]) for k in _fibers_density.keys()}
