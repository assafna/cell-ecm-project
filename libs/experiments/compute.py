import os

import numpy as np

from libs.experiments import paths
from libs.experiments.load import file_data, serieses


def normalized_fibers_density(_fibers_density, _normalization):
    return {k: list(np.array(_fibers_density[k]) / _normalization - 1) for k in _fibers_density.keys()}


def fibers_density_cut_edges(_fibers_density, _cut_amount=4):
    return {k: list(_fibers_density[k][_cut_amount:-_cut_amount]) for k in _fibers_density.keys()}


def fibers_density_cut_left_edge(_fibers_density, _cut_amount=4):
    return {k: list(_fibers_density[k][_cut_amount:]) for k in _fibers_density.keys()}


def group(_series_path, _group):
    _group_path = os.path.join(_series_path, _group)
    return file_data(_group_path)[0]


def series(_experiment_path, _series):
    _series_path = os.path.join(_experiment_path, _series)
    _series_groups = [_group for _group in os.listdir(_series_path) if _group.startswith('cell')]
    _series_data = {}
    for _group in _series_groups:
        _series_data[str(_group.split('.')[0])] = group(_series_path, _group)

    return _series_data


def experiment(_experiment):
    _experiment_path = paths.normalization(_experiment)
    _experiment_serieses = serieses(_experiment_path)
    _experiment_data = {}
    for _series in _experiment_serieses:
        _experiment_data[_series] = series(_experiment_path, _series)

    return _experiment_data