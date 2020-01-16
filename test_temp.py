from libs import save_lib
from libs.experiments import paths, load

if __name__ == '__main__':
    for _experiment in paths.folders(paths.STRUCTURED):
        for _tuple in load.experiment_groups_as_tuples(_experiment):
            _experiment, _series_id, _group = _tuple
            _group_properties = load.group_properties(_experiment, _series_id, _group)
            _group_properties['band'] = None
            _path = paths.group_properties(_experiment, 'Series ' + str(_series_id), _group)
            save_lib.to_json(_group_properties, _path)
