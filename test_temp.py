from libs import save_lib
from libs.experiments import paths, load

if __name__ == '__main__':
    for _experiment in ['SN16']:
        for _tuple in load.experiment_groups_as_tuples(_experiment):
            _experiment, _series_id, _group = _tuple
            _image_properties = load.image_properties(_experiment, 'Series ' + str(_series_id))
            _image_properties['experiment'] = 'SN16'
            _path = paths.image_properties(_experiment, 'series_' + str(_series_id) + '.json')
            save_lib.to_json(_image_properties, _path)
