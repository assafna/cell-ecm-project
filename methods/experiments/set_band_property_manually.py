import os

from libs import save_lib
from libs.experiments import load, paths


def main():
    _tuples = [
        ('SN45', 1, 'cells_0_1', False),
        ('SN45', 1, 'cells_0_2', False),
        ('SN45', 1, 'cells_0_3', False),
        ('SN45', 1, 'cells_0_4', False),
        ('SN45', 1, 'cells_0_5', False),
        ('SN45', 1, 'cells_0_6', False),
        ('SN45', 1, 'cells_1_2', False),
        ('SN45', 1, 'cells_1_3', False),
        ('SN45', 1, 'cells_1_4', False),
        ('SN45', 1, 'cells_1_5', True),
        ('SN45', 1, 'cells_1_6', True),
        ('SN45', 1, 'cells_2_3', False),
        ('SN45', 1, 'cells_2_4', False),
        ('SN45', 1, 'cells_2_5', False),
        ('SN45', 1, 'cells_2_6', False),
        ('SN45', 1, 'cells_3_4', True),
        ('SN45', 1, 'cells_3_5', True),
        ('SN45', 1, 'cells_3_6', False),
        ('SN45', 1, 'cells_4_5', True),
        ('SN45', 1, 'cells_4_6', False),
        ('SN45', 1, 'cells_5_6', False),
        ('SN45', 2, 'cells_0_1', False),
        ('SN45', 2, 'cells_0_2', False),
        ('SN45', 2, 'cells_0_3', False),
        ('SN45', 2, 'cells_1_2', False),
        ('SN45', 2, 'cells_1_3', False),
        ('SN45', 2, 'cells_2_3', True),
        ('SN45', 3, 'cells_0_1', False),
        ('SN45', 3, 'cells_0_2', False),
        ('SN45', 3, 'cells_1_2', True),
        ('SN45', 4, 'cells_0_1', False),
        ('SN45', 4, 'cells_0_2', True),
        ('SN45', 4, 'cells_0_3', False),
        ('SN45', 4, 'cells_1_2', False),
        ('SN45', 4, 'cells_1_3', False),
        ('SN45', 4, 'cells_2_3', False),
        ('SN45', 5, 'cells_0_1', False),
        ('SN45', 5, 'cells_0_2', False),
        ('SN45', 5, 'cells_0_3', False),
        ('SN45', 5, 'cells_0_4', False),
        ('SN45', 5, 'cells_0_5', False),
        ('SN45', 5, 'cells_0_6', False),
        ('SN45', 5, 'cells_1_2', False),
        ('SN45', 5, 'cells_1_3', False),
        ('SN45', 5, 'cells_1_4', False),
        ('SN45', 5, 'cells_1_5', False),
        ('SN45', 5, 'cells_1_6', False),
        ('SN45', 5, 'cells_2_3', False),
        ('SN45', 5, 'cells_2_4', False),
        ('SN45', 5, 'cells_2_5', False),
        ('SN45', 5, 'cells_2_6', False),
        ('SN45', 5, 'cells_3_4', False),
        ('SN45', 5, 'cells_3_5', False),
        ('SN45', 5, 'cells_3_6', False),
        ('SN45', 5, 'cells_4_5', False),
        ('SN45', 5, 'cells_4_6', False),
        ('SN45', 5, 'cells_5_6', False),
        ('SN45', 6, 'cells_0_1', False),
        ('SN45', 6, 'cells_0_2', False),
        ('SN45', 6, 'cells_0_3', False),
        ('SN45', 6, 'cells_1_2', False),
        ('SN45', 6, 'cells_1_3', False),
        ('SN45', 6, 'cells_2_3', False),
        ('SN45', 7, 'cells_0_1', False),
        ('SN45', 8, 'cells_0_1', False),
        ('SN45', 8, 'cells_0_2', False),
        ('SN45', 8, 'cells_0_3', False),
        ('SN45', 8, 'cells_1_2', False),
        ('SN45', 8, 'cells_1_3', False),
        ('SN45', 8, 'cells_2_3', False),
        ('SN45', 9, 'cells_0_1', False),
        ('SN45', 9, 'cells_0_2', False),
        ('SN45', 9, 'cells_1_2', False),
    ]

    for _tuple in _tuples:
        _e, _s, _g, _b = _tuple
        _p = load.group_properties(_e, _s, _g)
        _p['band'] = _b
        _group_structured_path = paths.structured(_e, 'Series ' + str(_s), _g)
        _properties_json_path = os.path.join(_group_structured_path, 'properties.json')
        save_lib.to_json(_p, _properties_json_path)
        print(_tuple)

        _g_fake = 'fake_' + _g.split('cells_')[1]
        _group_structured_path = paths.structured(_e, 'Series ' + str(_s), _g_fake)
        _properties_json_path = os.path.join(_group_structured_path, 'properties.json')
        if os.path.isfile(_properties_json_path):
            print(_g_fake)
            _p = load.group_properties(_e, _s, _g_fake)
            _p['band'] = _b
            save_lib.to_json(_p, _properties_json_path)


if __name__ == '__main__':
    main()
