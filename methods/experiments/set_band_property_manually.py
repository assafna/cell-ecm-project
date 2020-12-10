import os

from libs import save_lib
from libs.experiments import load, paths


def main():
    _tuples = [
        ('LiveLive_201207', 1, 'cells_1_3', True),
        ('LiveLive_201207', 1, 'cells_1_4', True),
        ('LiveLive_201207', 1, 'cells_2_3', True),
        ('LiveLive_201207', 1, 'cells_2_4', True),
        ('LiveLive_201207', 1, 'cells_3_4', True),

        ('LiveLive_201207', 3, 'cells_0_1', True),
        ('LiveLive_201207', 3, 'cells_0_2', True),
        ('LiveLive_201207', 3, 'cells_1_2', True),
        ('LiveLive_201207', 3, 'cells_1_3', True),
        ('LiveLive_201207', 3, 'cells_2_3', True),

        ('LiveLive_201207', 4, 'cells_0_1', True),
        ('LiveLive_201207', 4, 'cells_0_2', True),
        ('LiveLive_201207', 4, 'cells_1_2', True),

        ('LiveLive_201207', 5, 'cells_0_1', True),
        ('LiveLive_201207', 5, 'cells_1_4', True),

        ('LiveLive_201207', 6, 'cells_1_3', True),
        ('LiveLive_201207', 6, 'cells_1_4', True),
        ('LiveLive_201207', 6, 'cells_3_4', True),

        ('LiveLive_201207', 7, 'cells_0_1', True),
        ('LiveLive_201207', 7, 'cells_0_2', True),
        ('LiveLive_201207', 7, 'cells_1_2', True),
        ('LiveLive_201207', 7, 'cells_1_3', True),
        ('LiveLive_201207', 7, 'cells_2_3', True),

        ('LiveLive_201207', 8, 'cells_0_1', True),
        ('LiveLive_201207', 8, 'cells_0_3', True),
        ('LiveLive_201207', 8, 'cells_1_3', True),
        ('LiveLive_201207', 8, 'cells_2_3', True),

        ('LiveLive_201207', 9, 'cells_0_2', True),
        ('LiveLive_201207', 9, 'cells_0_3', True),
        ('LiveLive_201207', 9, 'cells_1_2', True),
        ('LiveLive_201207', 9, 'cells_1_3', True),
        ('LiveLive_201207', 9, 'cells_2_3', True),

        ('LiveLive_201207', 10, 'cells_2_6', True),
        ('LiveLive_201207', 10, 'cells_2_7', True),
        ('LiveLive_201207', 10, 'cells_4_6', True),
        ('LiveLive_201207', 10, 'cells_4_7', True),
        ('LiveLive_201207', 10, 'cells_4_8', True),
        ('LiveLive_201207', 10, 'cells_6_7', True),
        ('LiveLive_201207', 10, 'cells_7_8', True),

        ('LiveLive_201207', 12, 'cells_0_1', True),
        ('LiveLive_201207', 12, 'cells_0_2', True),
        ('LiveLive_201207', 12, 'cells_1_2', True),

        ('LiveLive_201207', 13, 'cells_1_3', True),
        ('LiveLive_201207', 13, 'cells_1_4', True),
        ('LiveLive_201207', 13, 'cells_2_3', True),
        ('LiveLive_201207', 13, 'cells_2_4', True),
        ('LiveLive_201207', 13, 'cells_3_4', True),

        ('LiveLive_201207', 14, 'cells_0_1', True),

        ('LiveLive_201207', 15, 'cells_0_1', True),
        ('LiveLive_201207', 15, 'cells_0_2', True),

        ('LiveLive_201207', 16, 'cells_2_4', True)
    ]

    for _tuple in _tuples:
        _e, _s, _g, _b = _tuple
        _p = load.group_properties(_e, _s, _g)
        _p['band'] = _b
        _group_structured_path = paths.structured(_e, _s, _g)
        _properties_json_path = os.path.join(_group_structured_path, 'properties.json')
        save_lib.to_json(_p, _properties_json_path)
        print(_tuple)

        _g_fake = 'fake_' + _g.split('cells_')[1]
        _group_structured_path = paths.structured(_e, _s, _g_fake)
        _properties_json_path = os.path.join(_group_structured_path, 'properties.json')
        if os.path.isfile(_properties_json_path):
            print(_g_fake)
            _p = load.group_properties(_e, _s, _g_fake)
            _p['band'] = _b
            save_lib.to_json(_p, _properties_json_path)


if __name__ == '__main__':
    main()
