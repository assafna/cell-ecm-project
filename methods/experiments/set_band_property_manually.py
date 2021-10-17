import os

from libs import save_lib
from libs.experiments import load, paths


def main():
    _tuples = [
        ('Metastasis_Airyscan_210401_sgAXL', 2, 'cells_1_2', True),
        ('Metastasis_Airyscan_210401_sgAXL', 3, 'cells_1_3', True),
        ('Metastasis_Airyscan_210401_sgAXL', 5, 'cells_0_1', True),
        ('Metastasis_Airyscan_210401_sgAXL', 9, 'cells_0_5', True),
        ('Metastasis_Airyscan_210401_sgAXL', 11, 'cells_1_2', True),
        ('Metastasis_Airyscan_210401_sgAXL', 14, 'cells_0_1', True),
        ('Metastasis_Airyscan_210401_sgAXL', 15, 'cells_1_2', True)
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
