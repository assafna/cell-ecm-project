from libs.experiments import load, save


def main():
    _tuples = [
        ('LiveDead_201115', 1, [0, 1]),
        ('LiveDead_201115', 2, [1]),
        ('LiveDead_201115', 3, [0]),
        ('LiveDead_201115', 4, [0, 3]),
        ('LiveDead_201115', 5, []),
        ('LiveDead_201115', 6, [0]),
        ('LiveDead_201115', 7, []),
        ('LiveDead_201115', 8, [1]),
        ('LiveDead_201115', 9, [0, 1, 2]),
        ('LiveDead_201115', 10, [0, 1]),
        ('LiveDead_201115', 11, [0, 1]),
        ('LiveDead_201115', 12, [1, 2]),
        ('LiveDead_201115', 13, [0, 1]),
        ('LiveDead_201115', 14, [1, 2]),
        ('LiveDead_201115', 15, [2]),
        ('LiveDead_201115', 16, [0]),
        ('LiveDead_201115', 17, [0])
    ]

    for _tuple in _tuples:
        _experiment, _series_id, _dead_cell_ids = _tuple
        _properties = load.image_properties(_experiment, _series_id)
        _properties['dead_cell_ids'] = _dead_cell_ids
        save.image_properties(_experiment, _series_id, _properties)
        print(_tuple)


if __name__ == '__main__':
    main()
