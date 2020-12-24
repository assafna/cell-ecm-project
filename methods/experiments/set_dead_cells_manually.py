from libs.experiments import load, save


def main():
    _tuples = [
        ('LiveDead_201220', 1, [1]),
        ('LiveDead_201220', 2, [1]),
        ('LiveDead_201220', 3, [0, 1]),
        ('LiveDead_201220', 4, [0, 1, 2]),
        ('LiveDead_201220', 5, [1, 2]),
        ('LiveDead_201220', 6, [0]),
        ('LiveDead_201220', 7, [0, 2]),
        ('LiveDead_201220', 8, [0, 1]),
        ('LiveDead_201220', 9, [0, 1]),
        ('LiveDead_201220', 10, [0]),
        ('LiveDead_201220', 11, [1, 2, 3]),
        ('LiveDead_201220', 12, [0, 1]),
        ('LiveDead_201220', 13, [1, 2, 3, 4]),
        ('LiveDead_201220', 14, [0, 2, 3]),
        ('LiveDead_201220', 15, [0, 1]),
        ('LiveDead_201220', 16, [0, 1, 2, 3]),
        ('LiveDead_201220', 17, [0]),
        ('LiveDead_201220', 18, [0, 1, 3]),
        ('LiveDead_201220', 19, [1, 3])
    ]

    for _tuple in _tuples:
        _experiment, _series_id, _dead_cell_ids = _tuple
        _properties = load.image_properties(_experiment, _series_id)
        _properties['dead_cell_ids'] = _dead_cell_ids
        save.image_properties(_experiment, _series_id, _properties)
        print(_tuple)


if __name__ == '__main__':
    main()
