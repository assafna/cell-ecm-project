from libs.experiments import load, save


def main():
    _tuples = [
        ('LiveDead_201117', 1, [1]),
        ('LiveDead_201117', 2, [0]),
        ('LiveDead_201117', 3, [0, 1]),
        ('LiveDead_201117', 4, [0, 1]),
        ('LiveDead_201117', 5, [0]),
        ('LiveDead_201117', 6, [0]),
        ('LiveDead_201117', 7, []),
        ('LiveDead_201117', 8, [0]),
        ('LiveDead_201117', 9, [0, 3]),
        ('LiveDead_201117', 10, [0]),
        ('LiveDead_201117', 11, [1]),
        ('LiveDead_201117', 12, []),
        ('LiveDead_201117', 13, [1, 2]),

        ('LiveDead_201117', 15, [1, 2]),
        ('LiveDead_201117', 16, [0]),
        ('LiveDead_201117', 17, [1, 2, 3]),
        ('LiveDead_201117', 18, [0])
    ]

    for _tuple in _tuples:
        _experiment, _series_id, _dead_cell_ids = _tuple
        _properties = load.image_properties(_experiment, _series_id)
        _properties['dead_cell_ids'] = _dead_cell_ids
        save.image_properties(_experiment, _series_id, _properties)
        print(_tuple)


if __name__ == '__main__':
    main()
