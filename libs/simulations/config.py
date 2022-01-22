CELL_DIAMETER = 0.08

CELLS_DIAMETERS_FILE_NAME = 'cells_diameters.txt'
ELEMENTS_FILE_NAME = 'elements.txt'

NET_DIMENSIONS = {
    'x1': -2,
    'y1': -2,
    'x2': 2,
    'y2': 2
}

ORIGIN_COORDINATES = {
    'x': 0,
    'y': 0
}

QUANTIFICATION_WINDOW_WIDTH_IN_CELL_DIAMETER = 1
QUANTIFICATION_WINDOW_HEIGHT_IN_CELL_DIAMETER = 1
EPSILON = 1.0e-8


def colors(_n):
    if _n == 1:
        return '#011f4b'
    elif _n == 2:
        return ['#011f4b', '#005b96']
    elif _n == 3:
        return ['#011f4b', '#005b96', '#74c2e8']
    elif _n == 4:
        return ['#011f4b', '#00417c', '#2e82bf', '#56caed']
    elif _n == 5:
        return ['#011f4b', '#013b66', '#1c6ea0', '#3ab1cc', '#73d8ef']
    else:
        raise Exception('No colors encoded')
