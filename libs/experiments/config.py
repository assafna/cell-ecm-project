SINGLE_CELL = [
    'Single_Cell_Ortal'
]

CELL_PAIRS = [
    'SN16',
    'SN41',
    'SN44',
    'SN45',
    'SN20_Bleb_fromStart',
    'SN26_BlebAdded',
    'Bleb_150',
    'LiveDead_201115',
    'LiveDead_201117',
    'LiveLive_201207',
    'DeadDead_201208'
]

DEAD_LIVE = [
    'LiveDead_201115',
    'LiveDead_201117',
    'DeadDead_201208'
]

BLEB = [
    'SN20_Bleb_fromStart',
    'SN26_BlebAdded',
    'Bleb_150'
]

BLEB_FROM_START = [
    'SN20_Bleb_fromStart',
    'Bleb_150'
]

AFTER_BLEB_INJECTION_FIRST_TIME_FRAME = {
    'SN20_Bleb_fromStart': 0,
    'SN26_BlebAdded': 9,
    'Bleb_150': 0
}

MINIMUM_TIME_FRAMES_CORRELATION = {
    'regular_temporal_resolution': 15,
    'high_temporal_resolution': 50
}

DENSITY_TIME_FRAME = {
    'regular_temporal_resolution': 18,
    'high_temporal_resolution': 52
}

HIGH_TEMPORAL_RESOLUTION_IN_MINUTES = 5

DERIVATIVE = 1

OUT_OF_BOUNDARIES = False
MAX_FRACTION_OUT_OF_BOUNDARIES_BLACK_PIXELS = 0.05

IMAGE_FIBER_CHANNEL_INDEX = 0
IMAGE_CELL_CHANNEL_INDEX = 1

AVERAGE_CELL_DIAMETER_IN_MICRONS = 15
QUANTIFICATION_WINDOW_START_BY_AVERAGE_CELL_DIAMETER = True

QUANTIFICATION_WINDOW_LENGTH_IN_CELL_DIAMETER = 1
QUANTIFICATION_WINDOW_WIDTH_IN_CELL_DIAMETER = 1
QUANTIFICATION_WINDOW_HEIGHT_IN_CELL_DIAMETER = 1


def all_experiments():
    return SINGLE_CELL + CELL_PAIRS


def colors(_n):
    if _n == 1:
        return '#ea8500'
    elif _n == 2:
        return ['#844b00', '#ea8500']
    elif _n == 3:
        return ['#844b00', '#ea8500', '#edbc80']
    else:
        raise Exception('No colors encoded')
