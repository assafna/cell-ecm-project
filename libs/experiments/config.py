SINGLE_CELL = [
    'Single_Cell_Ortal'
]

PAIRS = [
    'SN16',
    'SN18',
    'SN41',
    'SN44',
    'SN45',
    'SN20_Bleb_fromStart',
    'SN26_BlebAdded'
]

FIBERS_CHANNEL_INDEX = 0
CELLS_CHANNEL_INDEX = 1
MAX_DISTANCE_CHANGE = 20
AVERAGE_CELL_DIAMETER_IN_MICRONS = 15
ROI_START_BY_AVERAGE_CELL_DIAMETER = True

ROI_LENGTH = 1
ROI_WIDTH = 1
ROI_HEIGHT = 1

NO_RETURN = False


def experiments():
    return SINGLE_CELL + PAIRS
