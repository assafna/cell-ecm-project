import os

import paths


def simulations():
    return [_simulation for _simulation in os.listdir(paths.RAW)
            if os.path.isdir(os.path.join(paths.RAW, _simulation))]