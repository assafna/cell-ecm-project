import os
import platform


def working_directory():
    if platform.system() == 'Linux':
        return '/data/work/assafna/'
    else:
        return 'G:\\My Drive\\BGU\\Thesis\\Cell-ECM & Cell-ECM-Cell Project\\Data'


PLOTS = os.path.join(working_directory(), 'Plots')
