import platform


USER = 'Assaf'
# USER = 'Yoni'


def working_directory():
    if USER == 'Assaf':
        if platform.system() == 'Linux':
            return '/data/work/assafna/'
        else:
            return 'G:\\My Drive\\BGU\\Thesis\\Cell-ECM & Cell-ECM-Cell Project\\Data'
    elif USER == 'Yoni':
        if platform.system() == 'Linux':
            return '/data/work/assafna/'
        else:
            return None
