import fire
from hkvsobekpy.io.bui import __bui_class
from hkvsobekpy.io.his import __his_class
from hkvsobekpy.core.waterlevelstat import __waterlevelstat_class
from hkvsobekpy.core.plausibility import __plausibility_class

__doc__ = """package for water-statistics and plausibility checker using his- and bui-file"""
__version__ = '1.3.0'

# initiate class for .bui-files
__bui = __bui_class()
read_bui = __bui.read_bui
write_bui = __bui.write_bui

# initiate class for .his-files
read_his = __his_class()

# initiate class for waterlevelstats
waterlevelstat = __waterlevelstat_class()

# initiate class for plausibilitychecker
plausibility = __plausibility_class()

if __name__ == '__main__':
    fire.Fire()