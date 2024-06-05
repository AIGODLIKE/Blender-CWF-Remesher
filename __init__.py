bl_info = {
    'name': 'CWF Remesh',
    'author': '会飞的键盘侠',
    'version': (0, 0, 1),
    'blender': (3, 0, 0),
    'location': '3DView->Panel',
    'category': 'Mesh',
    'doc_url': "https://bing.com"
}

import sys
import bpy
from .src import register as reg
from .src import unregister as unreg
from .utils.logger import logger


def register():
    logger.debug(f'{bl_info["name"]}: register')
    reg()


def unregister():
    logger.debug(f'{bl_info["name"]}: unregister')
    unreg()
    modules_update()


def modules_update():
    from .utils.logger import logger
    logger.close()
    modules = []
    for i in sys.modules:
        if i.startswith(__package__) and i != __package__:
            modules.append(i)
    for i in modules:
        del sys.modules[i]
    del sys.modules[__package__]
