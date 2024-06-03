from pathlib import Path
import importlib
import traceback


def reg_ex(action="register"):
    package_dir = Path(__file__).parent
    for module_path in package_dir.iterdir():
        if not module_path.is_dir() or not module_path.joinpath('__init__.py').exists():
            continue

        module_name = module_path.name
        module = importlib.import_module(f'.{module_path.name}', package=__package__)
        if not hasattr(module, action):
            continue
        try:
            getattr(module, action)()
            print(f'Module {module_name} {action}ed.')
        except Exception:
            traceback.print_exc()


def register():
    """
    遍历当前目录下的所有模块，并尝试调用其 register 方法。
    """
    reg_ex("register")


def unregister():
    """
    遍历当前目录下的所有模块，并尝试调用其 unregister 方法。
    """
    reg_ex("unregister")

