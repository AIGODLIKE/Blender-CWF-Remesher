import platform
if platform.system() == "Darwin":
    from .macos import remesh
elif platform.system() == "Windows":
    from .windows import remesh
else:
    from .linux import remesh