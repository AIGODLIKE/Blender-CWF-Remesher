import bpy
from .operators import RemeshOperator, ExportXyz, ImportXyz
from .panels import TestPanel
from .props import Props


def register():
    bpy.utils.register_class(RemeshOperator)
    bpy.utils.register_class(ExportXyz)
    bpy.utils.register_class(ImportXyz)
    bpy.utils.register_class(TestPanel)
    bpy.utils.register_class(Props)
    bpy.types.Scene.cwf_prop = bpy.props.PointerProperty(type=Props)


def unregister():
    del bpy.types.Scene.cwf_prop
    bpy.utils.unregister_class(Props)
    bpy.utils.unregister_class(TestPanel)
    bpy.utils.unregister_class(ImportXyz)
    bpy.utils.unregister_class(ExportXyz)
    bpy.utils.unregister_class(RemeshOperator)

