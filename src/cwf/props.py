import bpy


class Props(bpy.types.PropertyGroup):
    xyz_path: bpy.props.StringProperty(name="XYZ Path", subtype="FILE_PATH")
    samples: bpy.props.IntProperty(name="Samples", default=1000, min=10, max=100000)