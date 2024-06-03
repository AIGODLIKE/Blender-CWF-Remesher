import bpy
from .operators import RemeshOperator, ExportXyz, ImportXyz


class TestPanel(bpy.types.Panel):
    bl_idname = "OBJECT_PT_CWF"
    bl_label = "CWF"
    bl_space_type = "VIEW_3D"
    bl_region_type = "UI"

    def draw(self, context):
        layout = self.layout
        layout.prop(context.scene.cwf_prop, "samples")
        layout.operator(RemeshOperator.bl_idname)
        # layout.prop(context.scene.cwf_prop, "xyz_path")
        # layout.operator(ImportXyz.bl_idname)
        # layout.operator(ExportXyz.bl_idname)
