import bpy
from .operators import RemeshOperator, ExportXyz, ImportXyz
from ..i18n.ctx import PCTX


class TestPanel(bpy.types.Panel):
    bl_idname = "OBJECT_PT_CWF"
    bl_label = "CWF"
    bl_category = "CWF"
    bl_space_type = "VIEW_3D"
    bl_region_type = "UI"

    def draw(self, context):
        layout = self.layout
        cfg_col = layout.column(align=True)
        cfg_col.prop(context.scene.cwf_prop, "p_samples")
        if context.scene.cwf_prop.pre_simplify:
            cfg_cr = cfg_col.row(align=True)
            cfg_cr.prop(context.scene.cwf_prop, "pre_simplify", toggle=True, icon="MOD_DECIM", text="")
            cfg_cr.prop(context.scene.cwf_prop, "pre_simplify_ratio")
        else:
            cfg_col.prop(context.scene.cwf_prop, "pre_simplify", toggle=True, icon="MOD_DECIM")
        # advanced: bpy.props.BoolProperty(default=False, name="Advanced", translation_context=PCTX)
        # p_samples: bpy.props.IntProperty(name="Samples", default=1000, min=10, max=100000, translation_context=PCTX)
        # p_fnum: bpy.props.IntProperty(name="Fnum", default=4, min=4, max=100, translation_context=PCTX)
        # p_alpha: bpy.props.FloatProperty(name="Alpha", default=1.0, min=0.01, max=1, translation_context=PCTX)
        # p_eplison: bpy.props.FloatProperty(name="Eplison", default=1.0, min=0.0001, max=1, translation_context=PCTX)
        # p_lambda: bpy.props.FloatProperty(name="Lambda", default=1.0, min=0.0001, max=1, translation_context=PCTX)
        # p_decay: bpy.props.FloatProperty(name="Decay", default=0.95, min=0.01, max=1, translation_context=PCTX)
        cfg_col.prop(context.scene.cwf_prop, "advanced", toggle=True, icon="OPTIONS")
        if context.scene.cwf_prop.advanced:
            cfg_col.prop(context.scene.cwf_prop, "p_fnum")
            cfg_col.prop(context.scene.cwf_prop, "p_alpha")
            cfg_col.prop(context.scene.cwf_prop, "p_eplison")
            cfg_col.prop(context.scene.cwf_prop, "p_lambda")
            cfg_col.prop(context.scene.cwf_prop, "p_decay")
        col = layout.column()
        col.enabled = not RemeshOperator.running
        if RemeshOperator.running:
            col.label(text="Remeshing...", icon="INFO")
        else:
            col.operator(RemeshOperator.bl_idname)
        # layout.prop(context.scene.cwf_prop, "xyz_path")
        # layout.operator(ImportXyz.bl_idname)
        # layout.operator(ExportXyz.bl_idname)
