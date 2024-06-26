import bpy
from ..i18n.ctx import PCTX

KWARGS = {"translation_context": PCTX}
if bpy.app.version < (4, 0):
    KWARGS.pop("translation_context")


class Props(bpy.types.PropertyGroup):
    xyz_path: bpy.props.StringProperty(name="XYZ Path", subtype="FILE_PATH", **KWARGS)
    advanced: bpy.props.BoolProperty(default=False, name="Advanced", **KWARGS)
    p_samples: bpy.props.IntProperty(name="Samples", default=1000, min=10, max=100000, **KWARGS)
    p_fnum: bpy.props.IntProperty(name="Fnum", default=4, min=4, max=100, **KWARGS)
    p_alpha: bpy.props.FloatProperty(name="Alpha", default=1.0, min=0.01, max=1, **KWARGS)
    p_eplison: bpy.props.FloatProperty(name="Eplison", default=1.0, min=0.0001, max=1, **KWARGS)
    p_lambda: bpy.props.FloatProperty(name="Lambda", default=1.0, min=0.0001, max=1, **KWARGS)
    p_decay: bpy.props.FloatProperty(name="Decay", default=0.95, min=0.01, max=1, **KWARGS)

    pre_simplify: bpy.props.BoolProperty(default=True, name="Pre Simplify", **KWARGS)
    pre_simplify_ratio: bpy.props.FloatProperty(default=.3, min=0.01, max=1, name="Pre Simplify Ratio", **KWARGS)

    def set_params(self, params):
        """
        "RemeshParams"
            "samples", samples
            "fnum", fnum
            "alpha", alpha
            "eplison", eplison
            "lambda", lambda
            "decay", decay
        """
        params._samples = self.p_samples
        params._fnum = self.p_fnum
        params._alpha = self.p_alpha
        params._eplison = self.p_eplison
        params._lambda = self.p_lambda
        params._decay = self.p_decay
