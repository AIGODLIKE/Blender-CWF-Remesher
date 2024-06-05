# formt: (t, translate_t, ctx) if not ctx => default "*"
from ..ctx import PCTX, OCTX

translations = (
    ("Samples", "采样", PCTX),
    ("Advanced", "高级设置", PCTX),
    ("XYZ Path", "XYZ路径", PCTX),
    ("Pre Simplify", "预精简", PCTX),
    ("Pre Simplify Ratio", "系数", PCTX),
    ("CWF Remesh", "CWF 网格重建", OCTX),
    ("Export XYZ", "导出XYZ", OCTX),
    ("Import XYZ", "导入XYZ", OCTX),
)
