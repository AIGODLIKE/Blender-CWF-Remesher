from typing import Set
from math import pi
import bpy

from pathlib import Path
from bpy.types import Context, Event
from ...utils.logger import logger

DESKTOP = Path().home().joinpath("Desktop")


def write_xyz(path: str, verts: list[bpy.types.MeshVertex]):
    with open(path, "w+", encoding="utf-8") as f:
        for v in verts:
            f.write(f"{v.co.x} {v.co.y} {v.co.z} {v.normal.x} {v.normal.y} {v.normal.z}\n")
        f.flush()


class RemeshOperator(bpy.types.Operator):
    bl_idname = "mesh.cwfremesh"
    bl_description = "CWF Remesh"
    bl_label = "CWF Remesh"
    bl_translation_context = "ops_ctx"

    @classmethod
    def poll(cls, context):
        return bpy.context.object and bpy.context.object.type == "MESH"

    def invoke(self, context: Context, event: Event) -> Set[int] | Set[str]:
        return self.execute(context)

    def execute(self, context):
        # 判断确实选择了网格
        obj = bpy.context.object
        if not obj.data.polygons:
            logger.info("No polygons found")
            return {"CANCELLED"}
        bpy.ops.object.select_all(action="DESELECT")
        obj.select_set(True)
        bpy.context.view_layer.objects.active = obj

        if not DESKTOP.joinpath("cwf_test").exists():
            DESKTOP.joinpath("cwf_test").mkdir()
        if not DESKTOP.joinpath("cwf_test/_LBFGSOUT").exists():
            DESKTOP.joinpath("cwf_test/_LBFGSOUT").mkdir()
        oname = obj.name
        objpath = DESKTOP.joinpath(f"cwf_test/{oname}.obj")
        bpy.ops.wm.obj_export(filepath=objpath.as_posix(),
                              export_selected_objects=True,
                              apply_modifiers=True,
                              export_triangulated_mesh=True,
                              export_uv=False,
                              export_materials=False,
                              )
        from ...cxxlibs import remesh
        p = objpath.as_posix()
        n = context.scene.cwf_prop.samples
        try:
            remesh.remesh(p, n)
            # 结束后得到 Remesh.obj 导入blender
            objpath2 = DESKTOP.joinpath(f"cwf_test/_LBFGSOUT/Ours_{n}_{oname}_Remesh.obj")
            self.import_obj(objpath2.as_posix())
            # 清理
            for file in DESKTOP.joinpath("cwf_test/_LBFGSOUT").iterdir():
                file.unlink()
        except RuntimeError:
            logger.error("Remesh failed")
            return {"CANCELLED"}
        return {"FINISHED"}

    def import_obj(self, filepath) -> list[bpy.types.Object]:
        rec_objs = set(bpy.context.scene.objects)
        bpy.ops.wm.obj_import(filepath=filepath)
        new_objs = set(bpy.context.scene.objects) - rec_objs
        for obj in new_objs:
            bpy.context.view_layer.objects.active = obj
            # bpy.ops.object.shade_smooth()
        return list(new_objs)

    def prepare_nodegroup(self) -> bpy.types.NodeGroup:
        if "._CWF Remesh" in bpy.data.node_groups:
            return bpy.data.node_groups["._CWF Remesh"]
        node_group = bpy.data.node_groups.new(name="._CWF Remesh", type="GeometryNodeTree")
        ndist = node_group.nodes.new(type="GeometryNodeDistributePointsOnFaces")
        ndist.location = 200, 0
        ndist.distribute_method = "POISSON"
        ndist.inputs["Distance Min"].default_value = 0.01
        ndist.inputs["Density Max"].default_value = 400
        ndist.inputs["Density Factor"].default_value = 1
        ntov = node_group.nodes.new(type="GeometryNodePointsToVertices")
        ntov.location = 400, 0

        ninput = node_group.nodes.new(type="NodeGroupInput")

        node_group.interface.new_socket("Geometry", socket_type="NodeSocketGeometry", in_out="INPUT")
        node_group.interface.new_socket("Geometry", socket_type="NodeSocketGeometry", in_out="OUTPUT")

        noutput = node_group.nodes.new(type="NodeGroupOutput")
        noutput.location = 600, 0

        node_group.links.new(ndist.inputs[0], ninput.outputs[0])
        node_group.links.new(ntov.inputs[0], ndist.outputs[0])
        node_group.links.new(noutput.inputs[0], ntov.outputs[0])
        return node_group


class ExportXyz(bpy.types.Operator):
    bl_idname = "mesh.export_xyz"
    bl_description = "Export XYZ"
    bl_label = "Export XYZ"
    bl_translation_context = "ops_ctx"

    def invoke(self, context: Context, event: Event) -> Set[int] | Set[str]:
        return self.execute(context)

    def execute(self, context):
        obj = bpy.context.object
        if not obj or not obj.data:
            return {"CANCELLED"}
        if obj.type != "MESH":
            return {"CANCELLED"}

        bpy.context.view_layer.update()
        eobj = obj.evaluated_get(context.view_layer.depsgraph)
        # 点云xyz
        xyzpath = context.scene.cwf_prop.xyz_path
        write_xyz(xyzpath, eobj.data.vertices)
        return {"FINISHED"}


class ImportXyz(bpy.types.Operator):
    bl_idname = "import_mesh.xyz"
    bl_label = "Import XYZ"
    bl_description = "Import XYZ"
    bl_options = {"REGISTER", "UNDO"}

    def execute(self, context):
        p = context.scene.cwf_prop.xyz_path
        if not Path(p).exists():
            return {"CANCELLED"}
        me: bpy.types.Mesh = bpy.data.meshes.new("test")
        me.from_pydata(vertices=self.read_xyz(p), edges=[], faces=[])
        obj = bpy.data.objects.new("xyz import", object_data=me)
        bpy.context.scene.collection.objects.link(obj)
        obj.rotation_euler = (pi * 0.5, 0, 0)
        bpy.context.view_layer.objects.active = obj
        obj.select_set(True)
        return {"FINISHED"}

    def read_xyz(self, path) -> list[tuple[float, float, float]]:
        vertices = []
        with open(path, "r", encoding="utf-8") as f:
            for line in f:
                if line.strip() == "":
                    continue
                v = line.strip().split(" ")
                vertices.append((float(v[0]), float(v[1]), float(v[2])))
        return vertices
