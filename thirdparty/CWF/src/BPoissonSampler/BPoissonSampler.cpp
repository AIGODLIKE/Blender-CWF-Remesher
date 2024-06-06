#include <array>
#include <span>
#include <stdint.h>
#include <vector>

#define STREQ(a, b) (strcmp(a, b) == 0)

using float3 = std::array<float, 3>;
using int3 = std::array<int, 3>;

/** #CustomData.type */
typedef enum eCustomDataType {
  /* Used by GLSL attributes in the cases when we need a delayed CD type
   * assignment (in the cases when we don't know in advance which layer
   * we are addressing).
   */
  CD_AUTO_FROM_NAME = -1,

#ifdef DNA_DEPRECATED_ALLOW
  CD_MVERT = 0,
  CD_MSTICKY = 1,
#endif
  CD_MDEFORMVERT = 2, /* Array of `MDeformVert`. */
#ifdef DNA_DEPRECATED_ALLOW
  CD_MEDGE = 3,
#endif
  CD_MFACE = 4,
  CD_MTFACE = 5,
  CD_MCOL = 6,
  CD_ORIGINDEX = 7,
  /**
   * Used for derived face corner normals on mesh `ldata`, since currently they are not computed
   * lazily. Derived vertex and polygon normals are stored in #Mesh_Runtime.
   */
  CD_NORMAL = 8,
#ifdef DNA_DEPRECATED_ALLOW
  CD_FACEMAP = 9,
#endif
  CD_PROP_FLOAT = 10,
  CD_PROP_INT32 = 11,
  CD_PROP_STRING = 12,
  CD_ORIGSPACE = 13, /* for modifier stack face location mapping */
  CD_ORCO = 14,      /* undeformed vertex coordinates, normalized to 0..1 range */
#ifdef DNA_DEPRECATED_ALLOW
  CD_MTEXPOLY = 15,
  CD_MLOOPUV = 16,
#endif
  CD_PROP_BYTE_COLOR = 17,
  CD_TANGENT = 18,
  CD_MDISPS = 19,
  CD_PREVIEW_MCOL = 20,           /* For displaying weight-paint colors. */
                                  /*  CD_ID_MCOL          = 21, */
  /* CD_TEXTURE_MLOOPCOL = 22, */ /* UNUSED */
  CD_CLOTH_ORCO = 23,
/* CD_RECAST = 24, */ /* UNUSED */

#ifdef DNA_DEPRECATED_ALLOW
  CD_MPOLY = 25,
  CD_MLOOP = 26,
#endif
  CD_SHAPE_KEYINDEX = 27,
  CD_SHAPEKEY = 28,
#ifdef DNA_DEPRECATED_ALLOW
  CD_BWEIGHT = 29,
  CD_CREASE = 30,
#endif
  CD_ORIGSPACE_MLOOP = 31,
  CD_PREVIEW_MLOOPCOL = 32,
  CD_BM_ELEM_PYPTR = 33,

  CD_PAINT_MASK = 34,
  CD_GRID_PAINT_MASK = 35,
  CD_MVERT_SKIN = 36,
  CD_FREESTYLE_EDGE = 37,
  CD_FREESTYLE_FACE = 38,
  CD_MLOOPTANGENT = 39,
  CD_TESSLOOPNORMAL = 40,
  CD_CUSTOMLOOPNORMAL = 41,
#ifdef DNA_DEPRECATED_ALLOW
  CD_SCULPT_FACE_SETS = 42,
#endif

  /* CD_LOCATION = 43, */ /* UNUSED */
  /* CD_RADIUS = 44, */   /* UNUSED */
  CD_PROP_INT8 = 45,
  /* Two 32-bit signed integers. */
  CD_PROP_INT32_2D = 46,

  CD_PROP_COLOR = 47,
  CD_PROP_FLOAT3 = 48,
  CD_PROP_FLOAT2 = 49,
  CD_PROP_BOOL = 50,

  /* CD_HAIRLENGTH = 51, */ /* UNUSED */

  CD_PROP_QUATERNION = 52,

  CD_NUMTYPES = 53,
} eCustomDataType;

/** Descriptor and storage for a custom data layer. */
typedef struct CustomDataLayer {
  /** Type of data in layer. */
  int type;
  /** In editmode, offset of layer in block. */
  int offset;
  /** General purpose flag. */
  int flag;
  /** Number of the active layer of this type. */
  int active;
  /** Number of the layer to render. */
  int active_rnd;
  /** Number of the layer to render. */
  int active_clone;
  /** Number of the layer to render. */
  int active_mask;
  /** Shape key-block unique id reference. */
  int uid;
  /** Layer name, MAX_CUSTOMDATA_LAYER_NAME. */
  char name[68];
  char _pad1[4];
  /** Layer data. */
  void *data;
  /**
   * Run-time identifier for this layer. Can be used to retrieve information about where this
   * attribute was created.
   */
  const void *anonymous_id;
  /**
   * Run-time data that allows sharing `data` with other entities (mostly custom data layers on
   * other geometries).
   */
  const void *sharing_info;
} CustomDataLayer;

typedef struct CustomDataExternal {
  /** FILE_MAX. */
  char filepath[1024];
} CustomDataExternal;

typedef struct CustomData {
  /** CustomDataLayers, ordered by type. */
  CustomDataLayer *layers;
  /**
   * runtime only! - maps types to indices of first layer of that type,
   * MUST be >= CD_NUMTYPES, but we can't use a define here.
   * Correct size is ensured in CustomData_update_typemap assert().
   */
  int typemap[53];
  /** Number of layers, size of layers array. */
  int totlayer, maxlayer;
  /** In editmode, total size of all data layers. */
  int totsize;
  /** (BMesh Only): Memory pool for allocation of blocks. */
  struct BLI_mempool *pool;
  /** External file storing custom-data layers. */
  CustomDataExternal *external;

  int CustomData_get_named_layer_index(const eCustomDataType type, const char *name) {
    for (int i = 0; i < totlayer; i++) {
      if (layers[i].type == type) {
        if (STREQ(layers[i].name, name)) {
          return i;
        }
      }
    }
    return -1;
  }

  const void* CustomData_get_layer_named(const eCustomDataType type, const char *name) {
    int layer_index = CustomData_get_named_layer_index(type, name);
    if (layer_index == -1) {
      return nullptr;
    }
    return layers[layer_index].data;
  }
} CustomData;

typedef struct ID_Runtime_Remap {
  /** Status during ID remapping. */
  int status;
  /** During ID remapping the number of skipped use cases that refcount the data-block. */
  int skipped_refcounted;
  /**
   * During ID remapping the number of direct use cases that could be remapped
   * (e.g. obdata when in edit mode).
   */
  int skipped_direct;
  /** During ID remapping, the number of indirect use cases that could not be remapped. */
  int skipped_indirect;
} ID_Runtime_Remap;

typedef struct ID_Runtime {
  ID_Runtime_Remap remap;
} ID_Runtime;

typedef struct ID {
  void *next, *prev;
  struct ID *newid;

  struct Library *lib;

  /** If the ID is an asset, this pointer is set. Owning pointer. */
  struct AssetMetaData *asset_data;

  /** MAX_ID_NAME. */
  char name[66];
  /**
   * LIB_... flags report on status of the data-block this ID belongs to
   * (persistent, saved to and read from .blend).
   */
  short flag;
  /**
   * LIB_TAG_... tags (runtime only, cleared at read time).
   */
  int tag;
  int us;
  int icon_id;
  unsigned int recalc;
  /**
   * Used by undo code. recalc_after_undo_push contains the changes between the
   * last undo push and the current state. This is accumulated as IDs are tagged
   * for update in the depsgraph, and only cleared on undo push.
   *
   * recalc_up_to_undo_push is saved to undo memory, and is the value of
   * recalc_after_undo_push at the time of the undo push. This means it can be
   * used to find the changes between undo states.
   */
  unsigned int recalc_up_to_undo_push;
  unsigned int recalc_after_undo_push;

  /**
   * A session-wide unique identifier for a given ID, that remain the same across potential
   * re-allocations (e.g. due to undo/redo steps).
   */
  unsigned int session_uuid;

  void *properties;

  /** Reference linked ID which this one overrides. */
  void *override_library;

  /**
   * Only set for data-blocks which are coming from copy-on-write, points to
   * the original version of it.
   * Also used temporarily during memfile undo to keep a reference to old ID when found.
   */
  struct ID *orig_id;

  /**
   * Holds the #PyObject reference to the ID (initialized on demand).
   *
   * This isn't essential, it could be removed however it gives some advantages:
   *
   * - Every time the #ID is accessed a #BPy_StructRNA doesn't have to be created & destroyed
   *   (consider all the polling and drawing functions that access ID's).
   *
   * - When this #ID is deleted, the #BPy_StructRNA can be invalidated
   *   so accessing it from Python raises an exception instead of crashing.
   *
   *   This is of limited benefit though, as it doesn't apply to non #ID data
   *   that references this ID (the bones of an armature or the modifiers of an object for e.g.).
   */
  void *py_instance;

  /**
   * Weak reference to an ID in a given library file, used to allow re-using already appended data
   * in some cases, instead of appending it again.
   *
   * May be NULL.
   */
  struct LibraryWeakReference *library_weak_reference;

  struct ID_Runtime runtime;
} ID;

typedef struct Mesh {
  ID id;
  /** Animation data (must be immediately after id for utilities to use it). */
  struct AnimData *adt;

  /** Old animation system, deprecated for 2.5. */
  void *ipo;
  struct Key *key;

  /**
   * An array of materials, with length #totcol. These can be overridden by material slots
   * on #Object. Indices in the "material_index" attribute control which material is used for every
   * face.
   */
  struct Material **mat;

  /** The number of vertices in the mesh, and the size of #vert_data. */
  int verts_num;
  /** The number of edges in the mesh, and the size of #edge_data. */
  int edges_num;
  /** The number of faces in the mesh, and the size of #face_data. */
  int faces_num;
  /** The number of face corners in the mesh, and the size of #corner_data. */
  int corners_num;

  /**
   * Array owned by mesh. See #Mesh::faces() and #OffsetIndices.
   *
   * This array is shared based on the bke::MeshRuntime::face_offsets_sharing_info.
   * Avoid accessing directly when possible.
   */
  int *face_offset_indices;

  CustomData vert_data;
  CustomData edge_data;
  CustomData face_data;
  CustomData corner_data;


  /**
   * List of vertex group (#bDeformGroup) names and flags only. Actual weights are stored in dvert.
   * \note This pointer is for convenient access to the #CD_MDEFORMVERT layer in #vert_data.
   */
  ListBase vertex_group_names;
  /** The active index in the #vertex_group_names list. */
  int vertex_group_active_index;

  /**
   * The index of the active attribute in the UI. The attribute list is a combination of the
   * generic type attributes from vertex, edge, face, and corner custom data.
   */
  int attributes_active_index;

  /**
   * Runtime storage of the edit mode mesh. If it exists, it generally has the most up-to-date
   * information about the mesh.
   * \note When the object is available, the preferred access method is #BKE_editmesh_from_object.
   */
  struct BMEditMesh *edit_mesh;

  /**
   * This array represents the selection order when the user manually picks elements in edit-mode,
   * some tools take advantage of this information. All elements in this array are expected to be
   * selected, see #BKE_mesh_mselect_validate which ensures this. For procedurally created meshes,
   * this is generally empty (selections are stored as boolean attributes in the corresponding
   * custom data).
   */
  struct MSelect *mselect;

  /** The length of the #mselect array. */
  int totselect;

  /**
   * In most cases the last selected element (see #mselect) represents the active element.
   * For faces we make an exception and store the active face separately so it can be active
   * even when no faces are selected. This is done to prevent flickering in the material properties
   * and UV Editor which base the content they display on the current material which is controlled
   * by the active face.
   *
   * \note This is mainly stored for use in edit-mode.
   */
  int act_face;

  /**
   * An optional mesh owned elsewhere (by #Main) that can be used to override
   * the texture space #loc and #size.
   * \note Vertex indices should be aligned for this to work usefully.
   */
  struct Mesh *texcomesh;

  /** Texture space location and size, used for procedural coordinates when rendering. */
  float texspace_location[3];
  float texspace_size[3];
  char texspace_flag;

  /** Various flags used when editing the mesh. */
  char editflag;
  /** Mostly more flags used when editing or displaying the mesh. */
  uint16_t flag;

  /**
   * The angle for auto smooth in radians. `M_PI` (180 degrees) causes all edges to be smooth.
   */
  float smoothresh;

  /** Per-mesh settings for voxel remesh. */
  float remesh_voxel_size;
  float remesh_voxel_adaptivity;

  int face_sets_color_seed;
  /* Stores the initial Face Set to be rendered white. This way the overlay can be enabled by
   * default and Face Sets can be used without affecting the color of the mesh. */
  int face_sets_color_default;

  /** The color attribute currently selected in the list and edited by a user. */
  char *active_color_attribute;
  /** The color attribute used by default (i.e. for rendering) if no name is given explicitly. */
  char *default_color_attribute;

  /**
   * User-defined symmetry flag (#eMeshSymmetryType) that causes editing operations to maintain
   * symmetrical geometry. Supported by operations such as transform and weight-painting.
   */
  char symmetry;

  /** Choice between different remesh methods in the UI. */
  char remesh_mode;

  /** The length of the #mat array. */
  short totcol;

  /**
   * Deprecated flag for choosing whether to store specific custom data that was built into #Mesh
   * structs in edit mode. Replaced by separating that data to separate layers. Kept for forward
   * and backwards compatibility.
   */
  char cd_flag DNA_DEPRECATED;
  char subdiv DNA_DEPRECATED;
  char subdivr DNA_DEPRECATED;
  char subsurftype DNA_DEPRECATED;

  /** Deprecated pointer to mesh polygons, kept for forward compatibility. */
  struct MPoly *mpoly DNA_DEPRECATED;
  /** Deprecated pointer to face corners, kept for forward compatibility. */
  struct MLoop *mloop DNA_DEPRECATED;

  /** Deprecated array of mesh vertices, kept for reading old files, now stored in #CustomData. */
  struct MVert *mvert DNA_DEPRECATED;
  /** Deprecated array of mesh edges, kept for reading old files, now stored in #CustomData. */
  struct MEdge *medge DNA_DEPRECATED;
  /** Deprecated "Vertex group" data. Kept for reading old files, now stored in #CustomData. */
  struct MDeformVert *dvert DNA_DEPRECATED;
  /** Deprecated runtime data for tessellation face UVs and texture, kept for reading old files. */
  struct MTFace *mtface DNA_DEPRECATED;
  /** Deprecated, use mtface. */
  struct TFace *tface DNA_DEPRECATED;
  /** Deprecated array of colors for the tessellated faces, kept for reading old files. */
  struct MCol *mcol DNA_DEPRECATED;
  /** Deprecated face storage (quads & triangles only). Kept for reading old files. */
  struct MFace *mface DNA_DEPRECATED;

  /**
   * Deprecated storage of old faces (only triangles or quads).
   *
   * \note This would be marked deprecated, however the particles still use this at run-time
   * for placing particles on the mesh (something which should be eventually upgraded).
   */
  CustomData fdata_legacy;
  /* Deprecated size of #fdata. */
  int totface_legacy;

  char _pad1[4];

  /**
   * Data that isn't saved in files, including caches of derived data, temporary data to improve
   * the editing experience, etc. The struct is created when reading files and can be accessed
   * without null checks, with the exception of some temporary meshes which should allocate and
   * free the data if they are passed to functions that expect run-time data.
   */
  MeshRuntimeHandle *runtime;
  
  inline float3 * vert_positions() {
    return static_cast<float3 *>(vert_data.CustomData_get_layer_named(CD_PROP_FLOAT3, "position"));
  }

  inline int* corner_verts()
  {
    return static_cast<const int *>(loop_data.CustomData_get_layer_named(CD_PROP_INT32, ".corner_vert"));
  }

  blender::Span<MLoopTri> looptris()
  {
      runtime->looptris_cache.ensure([&](blender::Array<MLoopTri> &r_data) {
      const Span<float3> positions = this->vert_positions();
      const blender::OffsetIndices faces = this->faces();
      const Span<int> corner_verts = this->corner_verts();

      r_data.reinitialize(poly_to_tri_count(faces.size(), corner_verts.size()));

      if (BKE_mesh_face_normals_are_dirty(this)) {
        blender::bke::mesh::looptris_calc(positions, faces, corner_verts, r_data);
      }
      else {
        blender::bke::mesh::looptris_calc_with_normals(
            positions, faces, corner_verts, this->face_normals(), r_data);
      }
    });

    return runtime->looptris_cache.data();
  }
} Mesh;

static void sample_mesh_surface(const Mesh &mesh, const float base_density, const std::span<float> density_factors, const int seed, std::vector<float3> &r_positions,
                                std::vector<float3> &r_bary_coords, std::vector<int> &r_tri_indices) {
  auto positions = mesh.vert_positions();
  auto corner_verts = mesh.corner_verts();
  auto corner_tris = static_cast<std::span<int3>>(mesh.corner_data.layers[0].data);

  for (const int tri_i : corner_tris.index_range()) {
    const int3 &tri = corner_tris[tri_i];
    const int v0_loop = tri[0];
    const int v1_loop = tri[1];
    const int v2_loop = tri[2];
    const float3 &v0_pos = positions[corner_verts[v0_loop]];
    const float3 &v1_pos = positions[corner_verts[v1_loop]];
    const float3 &v2_pos = positions[corner_verts[v2_loop]];

    float corner_tri_density_factor = 1.0f;
    if (!density_factors.is_empty()) {
      const float v0_density_factor = std::max(0.0f, density_factors[v0_loop]);
      const float v1_density_factor = std::max(0.0f, density_factors[v1_loop]);
      const float v2_density_factor = std::max(0.0f, density_factors[v2_loop]);
      corner_tri_density_factor = (v0_density_factor + v1_density_factor + v2_density_factor) / 3.0f;
    }
    const float area = area_tri_v3(v0_pos, v1_pos, v2_pos);

    const int corner_tri_seed = noise::hash(tri_i, seed);
    RandomNumberGenerator corner_tri_rng(corner_tri_seed);

    const int point_amount = corner_tri_rng.round_probabilistic(area * base_density * corner_tri_density_factor);

    for (int i = 0; i < point_amount; i++) {
      const float3 bary_coord = corner_tri_rng.get_barycentric_coordinates();
      float3 point_pos;
      interp_v3_v3v3v3(point_pos, v0_pos, v1_pos, v2_pos, bary_coord);
      r_positions.append(point_pos);
      r_bary_coords.append(bary_coord);
      r_tri_indices.append(tri_i);
    }
  }
}

static void distribute_points_poisson_disk(const Mesh &mesh, const float minimum_distance, const float max_density, const Field<float> &density_factor_field, const Field<bool> &selection_field,
                                           const int seed, Vector<float3> &positions, Vector<float3> &bary_coords, Vector<int> &tri_indices) {
  sample_mesh_surface(mesh, max_density, {}, seed, positions, bary_coords, tri_indices);

  Array<bool> elimination_mask(positions.size(), false);
  update_elimination_mask_for_close_points(positions, minimum_distance, elimination_mask);

  const Array<float> density_factors = calc_full_density_factors_with_selection(mesh, density_factor_field, selection_field);

  update_elimination_mask_based_on_density_factors(mesh, density_factors, bary_coords, tri_indices, elimination_mask.as_mutable_span());

  eliminate_points_based_on_mask(elimination_mask.as_span(), positions, bary_coords, tri_indices);
}