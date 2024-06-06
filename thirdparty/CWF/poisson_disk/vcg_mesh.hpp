/**
    @file vcg_mesh.hpp
    @brief vcg库mesh的处理和定义

    vcg库是用于网格处理（网格修复、孔填充、网格平滑、网格采样等）的强大工具.

    参考链接 http://vcg.sourceforge.net/index.php/Tutorial

    该头文件定义了 meshes/vertices/faces/...等可以被vcg处理的类.

    例子(处理vcg网格):
    @code
    #include "vcg_mesh.hpp"
    #include <vcg/complex/algorithms/smooth.h>
    {
        MyMesh vcg_mesh;
        // 填充网格 ...
        // 执行经典10次 laplacian 平滑迭代
        vcg::tri::Smooth<MyMesh>::VertexCoordLaplacian(vcg_mesh, 10);
    }
    @endcode
*/

#pragma once

#include <cassert>
#include <vector>

#include <vcg/simplex/face/base.h>
#include <vcg/simplex/face/component.h>
#include <vcg/simplex/vertex/base.h>
#include <vcg/simplex/vertex/component.h>

#include <vcg/complex/allocate.h>
#include <vcg/complex/complex.h>

#include <vcg/complex/algorithms/point_sampling.h>

#include <Utils/Common.hpp>

namespace vcg {
class MyVertex;
class MyFace;
class MyMesh;
class MySampler;

using SurfaceSampling = vcg::tri::SurfaceSampling<MyMesh, MySampler>;

class MyUsedTypes
    : public vcg::UsedTypes<vcg::Use<MyVertex>::AsVertexType, vcg::Use<MyFace>::AsFaceType> {};

/** vcg顶点定义
    使用案例:
    @code
    MyVertex v;
    // 访问位置
    vcg::vertex::Coord3f pos = v.P();
    pos.X() = pos.Y() = pos.Z() = 0.f;
    // 访问法线
    vcg::vertex::Normal3f nor = v.N();
    nor.X() = nor.Y() = nor.Z() = 0.f;
    @endcode
*/
class MyVertex : public vcg::Vertex<MyUsedTypes, vcg::vertex::Coord3f, vcg::vertex::Normal3f> {};

/* vcg面定义 */
class MyFace : public vcg::Face<MyUsedTypes, vcg::face::VertexRef> {};

/* vcg-mesh 类型 */
class MyMesh : public vcg::tri::TriMesh<std::vector<MyVertex>, std::vector<MyFace>> {
public:
  using NormalType = vcg::Point3f;

  MyMesh();

  /// 拷贝构造
  MyMesh(const MyMesh &m);

  /// 将三角形拼接到网格.
  /// @param verts: 顶点列表 verts[i*3+0~3] => xyz etc.
  /// @param triangles: 三角形顶点索引列表(3个一组).
  /// @param n_vert: 顶点数量. 'verts' array size == n_vert * 3
  /// @param n_face: 三角形数量. 'triangles' array size == n_face * 3
  /// @warning you might want to update the bounding box
  void concat(const float *verts, const int *triangles, int n_vert, int n_face);

  void concat(const Eigen::MatrixXf &V, const Eigen::MatrixXi &F);

  void fill_mesh(const float *verts, const int *triangles, int n_vert, int n_face);

  void fill_mesh(const Eigen::MatrixXf &V, const Eigen::MatrixXi &F);

  /// 设置法线
  /// @param normals: normals大小 == this->vert.size() * 3
  void set_normals(const float *normals);

  void set_normals(const Eigen::MatrixXf &N);

  /// 更新bounding box
  void update_bb();

  /// 使用vcg mesh的顶点填充 vertices
  template <class Vert_t> void get_vertices(std::vector<Vert_t> &vertices) {
    assert(sizeof(Vert_t) == (sizeof(float) * 3));

    const int n_vert = this->vert.size();
    vertices.clear();
    vertices.resize(n_vert);
    float *ptr = (float *)&(vertices[0]);
    MyMesh::ConstVertexIterator vi = this->vert.begin();
    for (int i = 0; i < n_vert; ++i, ++vi) {
      MyMesh::CoordType coord = (*vi).P();
      ptr[i * 3 + 0] = coord.X();
      ptr[i * 3 + 1] = coord.Y();
      ptr[i * 3 + 2] = coord.Z();
    }
  }

  int nb_vert() const { return this->vn; }
};

/// 存储采样到std::vector<MyMesh::CoordType>
class VecSampler : public vcg::tri::TrivialSampler<MyMesh> {};

/// 存储采样到vcg_mesh
class MySampler {
public:
  MyMesh *mesh;

  MySampler(MyMesh *_m) : mesh(_m) {}

  void reset() { mesh->Clear(); }

  /// 新增顶点
  void AddVert(const MyMesh::VertexType &p) {
    vcg::tri::Allocator<MyMesh>::AddVertices(*mesh, 1);
    mesh->vert.back().ImportData(p);
  }

  /// 新增面
  void AddFace(const MyMesh::FaceType &f, MyMesh::CoordType p) {
    vcg::tri::Allocator<MyMesh>::AddVertices(*mesh, 1);
    // mesh->vert.back().P() = f.P(0) * p[0] + f.P(1) * p[1] + f.P(2) * p[2];
    // mesh->vert.back().N() = f.V(0)->N() * p[0] + f.V(1)->N() * p[1] + f.V(2)->N() * p[2];

    // ref: meshLab
    mesh->vert.back().P() = f.cP(0) * p[0] + f.cP(1) * p[1] + f.cP(2) * p[2];
    mesh->vert.back().N() = f.cV(0)->N() * p[0] + f.cV(1)->N() * p[1] + f.cV(2)->N() * p[2];
    // Quality Sampling
    // mesh->vert.back().Q() = f.cV(0)->Q()*p[0] + f.cV(1)->Q()*p[1] + f.cV(2)->Q()*p[2];
  }

  /// 使用重心坐标p从面f检索样本的颜色，并将该颜色写入位置<tp[0]，texHeight-tp[1]>的纹理图像中.
  /// 如果 edgeDist > 0，则即使在面区域之外（在纹理空间中），对应点也会影响面颜色.
  /// 示例: vcg/complex/algorithms/point_sampling.h/: class TrivialSampler
  void AddTextureSample(const MyMesh::FaceType &,
                        const MyMesh::CoordType &,
                        const vcg::Point2i &,
                        float) {}
};

namespace MyAlgorithms {

struct PoisonParams {
  float _radius = 0.f;   // 最小 radius
  int _n_samples = 1000; // 样本最大数量
  int _montecarlo_rate = 20;
  bool _sub_sampling = false; // 降采样
  bool _approx_geodesic_dist = false;
  bool _refine = false; // 使用现有顶点进行采样
};

void poison_disk_sampling(MyMesh &vcg_mesh, const PoisonParams &params, MyMesh &samples);

} // namespace MyAlgorithms

} // namespace vcg
